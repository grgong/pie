"""GFF3/GTF parsing with longest-isoform selection."""

import hashlib
import logging
import os
from collections import defaultdict
from dataclasses import dataclass

import gffutils

log = logging.getLogger(__name__)


class NoGenesFoundError(ValueError):
    """Raised when annotation parsing yields zero usable genes."""


@dataclass
class GeneModel:
    gene_id: str
    transcript_id: str
    chrom: str
    start: int  # 0-based
    end: int  # 0-based half-open
    strand: str
    cds_exons: list[tuple[str, int, int]]  # [(chrom, start, end)] 0-based half-open

    @property
    def cds_length(self) -> int:
        return sum(end - start for _, start, end in self.cds_exons)


_CHECKSUM_TABLE = "pie_cache"
_CHECKSUM_KEY = "checksum"


def _file_checksum(path: str) -> str:
    """Content checksum of the annotation file.

    Hashes the whole file. Size + mtime are not enough on their own: ``rsync
    -a``, ``cp -p``, ``tar -x`` and ``touch -r`` all preserve mtime, so an
    edited GFF could pass validation and a stale DB be reused silently.
    Hashing only a prefix has the same hole for any edit past the prefix.
    Full MD5 costs well under a second even for a 178 MB GFF, which is
    negligible next to building the DB it guards.
    """
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _resolve_cache_path(gff_path: str) -> str | None:
    """Return a writable cache path for the DB, or None if not possible."""
    candidate = gff_path + ".pie.db"
    # Check if we can write next to the GFF file
    parent = os.path.dirname(candidate) or "."
    if os.access(parent, os.W_OK):
        return candidate
    return None


def _read_cached_checksum(db: gffutils.FeatureDB) -> str | None:
    """Read the pie checksum from a cached DB, or None if unavailable."""
    conn = db.conn
    # Check if our table exists (don't query gffutils' own 'meta' table)
    row = conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
        (_CHECKSUM_TABLE,),
    ).fetchone()
    if not row:
        return None
    stored = conn.execute(
        f"SELECT value FROM {_CHECKSUM_TABLE} WHERE key=?",
        (_CHECKSUM_KEY,),
    ).fetchone()
    return stored[0] if stored else None


def _write_cached_checksum(db: gffutils.FeatureDB, checksum: str) -> None:
    """Write the pie checksum into the DB for future cache validation."""
    conn = db.conn
    conn.execute(
        f"CREATE TABLE IF NOT EXISTS {_CHECKSUM_TABLE} "
        "(key TEXT PRIMARY KEY, value TEXT)"
    )
    conn.execute(
        f"INSERT OR REPLACE INTO {_CHECKSUM_TABLE} VALUES (?, ?)",
        (_CHECKSUM_KEY, checksum),
    )
    conn.commit()


def _load_or_create_db(path: str) -> gffutils.FeatureDB:
    """Load a cached gffutils database or create one.

    Uses create_unique strategy (fast) with fallback to merge (slow)
    for files with duplicate feature IDs. When the GFF directory is
    writable, caches the result as a SQLite file for reuse across runs.
    Falls back to in-memory DB for read-only locations.
    """
    checksum = _file_checksum(path)
    cache_path = _resolve_cache_path(path)

    # Try to reuse cached database
    if cache_path and os.path.exists(cache_path):
        try:
            db = gffutils.FeatureDB(cache_path)
            if _read_cached_checksum(db) == checksum:
                log.info("Reusing cached annotation DB: %s", cache_path)
                return db
        except Exception:
            pass  # stale or corrupt — rebuild

    # Determine target: file-backed if writable, else in-memory
    dbfn = cache_path if cache_path else ":memory:"

    # Build fresh database; try fast strategy first
    for strategy in ("create_unique", "merge"):
        try:
            db = gffutils.create_db(
                path,
                dbfn=dbfn,
                force=True,
                keep_order=True,
                merge_strategy=strategy,
                sort_attribute_values=True,
            )
            break
        except Exception:
            if strategy == "merge":
                raise
            log.info("create_unique failed for %s, falling back to merge", path)

    # Store checksum for cache validation (only meaningful for file-backed DBs)
    if cache_path:
        _write_cached_checksum(db, checksum)
        log.info("Built annotation DB (%s strategy): %s", strategy, cache_path)
    else:
        log.info("Built in-memory annotation DB (%s strategy) — "
                 "GFF directory is read-only", strategy)
    return db


def _merge_intervals(
    exons: list[tuple[str, int, int]],
) -> list[tuple[str, int, int]]:
    """Merge overlapping/adjacent (chrom, start, end) intervals, ascending."""
    merged: list[tuple[str, int, int]] = []
    for chrom, start, end in sorted(exons):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((chrom, start, end))
    return merged


def _canonical_exons(
    cds: list[tuple[str, int, int, str]], strand: str
) -> tuple[list[tuple[str, int, int]], int]:
    """Normalise a transcript's CDS features into pie's exon list.

    Every consumer of `GeneModel.cds_exons` depends on both normalisations:

    * **Overlap.** Overlapping or repeated CDS features (a hand-edited GFF,
      or one listing a CDS twice) would otherwise have their bases
      concatenated twice by `extract_codons` — shifting the reading frame —
      and their variants fetched once per covering exon.
    * **Phase.** GFF3 column 8 gives, for the first CDS of a transcript, the
      number of bases to remove from its 5' end to reach the first complete
      codon. It is non-zero for 5'-partial models (no annotated start codon,
      or one running off a contig edge), routine in NCBI/EGAPx output, and
      reading straight through from the first base instead translates the
      whole transcript out of frame.

    Only the *first* CDS's phase is honoured. Internal CDS phases are assumed
    consistent with the exon lengths and are not used to resynchronise the
    reading frame: an annotation whose internal phases disagree with its exon
    structure is an annotation problem, and letting it surface as an internal
    stop codon (already counted and reported per gene) is more useful than
    silently patching the frame at every exon boundary and hiding it.

    Returns the ascending exon list and the phase that was applied.
    """
    ordered = sorted(cds, key=lambda e: e[1], reverse=(strand == "-"))
    frame = ordered[0][3]
    phase = int(frame) if frame in ("1", "2") else 0

    exons = _merge_intervals([(c, s, e) for c, s, e, _f in ordered])
    if strand == "-":
        exons.reverse()

    trimmed = []
    remaining = phase
    for chrom, start, end in exons:
        if remaining:
            n = min(remaining, end - start)
            remaining -= n
            if strand == "-":
                end -= n
            else:
                start += n
            if end <= start:
                continue
        trimmed.append((chrom, start, end))
    return sorted(trimmed, key=lambda e: e[1]), phase


def parse_annotations(path: str) -> list[GeneModel]:
    """Parse GFF3 or GTF, select longest isoform per gene.

    CDS exons are normalised by `_canonical_exons`: overlaps merged and the
    GFF `phase` of the transcript's first CDS honoured.

    Returns list sorted by (chrom, start).
    """
    db = _load_or_create_db(path)

    genes = []
    n_phased = 0
    for gene in db.all_features(featuretype="gene"):
        gene_id = gene.id
        chrom = gene.seqid
        strand = gene.strand

        # Collect CDS features per transcript, keeping the GFF phase column
        transcript_cds: dict[str, list[tuple[str, int, int, str]]] = defaultdict(list)
        for tx in db.children(
            gene, featuretype=["mRNA", "transcript"], order_by="start"
        ):
            for cds in db.children(tx, featuretype="CDS", order_by="start"):
                transcript_cds[tx.id].append(
                    (cds.seqid, cds.start - 1, cds.end, cds.frame))

        # Fallback: CDS directly under gene (no transcript intermediate)
        if not transcript_cds:
            for cds in db.children(gene, featuretype="CDS", order_by="start"):
                transcript_cds[gene_id].append(
                    (cds.seqid, cds.start - 1, cds.end, cds.frame))

        if not transcript_cds:
            continue

        # Pick transcript with longest total CDS
        best_tid = max(
            transcript_cds,
            key=lambda tid: sum(e - s for _, s, e, _f in transcript_cds[tid]),
        )
        best_exons, phase = _canonical_exons(transcript_cds[best_tid], strand)
        n_phased += bool(phase)
        if not best_exons:
            continue

        genes.append(
            GeneModel(
                gene_id=gene_id,
                transcript_id=best_tid,
                chrom=chrom,
                start=gene.start - 1,
                end=gene.end,
                strand=strand,
                cds_exons=best_exons,
            )
        )

    if n_phased:
        log.info(
            "Trimmed a 5' phase offset on %d of %d gene model(s) "
            "(5'-partial CDS)", n_phased, len(genes),
        )

    genes.sort(key=lambda g: (g.chrom, g.start))
    return genes
