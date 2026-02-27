"""GFF3/GTF parsing with longest-isoform selection."""

import hashlib
import logging
import os
from collections import defaultdict
from dataclasses import dataclass

import gffutils

log = logging.getLogger(__name__)


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


def _file_checksum(path: str) -> str:
    """Fast checksum using file size + mtime + first 8KB hash."""
    stat = os.stat(path)
    h = hashlib.md5()
    h.update(f"{stat.st_size}:{stat.st_mtime_ns}".encode())
    with open(path, "rb") as f:
        h.update(f.read(8192))
    return h.hexdigest()


_CHECKSUM_TABLE = "pie_cache"
_CHECKSUM_KEY = "checksum"


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


def parse_annotations(path: str) -> list[GeneModel]:
    """Parse GFF3 or GTF, select longest isoform per gene.

    Returns list sorted by (chrom, start).
    """
    db = _load_or_create_db(path)

    genes = []
    for gene in db.all_features(featuretype="gene"):
        gene_id = gene.id
        chrom = gene.seqid
        strand = gene.strand

        # Collect CDS features per transcript
        transcript_cds: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
        for tx in db.children(
            gene, featuretype=["mRNA", "transcript"], order_by="start"
        ):
            for cds in db.children(tx, featuretype="CDS", order_by="start"):
                transcript_cds[tx.id].append((cds.seqid, cds.start - 1, cds.end))

        # Fallback: CDS directly under gene (no transcript intermediate)
        if not transcript_cds:
            for cds in db.children(gene, featuretype="CDS", order_by="start"):
                transcript_cds[gene_id].append((cds.seqid, cds.start - 1, cds.end))

        if not transcript_cds:
            continue

        # Pick transcript with longest total CDS
        best_tid = max(
            transcript_cds,
            key=lambda tid: sum(e - s for _, s, e in transcript_cds[tid]),
        )
        best_exons = sorted(transcript_cds[best_tid], key=lambda x: x[1])

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

    genes.sort(key=lambda g: (g.chrom, g.start))
    return genes
