"""GFF3/GTF parsing with longest-isoform selection."""

from collections import defaultdict
from dataclasses import dataclass

import gffutils


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


def _detect_format(path: str) -> str:
    """Detect whether a file is GFF3 or GTF."""
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                if "gff-version" in line.lower():
                    return "gff3"
                continue
            if "ID=" in line or "Parent=" in line:
                return "gff3"
            if 'gene_id "' in line or "gene_id '" in line:
                return "gtf"
    return "gff3"


def parse_annotations(path: str) -> list[GeneModel]:
    """Parse GFF3 or GTF, select longest isoform per gene.

    Returns list sorted by (chrom, start).
    """
    db = gffutils.create_db(
        path,
        dbfn=":memory:",
        force=True,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True,
    )

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
