"""Gene-level multiprocessing for parallel piN/piS computation."""

import logging
from multiprocessing import Pool

from pie.annotation import GeneModel, parse_annotations
from pie.reference import ReferenceGenome
from pie.vcf import VariantReader
from pie.diversity import compute_gene_diversity, GeneResult

log = logging.getLogger(__name__)


def _worker_init(fasta_path, vcf_path, min_freq, min_depth, min_qual, pass_only):
    """Initialize per-worker file handles (stored in globals)."""
    global _ref, _vcf
    _ref = ReferenceGenome(fasta_path)
    _vcf = VariantReader(vcf_path, min_freq=min_freq, min_depth=min_depth,
                         min_qual=min_qual, pass_only=pass_only)


def _worker_cleanup():
    """Close per-worker file handles."""
    global _ref, _vcf
    _ref.close()
    _vcf.close()


def _process_gene(gene: GeneModel) -> GeneResult:
    """Process a single gene using worker-local handles."""
    return compute_gene_diversity(gene, _ref, _vcf)


def run_parallel(
    fasta_path: str,
    gff_path: str,
    vcf_path: str,
    min_freq: float = 0.01,
    min_depth: int = 10,
    min_qual: float = 20.0,
    pass_only: bool = False,
    threads: int = 1,
) -> list[GeneResult]:
    """Run piN/piS analysis across all genes.

    Returns list of GeneResult sorted by (chrom, start).
    """
    genes = parse_annotations(gff_path)
    log.info("Parsed %d genes from %s", len(genes), gff_path)

    if threads <= 1:
        # Single-threaded: no multiprocessing overhead
        _worker_init(fasta_path, vcf_path, min_freq, min_depth, min_qual, pass_only)
        try:
            results = [_process_gene(g) for g in genes]
        finally:
            _worker_cleanup()
    else:
        with Pool(
            processes=threads,
            initializer=_worker_init,
            initargs=(fasta_path, vcf_path, min_freq, min_depth, min_qual, pass_only),
        ) as pool:
            results = pool.map(_process_gene, genes)

    results.sort(key=lambda r: (r.chrom, r.start))
    return results
