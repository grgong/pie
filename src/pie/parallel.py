"""Gene-level multiprocessing for parallel piN/piS computation."""

import logging
from multiprocessing import Pool

from pie.annotation import GeneModel, NoGenesFoundError, parse_annotations
from pie.reference import ReferenceGenome
from pie.vcf import IndividualVariantReader, VariantReader, get_vcf_contigs
from pie.diversity import compute_gene_diversity, GeneResult

log = logging.getLogger(__name__)


def _worker_init(fasta_path, vcf_path, min_freq, min_depth, min_qual,
                 pass_only, keep_multiallelic, exclude_stops, sample,
                 mode, samples, min_call_rate, min_an):
    """Initialize per-worker file handles (stored in globals)."""
    global _ref, _vcf, _exclude_stops, _n_samples
    _ref = ReferenceGenome(fasta_path)
    _exclude_stops = exclude_stops

    if mode == "individual":
        _vcf = IndividualVariantReader(
            vcf_path, samples=samples, min_freq=min_freq,
            min_qual=min_qual, pass_only=pass_only,
            keep_multiallelic=keep_multiallelic,
            min_call_rate=min_call_rate, min_an=min_an,
        )
        _n_samples = _vcf.n_samples
    else:
        _vcf = VariantReader(
            vcf_path, min_freq=min_freq, min_depth=min_depth,
            min_qual=min_qual, pass_only=pass_only,
            keep_multiallelic=keep_multiallelic, sample=sample,
        )
        _n_samples = None


def _worker_cleanup():
    """Close per-worker file handles."""
    global _ref, _vcf
    _ref.close()
    _vcf.close()


def _process_gene(gene: GeneModel) -> GeneResult:
    """Process a single gene using worker-local handles."""
    try:
        result = compute_gene_diversity(gene, _ref, _vcf, exclude_stops=_exclude_stops)
        result.n_samples = _n_samples
        return result
    except Exception as exc:
        raise RuntimeError(
            f"Failed processing gene {gene.gene_id} "
            f"({gene.chrom}:{gene.start}-{gene.end}): {exc!r}"
        ) from exc


def run_parallel(
    fasta_path: str,
    gff_path: str,
    vcf_path: str,
    min_freq: float = 0.01,
    min_depth: int = 10,
    min_qual: float = 20.0,
    pass_only: bool = False,
    keep_multiallelic: bool = False,
    exclude_stops: bool = True,
    threads: int = 1,
    sample: str | None = None,
    mode: str = "pool",
    samples: list[str] | None = None,
    min_call_rate: float = 0.8,
    min_an: int = 2,
) -> list[GeneResult]:
    """Run piN/piS analysis across all genes.

    Returns list of GeneResult sorted by (chrom, start).
    """
    genes = parse_annotations(gff_path)
    log.info("Parsed %d genes from %s", len(genes), gff_path)

    if not genes:
        raise NoGenesFoundError(
            f"No genes with CDS features found in {gff_path}. "
            "Ensure the annotation contains 'gene' features with child CDS records."
        )

    # Pre-flight: verify contig name overlap between annotation and VCF.
    # Skip when VCF has no ##contig headers (seqnames empty) — analysis
    # proceeds normally and fetch() returns no variants for each gene.
    gene_contigs = {g.chrom for g in genes}
    vcf_contigs = get_vcf_contigs(vcf_path)
    if vcf_contigs:
        shared = gene_contigs & vcf_contigs
        if not shared:
            missing_sample = sorted(gene_contigs)[:5]
            vcf_sample = sorted(vcf_contigs)[:5]
            raise ValueError(
                f"No contig names shared between annotation and VCF. "
                f"Annotation contigs (first 5): {missing_sample}; "
                f"VCF contigs (first 5): {vcf_sample}. "
                f"Check for naming mismatches (e.g. 'chr1' vs '1')."
            )
        missing = gene_contigs - vcf_contigs
        if missing:
            n_before = len(genes)
            genes = [g for g in genes if g.chrom in shared]
            log.warning(
                "%d of %d annotation contig(s) absent from VCF (%s); "
                "dropped %d genes on those contigs (%d genes remaining)",
                len(missing), len(gene_contigs),
                ", ".join(sorted(missing)[:10])
                + (" ..." if len(missing) > 10 else ""),
                n_before - len(genes), len(genes),
            )

    init_args = (fasta_path, vcf_path, min_freq, min_depth, min_qual,
                 pass_only, keep_multiallelic, exclude_stops, sample,
                 mode, samples, min_call_rate, min_an)

    if threads <= 1:
        # Single-threaded: no multiprocessing overhead
        _worker_init(*init_args)
        try:
            results = [_process_gene(g) for g in genes]
        finally:
            _worker_cleanup()
    else:
        with Pool(
            processes=threads,
            initializer=_worker_init,
            initargs=init_args,
        ) as pool:
            results = pool.map(_process_gene, genes)

    results.sort(key=lambda r: (r.chrom, r.start))
    return results
