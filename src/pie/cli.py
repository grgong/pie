"""pie CLI — piN/piS Estimator for pool-seq data."""
import logging
import os
import sys

import click

from pie import __version__

log = logging.getLogger("pie")


@click.group()
@click.version_option(version=__version__)
def main():
    """pie — piN/piS Estimator for pool-seq data."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )


@main.command()
@click.option("--vcf", required=True, help="Input VCF file (bgzipped or plain)")
@click.option("--gff", required=True, help="GFF3 or GTF annotation file")
@click.option("--fasta", required=True, help="Reference FASTA file (indexed)")
@click.option("--outdir", required=True, help="Output directory")
@click.option("--min-freq", default=0.01, show_default=True, help="Minimum allele frequency")
@click.option("--min-depth", default=10, show_default=True, help="Minimum read depth")
@click.option("--min-qual", default=20.0, show_default=True, help="Minimum variant quality")
@click.option("--pass-only", is_flag=True, help="Only use PASS-filtered variants")
@click.option("--keep-multiallelic", is_flag=True, help="Keep and merge multiallelic sites instead of skipping them")
@click.option("--exclude-stop-codons", is_flag=True, help="Exclude stop_gained mutations from piN (by default they count as nonsynonymous)")
@click.option("--window-size", default=1000, show_default=True, help="Sliding window size (bp)")
@click.option("--window-step", default=100, show_default=True, type=click.IntRange(min=1), help="Sliding window step (bp)")
@click.option("--threads", default=1, show_default=True, help="Number of threads")
@click.option("--sample", default=None, help="Sample name to analyse (required for multi-sample VCFs)")
def run(vcf, gff, fasta, outdir, min_freq, min_depth, min_qual, pass_only,
        keep_multiallelic, exclude_stop_codons, window_size, window_step,
        threads, sample):
    """Run piN/piS analysis."""
    from pie.vcf import ensure_indexed, get_sample_names
    from pie.parallel import run_parallel
    from pie.io import write_gene_results, write_window_results, write_summary

    # Validate inputs exist
    for path, name in [(vcf, "VCF"), (gff, "GFF"), (fasta, "FASTA")]:
        if not os.path.exists(path):
            click.echo(f"Error: {name} file not found: {path}", err=True)
            sys.exit(1)

    # Validate sample selection
    samples = get_sample_names(vcf)
    if sample is not None:
        if sample not in samples:
            click.echo(
                f"Error: sample '{sample}' not found in VCF. "
                f"Available samples: {', '.join(samples)}",
                err=True,
            )
            sys.exit(1)
        log.info("Using sample: %s", sample)
    elif len(samples) >= 2:
        click.echo(
            f"Error: VCF contains {len(samples)} samples. "
            f"Use --sample to specify one.\n"
            f"Available samples: {', '.join(samples)}",
            err=True,
        )
        sys.exit(1)

    # Create output directory
    os.makedirs(outdir, exist_ok=True)

    # Ensure VCF is indexed
    vcf = ensure_indexed(vcf)
    log.info("VCF ready: %s", vcf)

    # Run analysis
    log.info("Starting piN/piS analysis with %d thread(s)", threads)
    results = run_parallel(
        fasta_path=fasta, gff_path=gff, vcf_path=vcf,
        min_freq=min_freq, min_depth=min_depth, min_qual=min_qual,
        pass_only=pass_only, keep_multiallelic=keep_multiallelic,
        exclude_stops=exclude_stop_codons, threads=threads, sample=sample,
    )
    log.info("Processed %d genes", len(results))

    # Write outputs (prefix with sample name when --sample is given)
    prefix = f"{sample}." if sample else ""

    gene_path = os.path.join(outdir, f"{prefix}gene_results.tsv")
    write_gene_results(results, gene_path)
    log.info("Gene results: %s", gene_path)

    if window_size > 0:
        win_path = os.path.join(outdir, f"{prefix}window_results.tsv")
        write_window_results(results, win_path, window_size, window_step)
        log.info("Window results: %s", win_path)

    summary_path = os.path.join(outdir, f"{prefix}summary.tsv")
    write_summary(results, summary_path)
    log.info("Summary: %s", summary_path)

    click.echo(f"Done. Results written to {outdir}/")


@main.command()
@click.option("--gene-results", required=True, help="Path to gene_results.tsv")
@click.option("--output", required=True, help="Output plot path (PNG)")
@click.option("--width", default=16.0, show_default=True, help="Figure width (inches)")
@click.option("--height", default=6.0, show_default=True, help="Figure height (inches)")
def plot(gene_results, output, width, height):
    """Create Manhattan plot from results."""
    from pie.plot import manhattan_plot

    if not os.path.exists(gene_results):
        click.echo(f"Error: file not found: {gene_results}", err=True)
        sys.exit(1)

    manhattan_plot(gene_results, output, width=width, height=height)
    click.echo(f"Plot saved to {output}")


@main.command()
@click.argument("summary_file")
def summary(summary_file):
    """Print summary statistics from summary.tsv."""
    import pandas as pd

    if not os.path.exists(summary_file):
        click.echo(f"Error: file not found: {summary_file}", err=True)
        sys.exit(1)

    df = pd.read_csv(summary_file, sep="\t")
    for col in df.columns:
        val = df.iloc[0][col]
        if isinstance(val, float):
            click.echo(f"  {col}: {val:.6f}")
        else:
            click.echo(f"  {col}: {val}")
