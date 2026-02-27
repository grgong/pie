"""pie CLI — piN/piS Estimator for pool-seq and individual-sequencing data."""
import logging
import os
import sys

import click

from pie import __version__

log = logging.getLogger("pie")


@click.group()
@click.version_option(version=__version__)
def main():
    """pie — piN/piS Estimator for pool-seq and individual-sequencing data."""
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
@click.option("--mode", default="pool", show_default=True, help="Analysis mode: 'pool' (pool-seq) or 'individual'/'ind' (individual genotypes)")
@click.option("--min-freq", default=0.01, show_default=True, help="Minimum allele frequency")
@click.option("--min-depth", default=None, type=int, help="Minimum read depth [pool mode only, default: 10]")
@click.option("--min-qual", default=20.0, show_default=True, help="Minimum variant quality")
@click.option("--pass-only", is_flag=True, help="Only use PASS-filtered variants")
@click.option("--keep-multiallelic", is_flag=True, help="Keep and merge multiallelic sites instead of skipping them")
@click.option("--include-stop-codons", is_flag=True, help="Count stop_gained mutations as nonsynonymous (by default they are excluded, matching NG86/SNPGenie conventions)")
@click.option("--window-size", default=1000, show_default=True, help="Sliding window size (bp)")
@click.option("--window-step", default=100, show_default=True, type=click.IntRange(min=1), help="Sliding window step (bp)")
@click.option("--threads", default=1, show_default=True, help="Number of threads")
@click.option("--sample", default=None, help="Sample name to analyse [pool mode only, for multi-sample VCFs]")
@click.option("--samples", default=None, help="Comma-separated sample names [individual mode only]")
@click.option("--samples-file", default=None, type=click.Path(exists=True), help="File with one sample name per line [individual mode only]")
@click.option("--min-call-rate", default=None, type=float, help="Minimum genotype call rate [individual mode only, default: 0.8]")
@click.option("--min-an", default=None, type=int, help="Minimum allele number (AN) [individual mode only, default: 2]")
def run(vcf, gff, fasta, outdir, mode, min_freq, min_depth, min_qual,
        pass_only, keep_multiallelic, include_stop_codons, window_size,
        window_step, threads, sample, samples, samples_file, min_call_rate,
        min_an):
    """Run piN/piS analysis."""
    from pie.vcf import ensure_indexed, get_sample_names
    from pie.parallel import run_parallel
    from pie.io import write_gene_results, write_window_results, write_summary

    # --- Normalize mode ---
    if mode == "ind":
        mode = "individual"
    if mode not in ("pool", "individual"):
        click.echo(f"Error: --mode must be 'pool', 'individual', or 'ind', got '{mode}'", err=True)
        sys.exit(1)

    # --- Cross-option validation ---
    if mode == "pool":
        for opt_name, opt_val in [("--samples", samples), ("--samples-file", samples_file),
                                  ("--min-call-rate", min_call_rate), ("--min-an", min_an)]:
            if opt_val is not None:
                click.echo(
                    f"Error: {opt_name} is only valid in individual mode, not pool mode.",
                    err=True,
                )
                sys.exit(1)
        # Default min_depth for pool mode
        if min_depth is None:
            min_depth = 10
    else:  # individual
        if sample is not None:
            click.echo(
                "Error: --sample is only valid in pool mode. "
                "Use --samples for individual mode.",
                err=True,
            )
            sys.exit(1)
        if min_depth is not None:
            click.echo(
                "Error: --min-depth is only valid in pool mode, not individual mode.",
                err=True,
            )
            sys.exit(1)

    # --- Individual-mode specific validation ---
    if mode == "individual":
        if samples is not None and samples_file is not None:
            click.echo(
                "Error: --samples and --samples-file are mutually exclusive.",
                err=True,
            )
            sys.exit(1)

        if min_call_rate is not None and not (0.0 <= min_call_rate <= 1.0):
            click.echo(
                f"Error: --min-call-rate must be between 0 and 1, got {min_call_rate}",
                err=True,
            )
            sys.exit(1)

        # Apply defaults for individual mode
        if min_call_rate is None:
            min_call_rate = 0.8
        if min_an is None:
            min_an = 2
        min_depth = 0  # unused in individual mode but required by run_parallel

    # --- Validate inputs exist ---
    for path, name in [(vcf, "VCF"), (gff, "GFF"), (fasta, "FASTA")]:
        if not os.path.exists(path):
            click.echo(f"Error: {name} file not found: {path}", err=True)
            sys.exit(1)

    # --- Resolve sample list ---
    vcf_samples = get_sample_names(vcf)
    selected_samples = None

    if mode == "pool":
        # Pool mode: existing --sample validation
        if sample is not None:
            if sample not in vcf_samples:
                click.echo(
                    f"Error: sample '{sample}' not found in VCF. "
                    f"Available samples: {', '.join(vcf_samples)}",
                    err=True,
                )
                sys.exit(1)
            log.info("Using sample: %s", sample)
        elif len(vcf_samples) >= 2:
            click.echo(
                f"Error: VCF contains {len(vcf_samples)} samples. "
                f"Use --sample to specify one.\n"
                f"Available samples: {', '.join(vcf_samples)}",
                err=True,
            )
            sys.exit(1)
    else:  # individual
        if samples is not None:
            selected_samples = [s.strip() for s in samples.split(",") if s.strip()]
        elif samples_file is not None:
            with open(samples_file) as fh:
                selected_samples = [line.strip() for line in fh if line.strip()]
        else:
            selected_samples = list(vcf_samples)  # use all VCF samples
            log.info("Using all %d samples from VCF", len(selected_samples))

        if not selected_samples:
            click.echo("Error: no samples specified or resolved.", err=True)
            sys.exit(1)

        # Validate sample names exist in VCF
        if samples is not None or samples_file is not None:
            missing = [s for s in selected_samples if s not in vcf_samples]
            if missing:
                click.echo(
                    f"Error: sample(s) not found in VCF: {', '.join(missing)}. "
                    f"Available samples: {', '.join(vcf_samples)}",
                    err=True,
                )
                sys.exit(1)
            log.info("Using %d samples: %s", len(selected_samples), ", ".join(selected_samples))

    # Create output directory
    os.makedirs(outdir, exist_ok=True)

    # Ensure VCF is indexed
    vcf = ensure_indexed(vcf)
    log.info("VCF ready: %s", vcf)

    # Run analysis
    log.info("Starting piN/piS analysis with %d thread(s) in %s mode", threads, mode)
    results = run_parallel(
        fasta_path=fasta, gff_path=gff, vcf_path=vcf,
        min_freq=min_freq, min_depth=min_depth, min_qual=min_qual,
        pass_only=pass_only, keep_multiallelic=keep_multiallelic,
        exclude_stops=not include_stop_codons, threads=threads,
        sample=sample, mode=mode, samples=selected_samples,
        min_call_rate=min_call_rate, min_an=min_an,
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
