"""pie CLI — piN/piS Estimator for pool-seq and individual-sequencing data."""
import logging
import os
import sys

import click

from pie import __version__

log = logging.getLogger("pie")

_HELP_OPTS = {"help_option_names": ["-h", "--help"]}


class _OrderedGroup(click.Group):
    """A Click group that lists commands in registration order."""

    def list_commands(self, ctx):
        return list(self.commands)


def _apply_options(f, options):
    """Apply Click options while preserving the declared help order."""
    for args, kwargs in reversed(options):
        f = click.option(*args, **kwargs)(f)
    return f


@click.group(cls=_OrderedGroup, context_settings=_HELP_OPTS)
@click.version_option(__version__, "-V", "--version")
def main():
    """pie — piN/piS Estimator for pool-seq and individual-sequencing data.

    \b
    Compute per-gene and sliding-window piN/piS ratios from a VCF, a genome
    annotation (GFF3/GTF), and a reference FASTA.  Supports both pool-seq
    (allele-frequency based) and individual-sequencing (genotype based) data.

    \b
    Quick start:
      pie pool -v variants.vcf.gz -g genes.gff3 -f ref.fa -o results/
      pie ind  -v gatk.vcf.gz -g genes.gff3 -f ref.fa -o results/
      pie plot manhattan -i results/gene_results.tsv -o manhattan.png
      pie summary results/summary.tsv

    \b
    Use 'pie <command> -h' for command-specific help.
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )


# ---------------------------------------------------------------------------
# Shared run options
# ---------------------------------------------------------------------------

def _shared_run_options(f):
    """Decorator: shared options for both pool and ind subcommands."""
    return _apply_options(f, [
        (("-v", "--vcf"), {"required": True, "help": "Input VCF file (bgzipped or plain)."}),
        (("-g", "--gff"), {"required": True, "help": "GFF3 or GTF annotation file."}),
        (("-f", "--fasta"), {"required": True, "help": "Reference FASTA file (must be indexed with .fai)."}),
        (("-o", "--outdir"), {"required": True, "help": "Output directory (created if absent)."}),
        (("--min-freq",), {"default": 0.01, "show_default": True, "help": "Minimum allele frequency."}),
        (("-q", "--min-qual"), {"default": 20.0, "show_default": True, "help": "Minimum variant quality (QUAL)."}),
        (("--pass-only",), {"is_flag": True, "help": "Only use PASS-filtered variants."}),
        (
            ("--keep-multiallelic",),
            {"is_flag": True, "help": "Keep and merge multiallelic sites instead of skipping them."},
        ),
        (
            ("--include-stop-codons",),
            {
                "is_flag": True,
                "help": "Count stop_gained as nonsynonymous (excluded by default, matching NG86/SNPGenie).",
            },
        ),
        (("-w", "--window-size"), {"default": 1000, "show_default": True, "help": "Sliding window size in bp."}),
        (
            ("-W", "--window-step"),
            {
                "default": 100,
                "show_default": True,
                "type": click.IntRange(min=1),
                "help": "Sliding window step in bp.",
            },
        ),
        (("-t", "--threads"), {"default": 1, "show_default": True, "help": "Number of parallel threads."}),
        (("--quiet",), {"is_flag": True, "help": "Suppress progress messages (show only warnings and summary)."}),
    ])


# ---------------------------------------------------------------------------
# Shared analysis logic
# ---------------------------------------------------------------------------

def _run_analysis(*, vcf, gff, fasta, outdir, mode, min_freq, min_depth,
                  min_qual, pass_only, keep_multiallelic, include_stop_codons,
                  window_size, window_step, threads, quiet, sample=None,
                  samples=None, min_call_rate=None, min_an=None):
    """Core analysis logic shared by pool and ind subcommands."""
    from pie.vcf import ensure_indexed, get_sample_names
    from pie.parallel import run_parallel
    from pie.annotation import NoGenesFoundError
    from pie.io import write_gene_results, write_window_results, write_summary

    if quiet:
        logging.getLogger("pie").setLevel(logging.WARNING)

    # --- Validate inputs exist ---
    for path, name in [(vcf, "VCF"), (gff, "GFF"), (fasta, "FASTA")]:
        if not os.path.exists(path):
            click.echo(f"Error: {name} file not found: {path}", err=True)
            sys.exit(1)

    # --- Resolve and validate samples before any filesystem mutations ---
    vcf_samples = get_sample_names(vcf)
    selected_samples = None

    if mode == "pool":
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
            selected_samples = samples
        else:
            selected_samples = list(vcf_samples)
            log.info("Using all %d samples from VCF", len(selected_samples))

        if not selected_samples:
            click.echo("Error: no samples specified or resolved.", err=True)
            sys.exit(1)

        # Deduplicate sample names (preserving order)
        deduped = list(dict.fromkeys(selected_samples))
        if len(deduped) < len(selected_samples):
            log.warning("Duplicate sample names removed: %d -> %d unique",
                        len(selected_samples), len(deduped))
            selected_samples = deduped

        # Validate sample names exist in VCF
        if samples is not None:
            missing = [s for s in selected_samples if s not in vcf_samples]
            if missing:
                click.echo(
                    f"Error: sample(s) not found in VCF: {', '.join(missing)}. "
                    f"Available samples: {', '.join(vcf_samples)}",
                    err=True,
                )
                sys.exit(1)
            log.info("Using %d samples: %s", len(selected_samples), ", ".join(selected_samples))

    # Ensure VCF is indexed (bgzip + tabix if needed) — only after validation
    vcf = ensure_indexed(vcf)
    log.info("VCF ready: %s", vcf)

    # Create output directory
    os.makedirs(outdir, exist_ok=True)

    # Run analysis
    log.info("Starting piN/piS analysis with %d thread(s) in %s mode", threads, mode)
    try:
        results = run_parallel(
            fasta_path=fasta, gff_path=gff, vcf_path=vcf,
            min_freq=min_freq, min_depth=min_depth, min_qual=min_qual,
            pass_only=pass_only, keep_multiallelic=keep_multiallelic,
            exclude_stops=not include_stop_codons, threads=threads,
            sample=sample, mode=mode, samples=selected_samples,
            min_call_rate=min_call_rate, min_an=min_an,
        )
    except (NoGenesFoundError, ValueError) as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)
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


# ---------------------------------------------------------------------------
# pie pool
# ---------------------------------------------------------------------------

@main.command(no_args_is_help=True, context_settings=_HELP_OPTS)
@_shared_run_options
@click.option("-d", "--min-depth", default=10, show_default=True, type=int,
              help="Minimum read depth.")
@click.option("-s", "--sample", default=None,
              help="Sample name to analyse (for multi-sample VCFs).")
def pool(vcf, gff, fasta, outdir, min_freq, min_qual, pass_only,
         keep_multiallelic, include_stop_codons, window_size, window_step,
         threads, quiet, min_depth, sample):
    """Run pool-seq piN/piS analysis (allele-frequency based).

    \b
    Examples:
      # Basic usage:
      pie pool -v pool.vcf.gz -g genes.gff3 -f ref.fa -o results/

    \b
      # Stricter filters and 8 threads:
      pie pool -v pool.vcf.gz -g genes.gff3 -f ref.fa -o results/ \\
          -d 20 -q 30 --pass-only -t 8

    \b
      # Multi-sample pool VCF — select one sample:
      pie pool -v multi.vcf.gz -g genes.gff3 -f ref.fa -o out/ -s SampleA
    """
    _run_analysis(
        vcf=vcf, gff=gff, fasta=fasta, outdir=outdir, mode="pool",
        min_freq=min_freq, min_depth=min_depth, min_qual=min_qual,
        pass_only=pass_only, keep_multiallelic=keep_multiallelic,
        include_stop_codons=include_stop_codons, window_size=window_size,
        window_step=window_step, threads=threads, quiet=quiet, sample=sample,
    )


# ---------------------------------------------------------------------------
# pie ind
# ---------------------------------------------------------------------------

@main.command(no_args_is_help=True, context_settings=_HELP_OPTS)
@_shared_run_options
@click.option("-S", "--samples", default=None,
              help="Comma-separated sample names.")
@click.option("--samples-file", default=None, type=click.Path(exists=True, dir_okay=False),
              help="File with one sample name per line.")
@click.option("--min-call-rate", default=0.8, show_default=True, type=click.FloatRange(0.0, 1.0),
              help="Minimum genotype call rate.")
@click.option("--min-an", default=2, show_default=True, type=int,
              help="Minimum allele number (AN).")
def ind(vcf, gff, fasta, outdir, min_freq, min_qual, pass_only,
        keep_multiallelic, include_stop_codons, window_size, window_step,
        threads, quiet, samples, samples_file, min_call_rate, min_an):
    """Run individual-sequencing piN/piS analysis (genotype based).

    \b
    Examples:
      # All samples in VCF:
      pie ind -v gatk.vcf.gz -g genes.gff3 -f ref.fa -o out/

    \b
      # Subset of samples:
      pie ind -v gatk.vcf.gz -g genes.gff3 -f ref.fa -o out/ \\
          -S sampleA,sampleB,sampleC

    \b
      # Samples from a file:
      pie ind -v gatk.vcf.gz -g genes.gff3 -f ref.fa -o out/ \\
          --samples-file sample_list.txt --min-call-rate 0.9
    """
    # --- Individual-mode specific validation ---
    if samples is not None and samples_file is not None:
        click.echo(
            "Error: --samples and --samples-file are mutually exclusive.",
            err=True,
        )
        sys.exit(1)

    # Resolve samples to a list
    if samples is not None:
        resolved_samples = [s.strip() for s in samples.split(",") if s.strip()]
    elif samples_file is not None:
        with open(samples_file) as fh:
            resolved_samples = [line.strip() for line in fh if line.strip()]
    else:
        resolved_samples = None  # _run_analysis will use all VCF samples

    _run_analysis(
        vcf=vcf, gff=gff, fasta=fasta, outdir=outdir, mode="individual",
        min_freq=min_freq, min_depth=0, min_qual=min_qual,
        pass_only=pass_only, keep_multiallelic=keep_multiallelic,
        include_stop_codons=include_stop_codons, window_size=window_size,
        window_step=window_step, threads=threads, quiet=quiet,
        samples=resolved_samples, min_call_rate=min_call_rate, min_an=min_an,
    )


# ---------------------------------------------------------------------------
# pie plot (group + subcommands)
# ---------------------------------------------------------------------------

@main.group(cls=_OrderedGroup, invoke_without_command=True, no_args_is_help=True, context_settings=_HELP_OPTS)
def plot():
    """Create publication-ready plots from piN/piS results.
    """


def _shared_plot_options(f):
    """Decorator: shared options for all plot subcommands (except sliding-window)."""
    return _apply_options(f, [
        (("-i", "--input", "input_path"), {"required": True, "type": click.Path(exists=True, dir_okay=False), "help": "Input TSV file."}),
        (("-o", "--output", "output_path"), {"required": True, "help": "Output plot path (PNG/PDF/SVG)."}),
        (("--dpi",), {"default": 300, "show_default": True, "help": "Resolution in dots per inch."}),
        (("--min-codons",), {"default": None, "type": int, "help": "Exclude genes with fewer than N codons."}),
        (("--min-variants",), {"default": None, "type": int, "help": "Exclude genes with fewer than N variants."}),
    ])


def _run_plot(plot_fn, input_path, output_path, **kwargs):
    """Shared wrapper for all plot subcommands: run and report."""
    plot_fn(input_path, output_path, **kwargs)
    click.echo(f"Plot saved to {output_path}")


@plot.command(no_args_is_help=True, context_settings=_HELP_OPTS)
@_shared_plot_options
@click.option("-W", "--width", default=7.2, show_default=True, help="Figure width in inches.")
@click.option("-H", "--height", default=3.5, show_default=True, help="Figure height in inches.")
@click.option("--log-scale", is_flag=True, help="Use log2 scale for piN/piS y-axis.")
@click.option("--label-top", default=None, type=int, help="Label top N outlier genes.")
@click.option("--highlight-genes", default=None, help="Comma-separated gene IDs to label.")
@click.option("--max-ratio", default=2.0, show_default=True, help="Exclude genes with piN/piS above this value.")
@click.option("--exclude-zero-ratio", is_flag=True, help="Exclude genes with piN/piS = 0 (robustness check).")
def manhattan(input_path, output_path, width, height, dpi, min_codons, min_variants,
              log_scale, label_top, highlight_genes, max_ratio, exclude_zero_ratio):
    """Create Manhattan plot of per-gene piN/piS.

    \b
    Genes are ordered by chromosomal position along the x-axis.  A dashed
    line marks neutral expectation (piN/piS = 1).

    \b
    Examples:
      pie plot manhattan -i gene_results.tsv -o manhattan.png
      pie plot manhattan -i gene_results.tsv -o manhattan.png --log-scale
      pie plot manhattan -i gene_results.tsv -o manhattan.png --label-top 10
    """
    from pie.plot import manhattan_plot as _manhattan_plot

    genes = [g.strip() for g in highlight_genes.split(",")] if highlight_genes else None
    _run_plot(_manhattan_plot, input_path, output_path, width=width, height=height, dpi=dpi,
              log_scale=log_scale, label_top=label_top, highlight_genes=genes,
              max_ratio=max_ratio, exclude_zero_ratio=exclude_zero_ratio,
              min_codons=min_codons, min_variants=min_variants)


@plot.command(no_args_is_help=True, context_settings=_HELP_OPTS)
@_shared_plot_options
@click.option("-W", "--width", default=4.7, show_default=True, help="Figure width in inches.")
@click.option("-H", "--height", default=4.7, show_default=True, help="Figure height in inches.")
@click.option("--color-by-chrom", is_flag=True, help="Color points by chromosome.")
@click.option("--max-piN", default=2.0, show_default=True, help="Exclude genes with piN above this value.")
@click.option("--max-piS", default=2.0, show_default=True, help="Exclude genes with piS above this value.")
def scatter(input_path, output_path, width, height, dpi, min_codons, min_variants,
            color_by_chrom, max_pin, max_pis):
    """Create piN vs piS scatter plot.

    \b
    Each gene is a point with piS on the x-axis and piN on the y-axis.
    Point size reflects gene length (codons).  A diagonal line marks
    piN = piS (neutral expectation).

    \b
    Examples:
      pie plot scatter -i gene_results.tsv -o scatter.png
      pie plot scatter -i gene_results.tsv -o scatter.png --color-by-chrom
    """
    from pie.plot import scatter_plot as _scatter_plot

    _run_plot(_scatter_plot, input_path, output_path, width=width, height=height, dpi=dpi,
              color_by_chrom=color_by_chrom, max_piN=max_pin, max_piS=max_pis,
              min_codons=min_codons, min_variants=min_variants)


@plot.command(no_args_is_help=True, context_settings=_HELP_OPTS)
@_shared_plot_options
@click.option("-W", "--width", default=3.5, show_default=True, help="Figure width in inches.")
@click.option("-H", "--height", default=3.0, show_default=True, help="Figure height in inches.")
@click.option("--max-ratio", default=2.0, show_default=True, help="Exclude genes with piN/piS above this value.")
@click.option("--exclude-zero-ratio", is_flag=True, help="Exclude genes with piN/piS = 0 (robustness check).")
def histogram(input_path, output_path, width, height, dpi, min_codons, min_variants,
              max_ratio, exclude_zero_ratio):
    """Create histogram of piN/piS distribution.

    \b
    Shows the distribution of piN/piS ratios across all genes, with a
    density curve overlay and a vertical line at piN/piS = 1.

    \b
    Examples:
      pie plot histogram -i gene_results.tsv -o histogram.png
    """
    from pie.plot import histogram_plot as _histogram_plot

    _run_plot(_histogram_plot, input_path, output_path, width=width, height=height, dpi=dpi,
              max_ratio=max_ratio, exclude_zero_ratio=exclude_zero_ratio,
              min_codons=min_codons, min_variants=min_variants)


@plot.command(no_args_is_help=True, context_settings=_HELP_OPTS)
@_shared_plot_options
@click.option("-W", "--width", default=3.5, show_default=True, help="Figure width in inches.")
@click.option("-H", "--height", default=7.0, show_default=True, help="Figure height in inches.")
@click.option("--max-piN", default=2.0, show_default=True, help="Exclude genes with piN above this value.")
@click.option("--max-piS", default=2.0, show_default=True, help="Exclude genes with piS above this value.")
@click.option("--max-ratio", default=2.0, show_default=True, help="Exclude genes with piN/piS above this value.")
@click.option("--exclude-zero-ratio", is_flag=True, help="Exclude genes with piN/piS = 0 (robustness check).")
def boxplot(input_path, output_path, width, height, dpi, min_codons, min_variants,
            max_pin, max_pis, max_ratio, exclude_zero_ratio):
    """Create per-chromosome boxplots of piN, piS, and piN/piS.

    \b
    Three faceted panels showing the distribution of each metric per
    chromosome.  Useful for identifying chromosomes under different
    selection pressures.

    \b
    Examples:
      pie plot boxplot -i gene_results.tsv -o boxplot.png
    """
    from pie.plot import boxplot_plot as _boxplot_plot

    _run_plot(_boxplot_plot, input_path, output_path, width=width, height=height, dpi=dpi,
              max_piN=max_pin, max_piS=max_pis, max_ratio=max_ratio,
              exclude_zero_ratio=exclude_zero_ratio,
              min_codons=min_codons, min_variants=min_variants)


@plot.command("sliding-window", no_args_is_help=True, context_settings=_HELP_OPTS)
@click.option("-i", "--input", "input_path", required=True, type=click.Path(exists=True, dir_okay=False), help="Input TSV file.")
@click.option("-o", "--output", "output_path", required=True, help="Output plot path (PNG/PDF/SVG).")
@click.option("-W", "--width", default=7.2, show_default=True, help="Figure width in inches.")
@click.option("-H", "--height", default=None, type=float, help="Figure height in inches.  [default: 1.5 per chromosome]")
@click.option("--dpi", default=300, show_default=True, help="Resolution in dots per inch.")
@click.option("--max-ratio", default=2.0, show_default=True, help="Exclude windows with piN/piS above this value.")
def sliding_window(input_path, output_path, width, height, dpi,
                   max_ratio):
    """Create sliding window piN/piS line plot.

    \b
    Line plot of piN/piS along genomic coordinates, faceted by chromosome.
    Input should be window_results.tsv from 'pie pool' or 'pie ind'.

    \b
    Examples:
      pie plot sliding-window -i window_results.tsv -o sw.png
    """
    from pie.plot import sliding_window_plot as _sliding_window_plot

    _run_plot(_sliding_window_plot, input_path, output_path, width=width, height=height, dpi=dpi,
              max_ratio=max_ratio)


# ---------------------------------------------------------------------------
# pie summary
# ---------------------------------------------------------------------------

@main.command(context_settings=_HELP_OPTS)
@click.argument("summary_file")
def summary(summary_file):
    """Print summary statistics from summary.tsv.

    \b
    Display genome-wide weighted-average piN, piS, and piN/piS from a
    summary file produced by 'pie pool' or 'pie ind'.

    \b
    Examples:
      pie summary results/summary.tsv
    """
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
