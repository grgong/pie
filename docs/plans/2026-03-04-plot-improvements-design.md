# Plot Improvements Design

## Goal

Replace the minimal Manhattan-only plotting with a full suite of publication-ready visualizations using plotnine (ggplot2 for Python), exposed as `pie plot <subcommand>`.

## Plot Types

### Manhattan (`pie plot manhattan`)
- Genes ordered by cumulative chromosomal position on x-axis, piN/piS on y-axis
- Alternating Okabe-Ito colorblind-safe colors per chromosome
- Dashed neutral expectation line at piN/piS = 1
- `--log-scale`: log2 transform (neutral line at 0)
- `--label-top N`: annotate top N outlier genes
- `--highlight-genes gene1,gene2`: label specific genes
- Input: `gene_results.tsv`

### Scatter (`pie plot scatter`)
- piS on x-axis, piN on y-axis, one point per gene
- Diagonal line at piN = piS (neutral expectation)
- Point size scaled by n_codons
- `--color-by-chrom`: color points by chromosome
- Input: `gene_results.tsv`

### Histogram (`pie plot histogram`)
- Distribution of piN/piS across genes with KDE overlay
- Vertical line at piN/piS = 1
- Option for faceted piN/piS panels
- Input: `gene_results.tsv`

### Boxplot (`pie plot boxplot`)
- Grouped boxplots of piN, piS, and piN/piS per chromosome
- Faceted into 3 panels (one per metric) since scales differ
- Horizontal neutral line at piN/piS = 1 on the piN/piS facet
- Input: `gene_results.tsv`

### Sliding Window (`pie plot sliding-window`)
- piN/piS line plot along genomic position, faceted by chromosome
- Horizontal neutral line at 1
- Input: `window_results.tsv`

## CLI Interface

`pie plot` becomes a `click.Group` with subcommands.

### Shared options (all subcommands)
- `-i / --input` (required): input TSV file
- `-o / --output` (required): output path (PNG/PDF/SVG auto-detected from extension)
- `-W / --width`: figure width in inches (default: 12)
- `-H / --height`: figure height in inches (default: 6)
- `--dpi`: resolution (default: 300)

### manhattan-only options
- `--log-scale`: log2 transform y-axis
- `--label-top N`: label top N outlier genes
- `--highlight-genes gene1,gene2`: label specific genes

### scatter-only options
- `--color-by-chrom`: color points by chromosome

## Implementation Architecture

Single `plot.py` module with plotnine. Structure:

- `_OKABE_ITO`: colorblind-safe palette constant
- `_base_theme(width, height, dpi)`: publication defaults (theme_bw + font sizes + figure dimensions)
- `_load_gene_results(path)`: read TSV, drop NA piN_piS, add cumulative position columns
- `_sort_chroms(chroms)`: sort chromosomes numerically then alphabetically
- One function per plot type: `manhattan_plot()`, `scatter_plot()`, `histogram_plot()`, `boxplot_plot()`, `sliding_window_plot()`

## Dependencies

- Add `plotnine` to `environment.yml`
- Remove direct matplotlib usage (plotnine wraps it)

## Testing

- Smoke test per plot type: runs without error, produces output file
- Edge cases: empty dataframe, single chromosome, single gene
- CLI subcommands via CliRunner: exit code 0, output file exists
- Output format detection: `.png`, `.pdf`, `.svg`
- No visual regression testing
