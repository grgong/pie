# pie -- piN/piS Estimator

<p align="center">
  <img src="logo.png" alt="pie logo" width="400" />
</p>

Estimate nonsynonymous (piN) and synonymous (piS) nucleotide diversity
from **pool-seq** or **individual-sequencing** data on large eukaryotic
genomes.

Two analysis modes:

- **Pool mode** (default): Derives allele frequencies from AD (allele
  depth) fields in pool-seq VCFs.
- **Individual mode**: Derives allele frequencies from GT (genotype)
  fields across multiple diploid samples — a site-frequency based
  approximation consistent with pairwise diversity under the model
  assumptions.

## Comparison with SNPGenie

`pie` reimplements the Nei-Gojobori pool-seq method from
[SNPGenie](https://github.com/chasewnelson/SNPGenie) (Nelson et al., 2015)
in a numpy-vectorized Python stack with gene-level multiprocessing, and
extends it with an individual-sequencing mode.

**Stop codon site counting.**
SNPGenie excludes stop-producing mutations from site counting: for a
codon position where one of three possible changes introduces a stop
codon, only the two remaining sense→sense changes are counted (e.g.,
AGA pos 1: N = 1/2, S = 1/2). `pie` follows this same convention by
default, consistent with the classical NG86 method. Use
`--include-stop-codons` to instead count stop-gained mutations as
nonsynonymous (AGA pos 1: N = 2/3, S = 1/3), following the convention
used by tools such as PAML and libsequence.

**Multiallelic sites and frequency overlay.**
When a VCF is decomposed with `bcftools norm -m-`, a multiallelic site
(e.g., ALT=C,G) becomes two separate records at the same position.
SNPGenie processes each record independently with an incremental
subtract-from-ref / add-to-alt strategy, which can cause the reference
allele frequency to be double-subtracted. `pie` skips multiallelic sites
by default. With `--keep-multiallelic`, it groups records by genomic
position and recomputes frequencies from raw allele depths
(`freq(alt_i) = AD_i / (AD_ref + ΣAD_alt)`). For all sites, alt
frequencies are assigned directly and the reference is set to `1 − Σalt`.

**Simplified annotation input.**
SNPGenie requires a specific GTF format with a single isoform per gene,
which typically involves manual format conversion and isoform filtering.
`pie` accepts GFF3 or GTF directly (auto-detected) and automatically
selects the longest isoform per gene — no preprocessing needed.

**Auto-indexing.**
Plain VCF files are automatically bgzipped and tabix-indexed; FASTA files
are indexed via pysam if the `.fai` is missing. No manual preprocessing
required.

**Performance.**
Core computations are numpy-vectorized and genes are processed in parallel
via `multiprocessing.Pool`.

**Validation.**
A regression test suite pins piN/piS values for 400 real
*Acyrthosiphon pisum* genes, covering both strand orientations and
multiallelic sites. With unified inputs, synonymous/nonsynonymous site
counts match SNPGenie exactly (r = 1.0) and piN/piS values agree within
~1%, fully explained by the Nei (1987) sample-size correction that
SNPGenie applies (π = n/(n−1)·(1−Σfᵢ²)) vs. the population parameter
pie uses (π = 2p(1−p)). See [benchmark/](benchmark/) for details.

## Installation

```bash
git clone https://github.com/grgong/pie.git
cd pie

# Create conda environment
mamba env create -f environment.yml

# Activate and install
mamba activate pie
pip install -e .
```

Requires Python >= 3.12. Key dependencies: cyvcf2, gffutils, pysam, numpy,
pandas, plotnine, click.

## Quick start

```bash
# Pool-seq mode (default)
pie run -v variants.vcf.gz -g genes.gff3 -f reference.fa -o results/ -t 8

# Pool-seq with stricter filters
pie run -v variants.vcf.gz -g genes.gff3 -f reference.fa -o results/ \
  -d 20 -q 30 --pass-only -t 8

# Individual-sequencing mode (all samples)
pie run -v multi_sample.vcf.gz -g genes.gff3 -f reference.fa -o results/ \
  -m individual -t 8

# Individual mode with selected samples
pie run -v multi_sample.vcf.gz -g genes.gff3 -f reference.fa -o results/ \
  -m ind -S S1,S2,S3

# Individual mode with sample list file
pie run -v multi_sample.vcf.gz -g genes.gff3 -f reference.fa -o results/ \
  -m ind --samples-file samples.txt --min-call-rate 0.9

# View summary
pie summary results/summary.tsv

# Generate plots
pie plot manhattan -i results/gene_results.tsv -o manhattan.png
pie plot scatter -i results/gene_results.tsv -o scatter.png --color-by-chrom
pie plot histogram -i results/gene_results.tsv -o histogram.png
pie plot boxplot -i results/gene_results.tsv -o boxplot.png
pie plot sliding-window -i results/window_results.tsv -o sw.png
```

## Input requirements

| Input | Format | Notes |
|-------|--------|-------|
| VCF | `.vcf` or `.vcf.gz` | Single- or multi-sample. Multiallelic sites are skipped by default; use `--keep-multiallelic` to merge them. Accepts both decomposed (`bcftools norm -m-`) and non-decomposed VCFs. Plain VCF is auto-bgzipped; missing `.tbi` index is auto-created. See [Multi-sample VCFs](#multi-sample-vcfs) for sample selection. |
| Annotation | GFF3 or GTF | Must contain gene, mRNA/transcript, and CDS features. Format is auto-detected. |
| Reference | FASTA | Must match the VCF reference. Auto-indexed by pysam if `.fai` is missing. |

In pool mode, allele frequencies are extracted from the AD (allele depth)
format field, falling back to INFO AF/DP if AD is unavailable. In individual
mode, allele frequencies are derived from GT (genotype) fields:
`freq = alt_allele_count / (2 × called_samples)`.

## CLI reference

All commands support `-h`/`--help`.  Use `pie -V` to print the version.

### `pie run`

Run the piN/piS analysis.

```
pie run -v FILE -g FILE -f FILE -o DIR [OPTIONS]
```

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--vcf` | `-v` | required | Input VCF file (bgzipped or plain) |
| `--gff` | `-g` | required | GFF3 or GTF annotation file |
| `--fasta` | `-f` | required | Reference FASTA file (indexed) |
| `--outdir` | `-o` | required | Output directory (created if absent) |
| `--mode` | `-m` | pool | Analysis mode: `pool` or `individual` (alias `ind`) |
| `--min-freq` | | 0.01 | Minimum alt allele frequency |
| `--min-depth` | `-d` | 10 | Minimum read depth at variant site (pool mode only) |
| `--min-qual` | `-q` | 20.0 | Minimum variant QUAL score |
| `--pass-only` | | off | Only use PASS-filtered variants |
| `--keep-multiallelic` | | off | Keep and merge multiallelic sites instead of skipping them |
| `--include-stop-codons` | | off | Count stop-gained mutations as nonsynonymous (by default they are excluded, matching NG86/SNPGenie) |
| `--window-size` | `-w` | 1000 | Sliding window size in bp |
| `--window-step` | `-W` | 100 | Sliding window step in bp |
| `--threads` | `-t` | 1 | Number of parallel worker processes |
| `--sample` | `-s` | — | Sample name to analyse (pool mode only; required for multi-sample VCFs) |
| `--samples` | `-S` | all | Comma-separated sample names (individual mode only; defaults to all VCF samples) |
| `--samples-file` | | — | File with one sample name per line (individual mode only; mutually exclusive with `--samples`) |
| `--min-call-rate` | | 0.8 | Minimum genotype call rate per site (individual mode only; range 0–1) |
| `--min-an` | | 2 | Minimum allele number (AN) per site (individual mode only) |

### `pie plot`

Create publication-ready plots from piN/piS results.  Five subcommands are
available; all produce 300 DPI output in PNG, PDF, or SVG (auto-detected from
the file extension).  Default figure sizes follow NPG (Nature) journal specs.

**Shared options** (all subcommands except `sliding-window`):

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--input` | `-i` | required | Input TSV file |
| `--output` | `-o` | required | Output plot path |
| `--width` | `-W` | *(per subcommand)* | Figure width in inches |
| `--height` | `-H` | *(per subcommand)* | Figure height in inches |
| `--dpi` | | 300 | Resolution in dots per inch |
| `--min-codons` | | — | Exclude genes with fewer than N codons |
| `--min-variants` | | — | Exclude genes with fewer than N variants |

#### `pie plot manhattan`

Genome-wide Manhattan plot of per-gene piN/piS.  Input: `gene_results.tsv`.
Default size: 7.2 × 3.5 in (double column).

| Option | Default | Description |
|--------|---------|-------------|
| `--log-scale` | off | Use log₂ scale for piN/piS y-axis |
| `--label-top` | — | Label top N outlier genes |
| `--highlight-genes` | — | Comma-separated gene IDs to label |
| `--max-ratio` | 2.0 | Exclude genes with piN/piS above this value |
| `--exclude-zero-ratio` | off | Exclude genes with piN/piS = 0 |

#### `pie plot scatter`

piN vs piS scatter plot with point size reflecting gene length.
Input: `gene_results.tsv`.  Default size: 4.7 × 4.7 in (1.5 column).

| Option | Default | Description |
|--------|---------|-------------|
| `--color-by-chrom` | off | Color points by chromosome |
| `--max-piN` | 2.0 | Exclude genes with piN above this value |
| `--max-piS` | 2.0 | Exclude genes with piS above this value |

#### `pie plot histogram`

Distribution of piN/piS ratios with density overlay.
Input: `gene_results.tsv`.  Default size: 3.5 × 3.0 in (single column).

| Option | Default | Description |
|--------|---------|-------------|
| `--max-ratio` | 2.0 | Exclude genes with piN/piS above this value |
| `--exclude-zero-ratio` | off | Exclude genes with piN/piS = 0 |

#### `pie plot boxplot`

Per-chromosome boxplots of piN, piS, and piN/piS (three faceted panels).
Each metric is filtered independently.  Input: `gene_results.tsv`.
Default size: 3.5 × 7.0 in (single column).

| Option | Default | Description |
|--------|---------|-------------|
| `--max-piN` | 2.0 | Exclude genes with piN above this value |
| `--max-piS` | 2.0 | Exclude genes with piS above this value |
| `--max-ratio` | 2.0 | Exclude genes with piN/piS above this value |
| `--exclude-zero-ratio` | off | Exclude genes with piN/piS = 0 |

#### `pie plot sliding-window`

Sliding window piN/piS line plot, one row per chromosome.
Input: `window_results.tsv`.  Default size: 7.2 in wide × 1.5 in per chromosome.

| Option | Default | Description |
|--------|---------|-------------|
| `--max-ratio` | 2.0 | Exclude windows with piN/piS above this value |

> **Note:** `--min-codons` and `--min-variants` are not available for
> `sliding-window` because it operates on window-level data, not per-gene data.

### `pie summary`

Print genome-wide summary statistics to the terminal.

```
pie summary SUMMARY_FILE
```

Takes the path to `summary.tsv` as a positional argument.

## Multi-sample VCFs

### Pool mode

By default, `pie` expects a single-sample VCF. If the VCF contains two or
more samples, `pie` will abort and list the available sample names. Use
`--sample` to select one:

```bash
pie run -v multi.vcf.gz -g genes.gff3 -f ref.fa -o results/ -s pool_A
```

When `--sample` is given, output files are prefixed with the sample name
(e.g., `pool_A.gene_results.tsv`, `pool_A.window_results.tsv`,
`pool_A.summary.tsv`). Single-sample VCFs do not require `--sample` and
produce unprefixed filenames as before.

### Individual mode

In individual mode (`--mode individual`), `pie` uses all samples in the VCF
by default. Use `--samples` (comma-separated) or `--samples-file` (one name
per line) to select a subset. See the [Quick start](#quick-start) for
examples.

Allele frequencies are computed per site as `AC / AN`, where AC is the
alternate allele count and AN = 2 × called_samples (diploid). Sites with
call rate below `--min-call-rate` or AN below `--min-an` are skipped.

## Output files

All outputs are written to `--outdir`. In pool mode, when `--sample` is
specified, each filename is prefixed with `{sample}.`
(e.g., `pool_A.gene_results.tsv`).

### `gene_results.tsv`

One row per gene. Columns:

| Column | Description |
|--------|-------------|
| chrom | Chromosome or scaffold |
| gene_id | Gene identifier from annotation |
| transcript_id | Selected longest isoform |
| start, end | Gene genomic coordinates |
| strand | + or - |
| n_codons | Codons analyzed (excluding stop codons) |
| n_poly_codons | Codons containing polymorphic sites |
| N_sites | Total nonsynonymous sites |
| S_sites | Total synonymous sites |
| N_diffs | Frequency-weighted nonsynonymous differences |
| S_diffs | Frequency-weighted synonymous differences |
| piN | N_diffs / N_sites |
| piS | S_diffs / S_sites |
| piN_piS | piN / piS (NA when piS = 0) |
| mean_variant_depth | Pool mode: mean read depth; individual mode: mean AN across variant sites (0 when no variants) |
| n_variants | SNP count within CDS |
| n_samples | Number of selected samples (individual mode only) |
| mean_call_rate | Mean genotype call rate across variant sites (individual mode only) |

### `window_results.tsv`

Sliding windows within each gene's CDS. Columns: chrom, win_start,
win_end, gene_id, n_codons, N_sites, S_sites, N_diffs, S_diffs, piN,
piS, piN_piS.

### `summary.tsv`

Single-row genome-wide summary: total_genes, total_codons,
cds_snp_variants, genome_piN, genome_piS, genome_piN_piS,
mean_gene_piN, mean_gene_piS, median_gene_piN, median_gene_piS.
In individual mode, two additional columns are included:
n_samples_selected and mean_call_rate (variant-site-weighted average).

Genome-wide piN and piS are computed by summing N/S diffs and sites
across all genes before dividing (not averaging per-gene ratios).

## Algorithm

`pie` uses the Nei-Gojobori (1986) method adapted for allele frequency
data from pool-seq or individual-sequencing experiments.

### Site counting

For each of the 64 codons, a precomputed lookup table stores fractional
nonsynonymous (N) and synonymous (S) site counts at each of the three
codon positions. Each position is evaluated by considering all three
possible single-nucleotide changes and classifying them as synonymous or
nonsynonymous under the standard genetic code. Mutations that introduce a
premature stop codon (stop-gained) are excluded from site counting by
default, matching the classical NG86/SNPGenie convention. For example,
if one of 3 changes produces a stop codon and the other is synonymous,
that position contributes N = 1/2 and S = 1/2 (counting only the two
sense→sense changes). Use `--include-stop-codons` to instead count
stop-gained mutations as nonsynonymous.

### Allele frequency handling

In pool mode, per-site allele frequencies are extracted from the AD
(allele depth) format field. In individual mode, frequencies are derived
from diploid genotypes across selected samples: `freq = AC / AN`, where
AC is the alternate allele count and AN = 2 × called_samples. Missing
genotypes (`./.`) are excluded; only called samples contribute to AC
and AN.

### Codon haplotype enumeration

Regardless of mode, per-site allele frequencies are used to enumerate
all possible codon haplotypes via the product of per-position
frequencies (linkage equilibrium assumption). For a codon with variants
at positions 1 and 3:

```
p1: {A: 0.95, C: 0.05}   p2: {G: 1.0}   p3: {T: 0.80, A: 0.20}
  -> AGT (0.76), AGA (0.19), CGT (0.04), CGA (0.01)
```

Site counts are the frequency-weighted sum across codon haplotypes.
Pairwise differences are computed for all haplotype pairs, weighted by
`2 * freq_i * freq_j`.

### Diversity estimator and sample-size correction

`pie` computes per-site diversity as the population parameter:

```
π = 2p(1−p)          (biallelic)
π = Σ_{i<j} 2·f_i·f_j  (multiallelic / codon-level)
```

SNPGenie instead applies the Nei (1987) unbiased estimator, which
includes a finite-sample correction:

```
π = n/(n−1) · (1 − Σf_i²)
```

where *n* is the read depth at the site. The ratio between the two
estimators is exactly (*n*−1)/*n*, producing a systematic difference of
~1% at 93× depth.

**Why `pie` does not apply this correction.** In pool-seq, two layers
of sampling exist:

1. **Pool composition** — a finite number of individuals (2*N*
   chromosomes) is drawn from the population.
2. **Sequencing** — reads are drawn from the pooled DNA at depth *n*.

The Nei (1987) correction was designed for *n* = number of sampled gene
copies (chromosomes), not sequencing reads. SNPGenie substitutes read
depth for chromosome count, which corrects only the sequencing-sampling
layer and conflates it with the biological-sampling layer. A complete
unbiased estimator would require correcting for both levels
(Futschik & Schlötterer 2010):

```
π̂ = 2·p̂·(1−p̂) · n/(n−1) · 2N/(2N−1)
```

Since the pool size is rarely known with precision and the correction is
negligible at typical pool-seq depths (>50×), `pie` reports the plug-in
estimator directly. This is transparent — no implicit assumption about
sample size — and has no effect on piN/piS ratios or between-group
comparisons, where the systematic bias cancels. For a quantitative
concordance analysis, see [benchmark/](benchmark/).

### Multi-step substitutions

For codon pairs differing at 2 or 3 positions, all shortest mutational
pathways are enumerated. N and S differences are averaged across valid
pathways (those not passing through a stop codon intermediate).

### Edge cases

- **Stop-gained mutations:** By default, SNPs that introduce a premature
  stop codon are excluded from site counting and pairwise difference
  calculations, matching the classical NG86 method and SNPGenie. With
  `--include-stop-codons`, stop-gained mutations are instead counted as
  nonsynonymous changes, consistent with the convention used by PAML and
  libsequence.
- **Exon-spanning codons:** CDS exons are concatenated in genomic order
  before codon extraction. Codons split across exon boundaries are handled
  by tracking per-base genomic positions through the concatenated sequence.
- **Isoform selection:** The longest isoform (by total CDS length) is
  selected per gene; all others are ignored.

### Aggregation

```
gene piN  = sum(codon N_diffs) / sum(codon N_sites)
gene piS  = sum(codon S_diffs) / sum(codon S_sites)
gene piN/piS = piN / piS   (NA if piS = 0)
```

## Parallelization

Genes are processed independently via `multiprocessing.Pool`. Each worker
opens its own VCF and FASTA file handles to avoid shared state.

## Testing

```bash
# Unit and integration tests (~1.4s)
pytest tests/ -m "not slow"

# Regression tests against real Acyrthosiphon pisum dataset (~12s)
# Requires: tar xzf data/Acyrthosiphon_pisum.tar.gz -C data/
pytest tests/ -m slow

# All tests
pytest tests/
```

The regression tests pin piN/piS values for 400 genes, covering site
counting, multiallelic VCF handling, and both strand orientations.

## Citation

This tool implements the averaging method for counting synonymous and
nonsynonymous substitutions described in:

> Nei M, Gojobori T. Simple methods for estimating the numbers of
> synonymous and nonsynonymous nucleotide substitutions. *Mol Biol Evol*.
> 1986;3(5):418-426.

The pool-seq adaptation (frequency-weighted site and difference counting)
follows the approach used in SNPGenie:

> Nelson CW, Moncla LH, Hughes AL. SNPGenie: estimating evolutionary
> parameters to detect natural selection using pooled next-generation
> sequencing data. *Bioinformatics*. 2015;31(22):3709-3711.
