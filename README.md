# pie -- piN/piS Estimator for pool-seq data

<p align="center">
  <img src="logo.png" alt="pie logo" width="400" />
</p>

Estimate nonsynonymous (piN) and synonymous (piS) nucleotide diversity
from pooled sequencing data on large eukaryotic genomes. Reimplements the
Nei-Gojobori method from SNPGenie in a numpy-vectorized Python stack with
gene-level multiprocessing.

## Comparison with SNPGenie

`pie` reimplements the Nei-Gojobori pool-seq method from
[SNPGenie](https://github.com/chasewnelson/SNPGenie) (Nelson et al., 2015)
in a numpy-vectorized Python stack with gene-level multiprocessing.

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
multiallelic sites. On non-multiallelic sites the results match SNPGenie
with r = 1.0.

## Installation

```bash
# Create conda environment
mamba env create -f environment.yml

# Activate and install
mamba activate pie
pip install -e .
```

Requires Python >= 3.12. Key dependencies: cyvcf2, gffutils, pysam, numpy,
pandas, matplotlib, click.

## Quick start

```bash
# Run analysis
pie run \
  --vcf variants.vcf.gz \
  --gff genes.gff3 \
  --fasta reference.fa \
  --outdir results/ \
  --threads 8

# View summary
pie summary results/summary.tsv

# Generate Manhattan plot
pie plot \
  --gene-results results/gene_results.tsv \
  --output results/piN_piS_manhattan.png
```

## Input requirements

| Input | Format | Notes |
|-------|--------|-------|
| VCF | `.vcf` or `.vcf.gz` | Freebayes pool-seq mode. Single- or multi-sample; use `--sample` for multi-sample VCFs. Multiallelic sites are skipped by default; use `--keep-multiallelic` to merge them instead. Accepts both decomposed (`bcftools norm -m-`) and non-decomposed VCFs. Plain VCF is auto-bgzipped; missing `.tbi` index is auto-created. |
| Annotation | GFF3 or GTF | Must contain gene, mRNA/transcript, and CDS features. Format is auto-detected. |
| Reference | FASTA | Must match the VCF reference. Auto-indexed by pysam if `.fai` is missing. |

Allele frequencies are extracted from the AD (allele depth) format field.
Falls back to INFO AF/DP if AD is unavailable.

## CLI reference

### `pie run`

Run the piN/piS analysis.

```
pie run --vcf FILE --gff FILE --fasta FILE --outdir DIR [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--vcf` | required | Input VCF file (bgzipped or plain) |
| `--gff` | required | GFF3 or GTF annotation file |
| `--fasta` | required | Reference FASTA file |
| `--outdir` | required | Output directory (created if absent) |
| `--min-freq` | 0.01 | Minimum alt allele frequency |
| `--min-depth` | 10 | Minimum read depth at variant site |
| `--min-qual` | 20.0 | Minimum variant QUAL score |
| `--pass-only` | off | Only use PASS-filtered variants |
| `--keep-multiallelic` | off | Keep and merge multiallelic sites instead of skipping them |
| `--include-stop-codons` | off | Count stop-gained mutations as nonsynonymous (by default they are excluded, matching NG86/SNPGenie) |
| `--window-size` | 1000 | Sliding window size in bp |
| `--window-step` | 100 | Sliding window step in bp |
| `--threads` | 1 | Number of parallel worker processes |
| `--sample` | auto | Sample name to analyse (required for multi-sample VCFs) |

### `pie plot`

Generate a Manhattan plot from gene results.

```
pie plot --gene-results FILE --output FILE [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--gene-results` | required | Path to `gene_results.tsv` |
| `--output` | required | Output PNG path |
| `--width` | 16.0 | Figure width in inches |
| `--height` | 6.0 | Figure height in inches |

### `pie summary`

Print genome-wide summary statistics to the terminal.

```
pie summary SUMMARY_FILE
```

Takes the path to `summary.tsv` as a positional argument.

## Multi-sample VCFs

By default, `pie` expects a single-sample VCF. If the VCF contains two or
more samples, `pie` will abort and list the available sample names. Use
`--sample` to select one:

```bash
pie run --vcf multi.vcf.gz --sample pool_A --gff genes.gff3 --fasta ref.fa --outdir results/
```

When `--sample` is given, output files are prefixed with the sample name
(e.g., `pool_A.gene_results.tsv`, `pool_A.window_results.tsv`,
`pool_A.summary.tsv`). Single-sample VCFs do not require `--sample` and
produce unprefixed filenames as before.

## Output files

All outputs are written to `--outdir`. When `--sample` is specified, each
filename is prefixed with `{sample}.` (e.g., `pool_A.gene_results.tsv`).

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
| mean_variant_depth | Mean read depth across variant sites (0 when no variants) |
| n_variants | SNP count within CDS |

### `window_results.tsv`

Sliding windows within each gene's CDS. Columns: chrom, win_start,
win_end, gene_id, n_codons, N_sites, S_sites, N_diffs, S_diffs, piN,
piS, piN_piS.

### `summary.tsv`

Single-row genome-wide summary: total_genes, total_codons,
total_variants, genome_piN, genome_piS, genome_piN_piS,
mean_gene_piN, mean_gene_piS, median_gene_piN, median_gene_piS.

Genome-wide piN and piS are computed by summing N/S diffs and sites
across all genes before dividing (not averaging per-gene ratios).

### `piN_piS_manhattan.png`

Manhattan plot with genomic position on the x-axis, per-gene piN/piS on
the y-axis, chromosomes in alternating colors, and a dashed line at
piN/piS = 1 (neutral expectation).

## Algorithm

`pie` uses the Nei-Gojobori (1986) method adapted for pool-seq allele
frequencies.

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

### Handling pool-seq frequencies

At each codon, allele frequencies from the VCF are used to enumerate all
possible codon haplotypes via the product of per-position allele
frequencies. For a codon with variants at positions 1 and 3:

```
p1: {A: 0.95, C: 0.05}   p2: {G: 1.0}   p3: {T: 0.80, A: 0.20}
  -> AGT (0.76), AGA (0.19), CGT (0.04), CGA (0.01)
```

This assumes linkage equilibrium between positions within a codon, which
is inherent to pool-seq data where individual haplotypes are not observed.

Site counts are the frequency-weighted sum across codon haplotypes.
Pairwise differences are computed for all haplotype pairs, weighted by
`2 * freq_i * freq_j`.

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

The regression tests pin piN/piS values for 400 genes validated against
SNPGenie, covering site counting, multiallelic VCF handling, and both
strand orientations.

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
