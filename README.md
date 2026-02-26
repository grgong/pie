# pie -- piN/piS Estimator for pool-seq data

<p align="center">
  <img src="logo.png" alt="pie logo" width="220" />
</p>

Estimate nonsynonymous (piN) and synonymous (piS) nucleotide diversity
from pooled sequencing data on large eukaryotic genomes. Reimplements the
Nei-Gojobori method from SNPGenie in a numpy-vectorized Python stack with
gene-level multiprocessing.

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
| VCF | `.vcf` or `.vcf.gz` | Freebayes pool-seq mode. Multiallelic sites are skipped by default; use `--keep-multiallelic` to merge them instead. Accepts both decomposed (`bcftools norm -m-`) and non-decomposed VCFs. Plain VCF is auto-bgzipped; missing `.tbi` index is auto-created. |
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
| `--exclude-stop-codons` | off | Exclude stop-gained mutations from piN (by default they count as nonsynonymous) |
| `--window-size` | 1000 | Sliding window size in bp |
| `--window-step` | 100 | Sliding window step in bp |
| `--threads` | 1 | Number of parallel worker processes |

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

## Output files

All outputs are written to `--outdir`.

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
| mean_depth | Mean read depth across CDS variants |
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
premature stop codon (stop-gained) are counted as nonsynonymous by
default. For example, if 2 of 3 changes are nonsynonymous (including any
that produce a stop codon), that position contributes 2/3 N sites and 1/3
S sites. Use `--exclude-stop-codons` to exclude these mutations from site
counting instead.

### Handling pool-seq frequencies

At each codon, allele frequencies from the VCF are used to enumerate all
possible codon haplotypes via the product of per-position allele
frequencies. For a codon with variants at positions 1 and 3:

```
p1: {A: 0.95, C: 0.05}   p2: {G: 1.0}   p3: {T: 0.80, A: 0.20}
  -> AGT (0.76), AGA (0.19), CGT (0.04), CGA (0.01)
```

Site counts are the frequency-weighted sum across codon haplotypes.
Pairwise differences are computed for all haplotype pairs, weighted by
`2 * freq_i * freq_j`.

### Multi-step substitutions

For codon pairs differing at 2 or 3 positions, all shortest mutational
pathways are enumerated. N and S differences are averaged across valid
pathways (those not passing through a stop codon intermediate).

### Edge cases

- **Stop-gained mutations:** By default, SNPs that introduce a premature
  stop codon are treated as nonsynonymous changes, consistent with the
  standard approach in population genetics (e.g., PAML, libsequence).
  Both the site counting (N_sites) and pairwise difference (N_diffs)
  calculations include these mutations. With `--exclude-stop-codons`,
  stop-gained alleles are removed from polymorphic codons and excluded
  from site counting (legacy behavior).
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
opens its own VCF and FASTA file handles to avoid shared state. For a
300 Mb genome with ~30k genes on 8 threads, expect ~5-15 minutes and
~2-4 GB peak memory.

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
