# pie — piN/piS Estimator Design

**Date:** 2026-02-25
**Status:** Approved

## Overview

`pie` is a Python CLI tool for estimating nucleotide diversity metrics (piN, piS, piN/piS) from pooled sequencing data on large eukaryotic genomes (300+ Mb, multiple scaffolds/chromosomes). It reimplements the core Nei-Gojobori algorithm from SNPGenie in a high-performance Python stack.

## Requirements

- **Organism scope:** Large eukaryotic genomes (300+ Mb)
- **Primary VCF source:** Freebayes pool-seq mode (GATK nice-to-have)
- **Sample model:** Single-sample pooled VCF (multi-sample future consideration)
- **VCF assumption:** Biallelic, pre-decomposed (`bcftools norm -m-`)
- **Genetic code:** Standard only
- **Performance target:** 5-15 min for 300Mb genome (numpy), upgradeable to 1-3 min (Rust)

## Architecture

Hybrid approach: numpy-vectorized Python with a clean interface for a future Rust drop-in extension on the codon inner loop.

### Project Structure

```
pie/
├── pyproject.toml
├── src/
│   └── pie/
│       ├── __init__.py
│       ├── cli.py              # Click CLI: pie run / pie plot / pie summary
│       ├── vcf.py              # VCF parsing (cyvcf2), filtering, auto-indexing
│       ├── annotation.py       # GFF3/GTF parsing (gffutils), isoform selection
│       ├── reference.py        # FASTA handling (pysam.FastaFile), codon extraction
│       ├── codon.py            # Nei-Gojobori lookup tables, site/diff calculations
│       ├── diversity.py        # Core piN/piS engine: per-gene, per-window, genome-wide
│       ├── parallel.py         # Gene-level multiprocessing orchestration
│       ├── plot.py             # Manhattan plot, sliding window plots (matplotlib)
│       └── io.py               # TSV output writers
├── tests/
│   └── data/                   # Small test VCF, GFF, FASTA
└── docs/
    └── plans/
```

### Dependencies (conda env: `pie`)

- cyvcf2 (VCF parsing, htslib-backed)
- gffutils (GFF3/GTF parsing, SQLite-backed)
- pysam (FASTA indexing/access)
- numpy (vectorized codon computation)
- pandas (result aggregation)
- matplotlib (Manhattan plots)
- click (CLI framework)

## Data Flow

```
                    ┌─────────┐
                    │  FASTA  │
                    └────┬────┘
                         │
┌───────┐  filter   ┌────▼────┐  extract CDS    ┌──────────┐
│  VCF  │──────────►│  merge  │◄────────────────│  GFF/GTF │
└───────┘  (QUAL,   │  by gene│  (longest        └──────────┘
            depth,   └────┬────┘   isoform)
            freq,         │
            PASS)         ▼
                 ┌────────────────┐
                 │  Per-gene loop │ (parallel across genes)
                 │                │
                 │  1. Get CDS coords (exons, strand, phase)
                 │  2. Extract reference codons from FASTA
                 │  3. Overlay variants → allele freqs per codon position
                 │  4. Nei-Gojobori: count N_sites, S_sites per codon
                 │  5. Compute pairwise N_diffs, S_diffs (freq-weighted)
                 │  6. Accumulate per-codon → per-gene → per-window
                 └───────┬────────┘
                         │
            ┌────────────┼────────────┐
            ▼            ▼            ▼
     ┌────────────┐ ┌─────────┐ ┌──────────┐
     │ gene.tsv   │ │ win.tsv │ │summary.tsv│
     └────────────┘ └─────────┘ └──────────┘
```

### Processing Steps

1. **Parse GFF/GTF** → gffutils DB → select longest isoform per gene → extract ordered CDS intervals (strand, phase)
2. **Check/index VCF** → auto bgzip + tabix if needed
3. **Index FASTA** (pysam.FastaFile) for random access
4. **Per-gene processing** (parallelized via multiprocessing.Pool):
   - Fetch CDS exon coordinates, concatenate to coding sequence
   - Extract reference codons from FASTA (respecting strand/phase, joining split codons across exon boundaries)
   - Query VCF for variants overlapping CDS regions
   - Apply filters (min freq, min depth, QUAL, PASS)
   - Map each variant to its codon position within the coding sequence
   - For each codon: compute allele frequencies, calculate N/S sites and freq-weighted N/S diffs
   - Accumulate: per-codon, per-gene totals, bp-based sliding windows
5. **Aggregate** gene results → genome-wide summary
6. **Write** TSV outputs

## Core Algorithm: Nei-Gojobori for Pool-Seq

### Lookup Tables (precomputed at import, numpy arrays)

- **SITES[64, 3]:** For each of 64 codons × 3 positions, fractional N and S site counts.
  E.g., position 3 of AAA → 2/3 N, 1/3 S.
- **DIFFS[64, 64]:** For each codon pair, N_diffs and S_diffs for single-step substitutions.

### Per-Codon Calculation

```
For codon at positions [p1, p2, p3] in CDS:

1. Get allele frequencies at each position:
   p1: {A: 0.95, C: 0.05}    (variant site)
   p2: {G: 1.0}               (monomorphic)
   p3: {T: 0.80, A: 0.20}    (variant site)

2. Enumerate possible codons via itertools.product:
   AGT: 0.95 × 1.0 × 0.80 = 0.76
   AGA: 0.95 × 1.0 × 0.20 = 0.19
   CGT: 0.05 × 1.0 × 0.80 = 0.04
   CGA: 0.05 × 1.0 × 0.20 = 0.01

3. N_sites = Σ (freq_i × sites_N[codon_i])
   S_sites = Σ (freq_i × sites_S[codon_i])

4. N_diffs = Σ_i Σ_j>i (freq_i × freq_j × diffs_N[codon_i, codon_j])
   S_diffs = Σ_i Σ_j>i (freq_i × freq_j × diffs_S[codon_i, codon_j])
```

### Vectorization Strategy

- Build allele frequency array: `freqs` shape `(n_codons, 3, 4)` — 4 nucleotides per position
- Boolean mask for polymorphic codons:
  ```python
  is_poly = np.sum(freqs > 0, axis=2) > 1  # (n_codons, 3)
  poly_codons = np.any(is_poly, axis=1)     # (n_codons,)
  ```
- Batch site-count lookup for ALL codons (cheap array indexing)
- Expensive codon enumeration only on `poly_codons` subset
- For monomorphic codons: N_diffs = S_diffs = 0

### Edge Cases

**Stop codons:** After enumerating possible codons, check for stop codons (TAA, TAG, TGA). If present, set frequency to 0 and renormalize remaining codon frequencies. Log warning if stop codon frequency > 1%.

**Multi-allelic within codon:** Even with biallelic VCF records, two separate records at the same codon position can produce 3+ alleles. The generic `itertools.product` approach handles this naturally. Variant loading must merge multiple VCF records at the same position into a single allele frequency vector.

**Exon-spanning codons:** When a codon is split across exon boundaries, concatenate exon sequences in CDS order before extracting codons. The phase/frame field in GFF resolves reading frame offset.

### Aggregation

```
gene_piN = Σ(codon N_diffs) / Σ(codon N_sites)
gene_piS = Σ(codon S_diffs) / Σ(codon S_sites)
gene_piN_piS = gene_piN / gene_piS  (NA if piS == 0)
```

## Parallelization

- **Model:** Gene-level parallelism via `multiprocessing.Pool`
- **Worker isolation:** Each worker opens its own cyvcf2.VCF and pysam.FastaFile handles
- **Why gene-level:** Independent units, naturally balanced, no shared state
- **Estimated resources (300Mb, ~30k genes, 8 threads):** ~2-4 GB peak memory, ~5-15 min

## VCF Handling

- Auto-detect VCF state on input:
  - Plain `.vcf` → bgzip + tabix, inform user
  - `.vcf.gz` without `.tbi` → tabix index, inform user
  - `.vcf.gz` + `.tbi` → proceed
- Allele frequency extraction: AD (allele depth) and DP fields
- Filters: min allele frequency, min depth, min QUAL, PASS-only flag

## Isoform Selection

For each gene, select the longest isoform (by total CDS length) via gffutils parent-child relationships. All other isoforms ignored.

## Output Specification

### gene_results.tsv

| Column | Description |
|--------|-------------|
| chrom | Chromosome/scaffold |
| gene_id | Gene identifier |
| transcript_id | Selected (longest) isoform |
| start | Gene start (genomic) |
| end | Gene end (genomic) |
| strand | +/- |
| n_codons | Number of codons analyzed |
| n_poly_codons | Codons with polymorphic sites |
| N_sites | Total nonsynonymous sites |
| S_sites | Total synonymous sites |
| N_diffs | Freq-weighted nonsynonymous differences |
| S_diffs | Freq-weighted synonymous differences |
| piN | N_diffs / N_sites |
| piS | S_diffs / S_sites |
| piN_piS | piN / piS (NA if piS == 0) |
| mean_depth | Mean coverage across gene's CDS |
| n_variants | Number of SNPs in CDS |

### window_results.tsv

| Column | Description |
|--------|-------------|
| chrom | Chromosome/scaffold |
| win_start | Window start (genomic bp) |
| win_end | Window end (genomic bp) |
| gene_id | Gene this window belongs to |
| n_codons | Codons in this window |
| N_sites | Nonsynonymous sites in window |
| S_sites | Synonymous sites in window |
| N_diffs | Nonsynonymous differences |
| S_diffs | Synonymous differences |
| piN | Window piN |
| piS | Window piS |
| piN_piS | Window piN/piS |

### summary.tsv

- Total genes analyzed, total codons, total variants in CDS
- Genome-wide piN, piS, piN/piS (summed across all genes)
- Mean/median per-gene piN, piS, piN/piS

### Manhattan plot (piN_piS_manhattan.png)

- X-axis: genomic position (chromosomes colored alternately)
- Y-axis: per-gene piN/piS
- Horizontal line at piN/piS = 1 (neutral expectation)

## CLI Interface

```bash
# Main analysis
pie run \
  --vcf variants.vcf.gz \
  --gff genes.gff3 \
  --fasta reference.fa \
  --outdir results/ \
  --min-freq 0.01 \
  --min-depth 10 \
  --min-qual 20 \
  --pass-only \
  --window-size 1000 \
  --window-step 100 \
  --threads 8

# Plot from existing results
pie plot \
  --gene-results results/gene_results.tsv \
  --output results/piN_piS_manhattan.png \
  --width 16 --height 6

# Print summary
pie summary results/summary.tsv
```

### Default Values

| Flag | Default |
|------|---------|
| --min-freq | 0.01 |
| --min-depth | 10 |
| --min-qual | 20 |
| --pass-only | off |
| --window-size | 1000 |
| --window-step | 100 |
| --threads | 1 |

## Future Considerations

- Multi-sample VCF support (per-population piN/piS)
- Rust extension for codon inner loop (drop-in replacement)
- GATK VCF format support
- Alternative genetic codes
