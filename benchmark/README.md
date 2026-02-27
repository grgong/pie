# Benchmark: pie vs SNPGenie

Performance comparison between pie and [SNPGenie](https://github.com/chasewnelson/snpgenie) on the same pool-seq dataset.

## Dataset

*Acyrthosiphon pisum* (pea aphid) pool-seq data from `data/Acyrthosiphon_pisum/`:

| | |
|---|---|
| Chromosomes | 4 (X, 1, 2, 3) |
| Genes | 400 protein-coding (longest isoform per gene) |
| SNPs | 91,787 |
| Reference size | ~9.9 Mb |
| VCF source | freebayes `--pooled-continuous` |

## Results

3 replicates per condition, login node with 8 CPUs. "Cold" = annotation DB rebuilt each run; "warm" = cached `.pie.db` reused.

### pie vs SNPGenie

| Tool | Threads | Wall time (s) | Speedup |
|------|--------:|-------------:|--------:|
| SNPGenie | 1 | 180.2 ± 0.8 | ref |
| pie cold | 1 | 4.7 ± 0.9 | **38.7x** |
| pie cold | 2 | 3.7 ± 0.0 | **49.3x** |
| pie cold | 4 | 3.7 ± 0.0 | **48.4x** |
| pie cold | 8 | 3.7 ± 0.0 | **48.2x** |
| pie warm | 1 | 3.0 ± 0.0 | **60.8x** |
| pie warm | 2 | 2.4 ± 0.0 | **73.7x** |
| pie warm | 4 | 2.6 ± 0.0 | **70.2x** |
| pie warm | 8 | 2.6 ± 0.0 | **70.3x** |

pie is **~39x faster** than SNPGenie single-threaded (cold), or **~61x** with cached annotations (warm).

### Optimizations (vs previous pie version)

| Optimization | Effect |
|---|---|
| `gffutils` strategy: `merge` → `create_unique` with fallback | ~12x faster DB creation |
| Annotation DB cache (`.pie.db` with checksum validation) | near-instant reuse on repeated runs |
| Window output: O(n·m) → O(n log n) via prefix sums + bisect | faster `window_results.tsv` generation |
| Stop-codon warnings: per-codon → per-gene summary | reduced I/O noise |

### pie thread scaling (cold)

| Threads | Wall time (s) | Speedup | Efficiency |
|--------:|-------------:|--------:|-----------:|
| 1 | 4.7 | 1.00x | 100% |
| 2 | 3.7 | 1.28x | 64% |
| 4 | 3.7 | 1.25x | 31% |
| 8 | 3.7 | 1.25x | 16% |

Thread scaling is minimal on this dataset. The serial portion (VCF/GFF parsing, I/O) dominates at this scale (~400 genes). Larger genomes with tens of thousands of genes should benefit more from multi-threading.

## Concordance analysis

Detailed comparison of piN/piS values between pie (`--keep-multiallelic`) and SNPGenie, using unified input (same isoform selection, same GTF derived from pie's GFF3 parsing).

### Genome-wide agreement

| | piN | piS | piN/piS |
|---|---|---|---|
| pie | 0.00101382 | 0.00416695 | 0.243300 |
| SNPGenie | 0.00103617 | 0.00420721 | 0.246284 |
| Relative diff | **2.2%** | **1.0%** | **1.2%** |

### Per-gene agreement (210 matched genes)

| Metric | Pearson r | Spearman r |
|--------|:---------:|:----------:|
| N_sites / S_sites | 1.0000 | 1.0000 |
| piN | 0.996 | 0.993 |
| piS | 0.998\* | 0.999 |
| N_diffs / S_diffs | 0.996 | 0.993 |

\*Excluding one outlier gene (see below).

- Median per-gene relative difference: piN **0.97%**, piS **0.85%**
- 84% of genes have piN within 5%; 92% have piS within 5%
- pie reports slightly lower values due to the Nei correction (see below)

### Sources of difference

Three factors explain all observed differences between pie and SNPGenie:

#### 1. Multi-allelic site handling (resolved with `--keep-multiallelic`)

pie defaults to skipping positions with >1 ALT allele. The freebayes VCF contains 3.85% multi-allelic sites (mostly `ALT=X,*` with spanning deletions), affecting 4,634 variants in 198 genes. Using `--keep-multiallelic` recovers these sites and reduces the genome-wide piN difference from ~8% to ~1%.

#### 2. Nei (1987) sample-size correction (~1% systematic bias)

The remaining ~1% systematic difference is fully explained by a statistical choice:

| Tool | Per-site diversity formula |
|------|--------------------------|
| pie | π = 2p(1−p) (population parameter) |
| SNPGenie | π = n/(n−1) · (1 − Σfᵢ²) (Nei 1987 unbiased estimator) |

At depth *n*, the ratio is exactly (n−1)/n. At the dataset's mean variant depth of ~93×, the expected bias is 1/92 = **1.09%**, matching the observed ~1%. Verified to 6 decimal places on individual sites (e.g., depth=133: pie π=0.4016, SNPGenie π=0.4046, ratio=0.9925 = 132/133).

Both approaches are standard; the choice depends on whether reads are treated as exact observations or as a sample from an underlying population.

#### 3. SNPGenie bug with decomposed multi-allelic VCF records (one outlier gene)

Apisum_003668 shows a 3.5× piS discrepancy (pie=0.018, SNPGenie=0.063). The cause is a single codon (GGG at positions 111787–111789) where freebayes decomposed a 4-allele complex variant into 6 VCF records, all with 0 reference reads. SNPGenie cannot properly handle this input: it produces negative nucleotide counts, S_diffs=17.96 for a single codon (maximum possible is 1.0), and N_sites=16 (maximum for one codon is 3). pie correctly merges the decomposed records and computes the expected synonymous diversity. Excluding this gene, piS Pearson r > 0.999.

## Notes

- SNPGenie requires single-sequence FASTA, so it was run separately per chromosome (total wall time summed).
- `prepare_snpgenie_input.py` uses pie's annotation parser (`pie.annotation.parse_annotations`) to select the longest isoform per gene, ensuring both tools use the exact same gene models.
- Both tools used `--minfreq 0.01`.
- pie was run with `--keep-multiallelic` to include multi-allelic sites.
- Performance numbers include `--keep-multiallelic`.

## Reproduce

```bash
# Prepare SNPGenie-compatible input (split by chromosome, GFF3→GTF)
python benchmark/prepare_snpgenie_input.py \
  data/Acyrthosiphon_pisum/Acyrthosiphon_pisum.gff \
  benchmark/Acyrthosiphon_pisum.gtf
bash benchmark/split_by_chrom.sh

# Run on login node (1–2 threads)
BENCH_THREADS="1 2" REPS=3 bash benchmark/run_benchmark.sh

# Or submit to SLURM for more threads (1–8)
sbatch benchmark/run_benchmark.sbatch
```
