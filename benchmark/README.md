# Benchmark: pie vs SNPGenie

Performance comparison between pie and [SNPGenie](https://github.com/chasewnelson/snpgenie) on the same pool-seq dataset.

## Dataset

*Acyrthosiphon pisum* (pea aphid) pool-seq data from `data/Acyrthosiphon_pisum/`:

| | |
|---|---|
| Chromosomes | 4 (X, 1, 2, 3) |
| Genes | 439 protein-coding |
| SNPs | 91,787 |
| Reference size | ~9.9 Mb |
| VCF source | freebayes `--pooled-continuous` |

## Results

3 replicates per condition, run on a SLURM node with 8 CPUs.

### pie vs SNPGenie

| Tool | Threads | Wall time (s) | Speedup |
|------|--------:|-------------:|--------:|
| SNPGenie | 1 | 184.8 ± 5.7 | ref |
| pie | 1 | 15.8 ± 2.0 | **11.7x** |
| pie | 2 | 14.6 ± 0.0 | **12.7x** |
| pie | 4 | 14.3 ± 0.1 | **13.0x** |
| pie | 8 | 14.5 ± 0.1 | **12.8x** |

pie is **~12x faster** than SNPGenie in single-threaded mode.

### pie thread scaling

| Threads | Wall time (s) | Speedup | Efficiency |
|--------:|-------------:|--------:|-----------:|
| 1 | 15.8 | 1.00x | 100% |
| 2 | 14.6 | 1.08x | 54% |
| 4 | 14.3 | 1.11x | 28% |
| 8 | 14.5 | 1.09x | 14% |

Thread scaling is minimal on this dataset. The serial portion (VCF/GFF parsing, I/O) dominates at this scale (~400 genes). Larger genomes with tens of thousands of genes should benefit more from multi-threading.

## Notes

- SNPGenie requires single-sequence FASTA, so it was run separately per chromosome (total wall time summed).
- SNPGenie input was converted from GFF3 to GTF (CDS records only) using `prepare_snpgenie_input.py`.
- Both tools used `--minfreq 0.01`.

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
