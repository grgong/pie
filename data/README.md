# Test Data

## Acyrthosiphon pisum pool-seq dataset

Mini test dataset extracted from pea aphid (*Acyrthosiphon pisum*) pool-seq
data (sample SRR27175631). Contains 4 chromosomes (X, 1, 2, 3) with ~100
protein-coding genes each (400 genes total, ~92k SNPs).

### Setup

```bash
tar xzf data/Acyrthosiphon_pisum.tar.gz -C data/
```

### Contents

| File | Description |
|------|-------------|
| `Acyrthosiphon_pisum.fa` + `.fai` | Reference genome (coordinate-shifted subregions) |
| `Acyrthosiphon_pisum.gff` | GFF3 gene annotations (400 protein-coding genes) |
| `SRR27175631.filtered.snps.vcf.gz` + `.tbi` | Filtered SNP variants from pool-seq |

### Run pie on this dataset

```bash
pie run \
  --vcf data/Acyrthosiphon_pisum/SRR27175631.filtered.snps.vcf.gz \
  --gff data/Acyrthosiphon_pisum/Acyrthosiphon_pisum.gff \
  --fasta data/Acyrthosiphon_pisum/Acyrthosiphon_pisum.fa \
  --outdir results/Acyrthosiphon_pisum \
  --threads 4
```

### Source

Created by `scripts/create_real_testset.sh` from:
- VCF: `SRR27175631.filtered.snps.vcf.gz` (freebayes, filtered)
- GFF: NCBI *A. pisum* annotation (`Acyrthosiphon_pisum.rename.gff`)
- FASTA: NCBI *A. pisum* genome assembly

Coordinates are shifted so each chromosome subregion starts at position 1.
