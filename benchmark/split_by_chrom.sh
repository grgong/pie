#!/usr/bin/env bash
set -euo pipefail

# Split FASTA, VCF, and GTF by chromosome for SNPGenie
BENCH_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$(dirname "$BENCH_DIR")/data/Acyrthosiphon_pisum"
SNPGENIE_DIR="$BENCH_DIR/snpgenie_input"

mkdir -p "$SNPGENIE_DIR"

FASTA="$DATA_DIR/Acyrthosiphon_pisum.fa"
VCF="$BENCH_DIR/SRR27175631.filtered.snps.vcf"
GTF="$BENCH_DIR/Acyrthosiphon_pisum.gtf"

# Get chromosome names from FAI
chroms=($(cut -f1 "$DATA_DIR/Acyrthosiphon_pisum.fa.fai"))

for chrom in "${chroms[@]}"; do
    dir="$SNPGENIE_DIR/$chrom"
    mkdir -p "$dir"

    # Split FASTA: extract single chromosome
    samtools faidx "$FASTA" "$chrom" > "$dir/ref.fa"

    # Split VCF: header + matching lines
    grep "^#" "$VCF" > "$dir/variants.vcf"
    grep -P "^${chrom}\t" "$VCF" >> "$dir/variants.vcf" || true

    # Split GTF: matching lines
    grep -P "^${chrom}\t" "$GTF" > "$dir/genes.gtf" || true

    n_snps=$(grep -cv "^#" "$dir/variants.vcf" || echo 0)
    n_cds=$(wc -l < "$dir/genes.gtf")
    echo "  $chrom: $n_snps SNPs, $n_cds CDS features"
done

echo "Done. Per-chromosome input in $SNPGENIE_DIR/"
