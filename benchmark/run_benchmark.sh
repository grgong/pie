#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Benchmark: pie vs SNPGenie
# Dataset: Acyrthosiphon pisum (4 chromosomes, ~400 genes, ~92k SNPs)
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="$PROJECT_DIR/data/Acyrthosiphon_pisum"
SNPGENIE="$PROJECT_DIR/../SNPGenie/snpgenie.pl"
BENCH_DIR="$SCRIPT_DIR"

# Thread counts to benchmark for pie (adjust for available cores)
THREADS=(${BENCH_THREADS:-1 2})

# Number of replicates
REPS=${REPS:-3}

# Results file
RESULTS="$BENCH_DIR/benchmark_results.tsv"

# Input files
VCF="$DATA_DIR/SRR27175631.filtered.snps.vcf.gz"
GFF="$DATA_DIR/Acyrthosiphon_pisum.gff"
FASTA="$DATA_DIR/Acyrthosiphon_pisum.fa"
SNPGENIE_INPUT="$BENCH_DIR/snpgenie_input"

# Chromosomes
CHROMS=($(cut -f1 "$DATA_DIR/Acyrthosiphon_pisum.fa.fai"))
N_SNPS=$(zcat "$VCF" | grep -cv '^#')

echo "============================================================"
echo "Benchmark: pie vs SNPGenie"
echo "Dataset: Acyrthosiphon pisum"
echo "  Chromosomes: ${#CHROMS[@]}"
echo "  Genes: $(grep -cP '\tgene\t' "$GFF" || true)"
echo "  SNPs: $N_SNPS"
echo "  Replicates: $REPS"
echo "  pie threads: ${THREADS[*]}"
echo "============================================================"
echo ""

# Header
echo -e "tool\tthreads\trep\treal_sec\tuser_sec\tsys_sec" > "$RESULTS"

# --------------------------------------------------------------------------
# SNPGenie benchmark (run per-chromosome, sum total time)
# SNPGenie requires single-sequence FASTA and creates temp files in cwd,
# so we cd into a per-chromosome working directory for each run.
# --------------------------------------------------------------------------
echo "--- SNPGenie (single-threaded, per-chromosome) ---"

for rep in $(seq 1 "$REPS"); do
    total_real=0
    total_user=0
    total_sys=0

    for chrom in "${CHROMS[@]}"; do
        workdir="$BENCH_DIR/snpgenie_work/${chrom}_rep${rep}"
        rm -rf "$workdir"
        mkdir -p "$workdir"

        # Copy input to work directory (SNPGenie expects local files)
        cp "$SNPGENIE_INPUT/$chrom/ref.fa" "$workdir/"
        cp "$SNPGENIE_INPUT/$chrom/variants.vcf" "$workdir/"
        cp "$SNPGENIE_INPUT/$chrom/genes.gtf" "$workdir/"

        timefile=$(mktemp)
        (
            cd "$workdir"
            /usr/bin/time -f "%e\t%U\t%S" -o "$timefile" \
                perl "$SNPGENIE" \
                    --vcfformat=4 \
                    --snpreport=variants.vcf \
                    --fastafile=ref.fa \
                    --gtffile=genes.gtf \
                    --minfreq=0.01 \
                > /dev/null 2>&1
        )

        read -r real user sys < "$timefile"
        total_real=$(echo "$total_real + $real" | bc)
        total_user=$(echo "$total_user + $user" | bc)
        total_sys=$(echo "$total_sys + $sys" | bc)
        rm -f "$timefile"

        printf "    chrom=%-2s  %.1fs\n" "$chrom" "$real"
    done

    echo -e "SNPGenie\t1\t${rep}\t${total_real}\t${total_user}\t${total_sys}" >> "$RESULTS"
    printf "  %-12s rep=%d  total=%.1fs\n" "SNPGenie" "$rep" "$total_real"
done
echo ""

# --------------------------------------------------------------------------
# pie benchmark (varying threads)
# --------------------------------------------------------------------------
for t in "${THREADS[@]}"; do
    echo "--- pie (threads=$t) ---"
    for rep in $(seq 1 "$REPS"); do
        outdir="$BENCH_DIR/pie_t${t}_rep${rep}"
        rm -rf "$outdir"

        timefile=$(mktemp)
        /usr/bin/time -f "%e\t%U\t%S" -o "$timefile" \
            pie run \
                --vcf "$VCF" \
                --gff "$GFF" \
                --fasta "$FASTA" \
                --outdir "$outdir" \
                --threads "$t" \
                --min-freq 0.01 \
                --keep-multiallelic \
            > /dev/null 2>&1

        read -r real user sys < "$timefile"
        rm -f "$timefile"

        echo -e "pie\t${t}\t${rep}\t${real}\t${user}\t${sys}" >> "$RESULTS"
        printf "  %-12s rep=%d  %.1fs (user: %.1fs)\n" "pie -t $t" "$rep" "$real" "$user"
    done
    echo ""
done

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
echo "============================================================"
echo "Results summary (mean ± sd of $REPS replicates):"
echo "============================================================"
echo ""

python3 - "$RESULTS" <<'PYEOF'
import sys
from collections import defaultdict
import math

results = defaultdict(list)
with open(sys.argv[1]) as f:
    next(f)  # skip header
    for line in f:
        tool, threads, rep, real, user, sys_ = line.strip().split("\t")
        key = (tool, int(threads))
        results[key].append((float(real), float(user), float(sys_)))

def mean(xs): return sum(xs) / len(xs)
def sd(xs):
    m = mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / max(len(xs) - 1, 1)) if len(xs) > 1 else 0.0

# Print table
print(f"{'Tool':<12} {'Threads':>7} {'Real (s)':>14} {'User (s)':>14} {'Speedup':>8}")
print("-" * 60)

snpgenie_mean = None
pie_1t_mean = None

for key in sorted(results.keys(), key=lambda k: (0 if k[0] == "SNPGenie" else 1, k[1])):
    tool, threads = key
    reals = [r[0] for r in results[key]]
    users = [r[1] for r in results[key]]
    m_real, s_real = mean(reals), sd(reals)
    m_user, s_user = mean(users), sd(users)

    if tool == "SNPGenie":
        snpgenie_mean = m_real
        speedup_str = "ref"
    else:
        if pie_1t_mean is None:
            pie_1t_mean = m_real
        speedup_vs_snpgenie = snpgenie_mean / m_real if m_real > 0 else float("inf")
        speedup_str = f"{speedup_vs_snpgenie:.1f}x"

    print(f"{tool:<12} {threads:>7} {m_real:>7.1f} ± {s_real:<4.1f} {m_user:>7.1f} ± {s_user:<4.1f} {speedup_str:>8}")

# pie scaling efficiency
if pie_1t_mean and len([k for k in results if k[0] == "pie"]) > 1:
    print()
    print("pie thread scaling (vs pie -t 1):")
    for key in sorted(results.keys(), key=lambda k: k[1]):
        tool, threads = key
        if tool != "pie":
            continue
        m_real = mean([r[0] for r in results[key]])
        scaling = pie_1t_mean / m_real if m_real > 0 else 0
        efficiency = scaling / threads * 100
        print(f"  -t {threads}: {m_real:.1f}s  speedup={scaling:.2f}x  efficiency={efficiency:.0f}%")

print()
PYEOF

echo "Full results: $RESULTS"
echo "Done."
