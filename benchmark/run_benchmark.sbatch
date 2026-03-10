#!/usr/bin/env bash
#SBATCH --job-name=pie-bench
#SBATCH --partition=workq
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=benchmark/slurm-%x-%j.out
#SBATCH --error=benchmark/slurm-%x-%j.err

set -euo pipefail

source /work/user/ggong/miniforge3/etc/profile.d/conda.sh
source /work/user/ggong/miniforge3/etc/profile.d/mamba.sh
mamba activate pie

export BENCH_THREADS="1 2 4 8"
export REPS=3

bash benchmark/run_benchmark.sh
