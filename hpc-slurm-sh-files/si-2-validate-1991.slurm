#!/bin/bash
### Step 2: validation checks on the saved fit
### Usage: sbatch si-2-validate-1991.slurm <index>     (leaf, bloom, frz or damage; default bloom)
### Normally submitted by submit-index-1991.sh, which chains steps 1-3.

#SBATCH --job-name=si-validate-1991
#SBATCH --account=klprudic
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
## Single-threaded; CPUs are requested only for memory (5 GB each on Puma).
## 8 CPUs = 40 GB. Earlier bloom run (70+ years) peaked at 25 GB, took 7 min;
## this model has about half the random effects, so expect less.
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ezylstra@arizona.edu

INDEX=${1:-bloom}

module load gnu13
module load gdal/3.8.5
module load R/4.5.2

# Use OpenBLAS in place of R's built-in BLAS (~14x faster for these models).
# One thread: more gave no gain. si-setup-1991.R logs a speed test to confirm it.
export LD_PRELOAD=/opt/ohpc/pub/libs/gnu13/openblas/0.3.21/lib/libopenblas.so
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

Rscript si-2-validate-1991.R $INDEX
