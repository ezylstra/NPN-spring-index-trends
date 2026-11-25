#!/bin/bash
### Step 1: fit the model for one index
### Usage: sbatch si-1-fit-1991.slurm <index>     (leaf, bloom, frz or damage; default bloom)
### Normally submitted by submit-index-1991.sh, which chains steps 1-3.

#SBATCH --job-name=si-fit-1991
#SBATCH --account=klprudic
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
## Single-threaded; CPUs are requested only for memory (5 GB each on Puma).
## 6 CPUs = 30 GB. Earlier bloom run (70+ years) peaked at 24 GB, took 8.1 h;
## this model has about half the random effects, so expect less.
#SBATCH --cpus-per-task=6
#SBATCH --time=24:00:00
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

Rscript si-1-fit-1991.R $INDEX
