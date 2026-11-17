#!/bin/bash
# Submit fit -> validate -> run for one index, chained so each step starts
# only if the previous one succeeded (validation exits with an error if any
# check fails, so the run step will not start).
#
# Usage:  bash submit-index.sh leaf
#
# If a step fails, the later steps are cancelled automatically.

set -e
INDEX=${1:?Usage: bash submit-index.sh <index>   (leaf, bloom, frz or damage)}
mkdir -p logs output

# --parsable prints "jobid" or "jobid;cluster"; keep only the job ID
fit=$(sbatch --parsable -J si-fit-$INDEX si-1-fit.slurm $INDEX); fit=${fit%%;*}
val=$(sbatch --parsable -J si-validate-$INDEX \
        --dependency=afterok:$fit --kill-on-invalid-dep=yes \
        si-2-validate.slurm $INDEX); val=${val%%;*}
run=$(sbatch --parsable -J si-run-$INDEX \
        --dependency=afterok:$val --kill-on-invalid-dep=yes \
        si-3-run.slurm $INDEX); run=${run%%;*}

echo "Submitted $INDEX: fit $fit -> validate $val -> run $run"
echo "Logs: logs/si-fit-$INDEX-$fit.out, logs/si-validate-$INDEX-$val.out, logs/si-run-$INDEX-$run.out"
