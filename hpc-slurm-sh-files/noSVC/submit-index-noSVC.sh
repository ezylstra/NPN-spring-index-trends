#!/bin/bash
# Submit fit -> validate -> run for one index, NO-SVC MODEL, chained so each
# step starts only if the previous one succeeded (the fit exits with an error
# if the Hessian is not positive definite, and validation exits with an error
# if any check fails).
#
# Usage:  bash submit-index-noSVC.sh leaf
#         bash submit-index-noSVC.sh leaf skip-fit   # fit file already exists:
#                                                    # submit validate -> run only
#
# If a step fails, the later steps are cancelled automatically.

set -e
INDEX=${1:?Usage: bash submit-index-noSVC.sh <index> [skip-fit]   (leaf, bloom, frz or damage)}
MODE=${2:-}
mkdir -p logs output

if [ "$MODE" = "skip-fit" ]; then
  FIT_FILE=output/si-$INDEX-1951-2024-northhemis-noSVC-fit.rds
  [ -f "$FIT_FILE" ] || { echo "Not found: $FIT_FILE"; exit 1; }
  val=$(sbatch --parsable -J si-validate-noSVC-$INDEX \
          si-2-validate-noSVC.slurm $INDEX); val=${val%%;*}
  fit="(skipped)"
else
  # --parsable prints "jobid" or "jobid;cluster"; keep only the job ID
  fit=$(sbatch --parsable -J si-fit-noSVC-$INDEX si-1-fit-noSVC.slurm $INDEX); fit=${fit%%;*}
  val=$(sbatch --parsable -J si-validate-noSVC-$INDEX \
          --dependency=afterok:$fit --kill-on-invalid-dep=yes \
          si-2-validate-noSVC.slurm $INDEX); val=${val%%;*}
fi
run=$(sbatch --parsable -J si-run-noSVC-$INDEX \
        --dependency=afterok:$val --kill-on-invalid-dep=yes \
        si-3-run-noSVC.slurm $INDEX); run=${run%%;*}

echo "Submitted $INDEX (no SVC): fit $fit -> validate $val -> run $run"
echo "Logs: logs/si-fit-noSVC-$INDEX-$fit.out, logs/si-validate-noSVC-$INDEX-$val.out, logs/si-run-noSVC-$INDEX-$run.out"
