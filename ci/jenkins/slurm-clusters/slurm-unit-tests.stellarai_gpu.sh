#!/bin/bash -l
# StellarAI GPU unit-test payload for the shared Slurm Pipeline.
#
# This script deliberately has no #SBATCH directives. The trusted Jenkins
# pipeline owns the resource request, while this script owns only the runtime
# environment and the test command.
set -euo pipefail

# Jenkins supplies the shared output directory and candidate checkout.
: "${CI_WORKSPACE:?CI_WORKSPACE must name the shared Jenkins workspace}"
: "${CI_CANDIDATE_DIR:?CI_CANDIDATE_DIR must name the candidate checkout}"

cd "$CI_CANDIDATE_DIR"

# Load the build environment explicitly in this batch shell.
. machines/module_load.stellarai-gpu.sh
# Anaconda is needed to download ADAS during builds, but interferes with Open MPI.
module unload anaconda3/2026.7
# These tests stay on one node; allow shared-memory and CUDA-copy MPI transports.
export UCX_TLS=sm,cuda_copy,self

# The job allocation supplies GPUs through GRES; srun inherits that allocation.
# Run the unit suite and write per-test results for the Jenkins report.
started="$(date +%s)"
GKYL_UNIT_RESULTS="$CI_WORKSPACE/candidate-unit-results.txt" srun --ntasks=1 --cpus-per-task="${CI_SLURM_CPUS_PER_TASK:-4}" --cpu-bind=cores --export=ALL make unit-run
# Save elapsed time alongside the test results.
elapsed="$(( $(date +%s) - started ))"
printf '%s\n' "$elapsed" > "$CI_WORKSPACE/unit-test-seconds.txt"
echo "Unit-test runtime: $elapsed seconds"
