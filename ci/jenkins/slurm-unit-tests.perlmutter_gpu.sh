#!/bin/bash -l
# Unit-test payload for ci/jenkins/jenkinsfile.perlmutter_gpu.
#
# This script deliberately has no #SBATCH directives. The trusted Jenkins
# pipeline owns the resource request, while this script owns only the runtime
# environment and the test command.
set -euo pipefail

: "${CI_WORKSPACE:?CI_WORKSPACE must name the shared Jenkins workspace}"

cd "$CI_WORKSPACE"

# A batch shell does not inherit a user's interactive module selection in a
# reliable or reproducible way. Source the same environment used to configure
# and build this checkout.
. machines/module_load.perlmutter-gpu.sh

export SLURM_CPU_BIND=cores
started="$(date +%s)"
srun --ntasks=1 --cpus-per-task=32 --gpus-per-task=1 --cpu-bind=cores make unit-run
elapsed="$(( $(date +%s) - started ))"
printf '%s\n' "$elapsed" > unit-test-seconds.txt
echo "Unit-test runtime: $elapsed seconds"
