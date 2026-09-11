#!/bin/bash -l
# Unit-test payload for ci/jenkins/Jenkinsfile.stellar-intel.
#
# This script deliberately has no #SBATCH directives. The trusted Jenkins
# pipeline owns the resource request, while this script owns only the runtime
# environment and the test command.
set -euo pipefail

: "${CI_WORKSPACE:?CI_WORKSPACE must name the shared Jenkins workspace}"

cd "$CI_WORKSPACE"

# A batch shell does not inherit a user's interactive module selection in a
# reliable or reproducible way. Keep these in sync with
# machines/configure.stellar-intel.sh.
module purge
module load intel/2022.2.0
module load openmpi/intel-2022.0/4.1.8
module load openblas/0.3.x

make unit-run
