#!/bin/bash -l
# Unit-test payload for ci/jenkins/Jenkinsfile.stellar_cpu.
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
. machines/module_load.stellar-intel.sh

make unit-run
