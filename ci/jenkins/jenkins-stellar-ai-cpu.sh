#!/usr/bin/env bash
# Stellar AI defaults; share the authenticated CLI and controller lifecycle.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export GKEYLL_CI_ROOT="${GKEYLL_CI_ROOT:-/scratch/gpfs/CSD/${USER:?}/gkeyll_ci_stellar_ai_cpu}"
export JAVA_HOME="${JAVA_HOME:-/usr/lib/jvm/jre-21-openjdk}"
export JENKINS_JOB="${JENKINS_JOB:-gkeyll-ci-stellar-ai-cpu}"
export JENKINS_SESSION="${JENKINS_SESSION:-gkeyll_ci_stellar_ai_cpu}"
export GKEYLL_CI_CLIENT_NAME=jenkins-stellar-ai-cpu.sh
export GKEYLL_CI_CONTROLLER_LAUNCHER="$script_dir/jenkins-stellar-ai-cpu.sh"
exec "$script_dir/jenkins-stellar_cpu.sh" "$@"
