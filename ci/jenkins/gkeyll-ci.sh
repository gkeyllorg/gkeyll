#!/usr/bin/env bash
# Dispatch to the Jenkins client for a locally selected CI platform.

set -euo pipefail

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<'EOF'
Usage:
  gkeyll-ci.sh <platform> <command> [flags]

Platforms:
  personal       Local personal-computer Jenkins CI.
  stellar_cpu    Stellar CPU Slurm Jenkins CI.
  perlmutter_gpu Perlmutter GPU Slurm Jenkins CI.
  team           Local team-workstation Jenkins CI.

Commands:
  scan           Discover branches and report scheduled builds (team only).
  start          Start the Jenkins controller (stellar_cpu, perlmutter_gpu).
  run            Queue a pull-request or candidate/baseline comparison build.
  active         List queued and running work.
  recent         List retained builds.
  status         Show the state of a queued or known build.
  follow         Wait for and stream a queued or known build.
  abort          Cancel a queued or running build.

Examples:
  gkeyll-ci.sh personal run --pr 1128 --follow
  gkeyll-ci.sh stellar_cpu recent
  gkeyll-ci.sh perlmutter_gpu status --queue 42
  gkeyll-ci.sh team scan

This wrapper runs the selected platform client on the current machine. It
forwards all remaining arguments and environment variables unchanged; use the
selected client's --help for its complete command reference and setup needs.
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 2
}

main() {
    (($# > 0)) || { usage >&2; exit 2; }

    local platform="$1"
    case "$platform" in
        -h|--help|help)
            (($# == 1)) || die 'usage: gkeyll-ci.sh <platform> <command> [flags]'
            usage
            return
            ;;
        personal) client='jenkins-personal.sh' ;;
        stellar_cpu) client='jenkins-stellar_cpu.sh' ;;
        perlmutter_gpu) client='jenkins-perlmutter_gpu.sh' ;;
        team) client='jenkins-team-workstation.sh' ;;
        *) usage >&2; die "unknown CI platform: $platform" ;;
    esac
    shift

    (($# > 0)) || die "provide a command for the $platform CI client (try: gkeyll-ci.sh $platform --help)"
    exec "$SCRIPT_DIR/$client" "$@"
}

main "$@"
