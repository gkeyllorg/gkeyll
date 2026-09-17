#!/usr/bin/env bash
# Dispatch to the Jenkins client for a locally selected CI platform.

set -euo pipefail

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<'EOF'
Usage:
  gkeyll-ci.sh <platform> <command> [flags]
  gkeyll-ci.sh <command> -h

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

Use `gkeyll-ci.sh <command> -h` for flags shared by CI platforms, or
`gkeyll-ci.sh <platform> <command> -h` for platform-specific command help.
This wrapper forwards platform-qualified commands and environment variables
unchanged to the selected client.
EOF
}

command_usage() {
    case "$1" in
        run) cat <<'EOF'
Usage: gkeyll-ci.sh run (--pr NUMBER | --candidate-ref REF --baseline-ref REF) [--follow]

Flags:
  --pr NUMBER           Build GitHub pull request NUMBER.
  --candidate-ref REF   Candidate branch or commit; requires --baseline-ref.
  --baseline-ref REF    Baseline branch or commit; requires --candidate-ref.
  --follow              Stream the build console after Jenkins queues it.

Specify a platform to run the build: gkeyll-ci.sh <platform> run ...
EOF
        ;;
        active) cat <<'EOF'
Usage: gkeyll-ci.sh active

This command takes no flags. Specify a platform: gkeyll-ci.sh <platform> active
EOF
        ;;
        recent) cat <<'EOF'
Usage: gkeyll-ci.sh recent [--limit NUMBER]

Flags:
  --limit NUMBER        Number of retained builds to list (default: 10).

Specify a platform to list builds: gkeyll-ci.sh <platform> recent ...
EOF
        ;;
        status) cat <<'EOF'
Usage: gkeyll-ci.sh status (--queue ID | --build NUMBER)

Flags:
  --queue ID            Inspect a Jenkins queue item.
  --build NUMBER        Inspect a known Jenkins build.

Specify a platform to inspect the build: gkeyll-ci.sh <platform> status ...
EOF
        ;;
        follow) cat <<'EOF'
Usage: gkeyll-ci.sh follow (--queue ID | --build NUMBER)

Flags:
  --queue ID            Wait for this queue item, then stream its build.
  --build NUMBER        Stream this known Jenkins build.

Specify a platform to follow the build: gkeyll-ci.sh <platform> follow ...
EOF
        ;;
        abort) cat <<'EOF'
Usage: gkeyll-ci.sh abort (--queue ID | --build NUMBER)

Flags:
  --queue ID            Cancel this Jenkins queue item.
  --build NUMBER        Abort this running Jenkins build.

Specify a platform to cancel the build: gkeyll-ci.sh <platform> abort ...
EOF
        ;;
        start) cat <<'EOF'
Usage: gkeyll-ci.sh start

This command takes no flags and is available on stellar_cpu and perlmutter_gpu.
Specify a platform: gkeyll-ci.sh <platform> start
EOF
        ;;
        scan) cat <<'EOF'
Usage: gkeyll-ci.sh scan

This command takes no flags and is available on team.
Specify a platform: gkeyll-ci.sh team scan
EOF
        ;;
    esac
}

die() {
    echo "ERROR: $*" >&2
    exit 2
}

main() {
    (($# > 0)) || { usage >&2; exit 2; }

    local platform="$1"
    case "$platform" in
        scan|start|run|active|recent|status|follow|abort)
            if (($# == 2)) && [[ "$2" == -h || "$2" == --help ]]; then
                command_usage "$platform"
                return
            fi
            usage >&2
            die "specify a platform before the $platform command"
            ;;
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
