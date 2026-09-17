#!/usr/bin/env bash
#
# Start and use the private, loopback-only Jenkins controller on Stellar CPU.
# This is intentionally a controller launcher, not a replacement for the
# Jenkins UI: both clients submit the same parameterized Pipeline job.

set -euo pipefail

readonly SCRIPT_PATH="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/$(basename -- "${BASH_SOURCE[0]}")"
readonly DEFAULT_JAVA_HOME=/usr/lib/jvm/java-21-openjdk-21.0.12.1.1-1.1.el8.x86_64

CI_ROOT="${GKEYLL_CI_ROOT:-/scratch/gpfs/${USER:?USER must be set}/gkeyll_ci}"
JENKINS_HOME="${JENKINS_HOME:-$CI_ROOT/jenkins_home}"
JENKINS_WEBROOT="${JENKINS_WEBROOT:-$CI_ROOT/jenkins_webroot}"
JENKINS_TMPDIR="${TMPDIR:-$CI_ROOT/tmp}"
JENKINS_PORT="${JENKINS_PORT:-8080}"
JENKINS_URL="${JENKINS_URL:-http://127.0.0.1:$JENKINS_PORT}"
JENKINS_JOB="${JENKINS_JOB:-gkeyll-ci-stellar_cpu}"
JENKINS_SESSION="${JENKINS_SESSION:-gkeyll_ci}"
JENKINS_CLI_AUTH_FILE="${JENKINS_CLI_AUTH_FILE:-$JENKINS_HOME/jenkins-cli.auth}"
JAVA_HOME="${JAVA_HOME:-$DEFAULT_JAVA_HOME}"
CLI_JAR="$CI_ROOT/jenkins-cli.jar"

CURL_CONFIG=''
SUBMITTED_QUEUE_ID=''
RESOLVED_BUILD_NUMBER=''

die() {
    echo "ERROR: $*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage:
  jenkins-stellar_cpu.sh <command> [flags]

Commands:
  start    Start Jenkins in detached tmux.
  run      Queue a pull-request or comparison build.
  active   List this job's queued and running work.
  recent   List retained builds.
  status   Show a queued or known build's state.
  follow   Wait for and stream a queued or known build.
  abort    Cancel a queued or running build.

The run command returns after Jenkins accepts the request. --follow streams
the build console and returns its final Jenkins result. Press Ctrl-C to stop
following without aborting the Jenkins build.

GKEYLL_CI_ROOT defaults to /scratch/gpfs/$USER/gkeyll_ci. The Jenkins API
credential file defaults to $JENKINS_HOME/jenkins-cli.auth and must contain
one line in the form: jenkins-user:api-token
EOF
}

command_usage() {
    local script='jenkins-stellar_cpu.sh'
    case "$1" in
        start) cat <<EOF
Usage: $script start

This command takes no options.
EOF
        ;;
        run) cat <<EOF
Usage: $script run (--pr NUMBER | --candidate-ref REF --baseline-ref REF) [--follow]

Flags:
  --pr NUMBER           Build GitHub pull request NUMBER.
  --candidate-ref REF   Candidate branch or commit; requires --baseline-ref.
  --baseline-ref REF    Baseline branch or commit; requires --candidate-ref.
  --follow              Stream the build console after Jenkins queues it.
EOF
        ;;
        active) cat <<EOF
Usage: $script active

This command takes no options.
EOF
        ;;
        recent) cat <<EOF
Usage: $script recent [--limit NUMBER]

Flags:
  --limit NUMBER        Number of retained builds to list (default: 10).
EOF
        ;;
        status) cat <<EOF
Usage: $script status (--queue ID | --build NUMBER)

Flags:
  --queue ID            Inspect a Jenkins queue item.
  --build NUMBER        Inspect a known Jenkins build.
EOF
        ;;
        follow) cat <<EOF
Usage: $script follow (--queue ID | --build NUMBER)

Flags:
  --queue ID            Wait for this queue item, then stream its build.
  --build NUMBER        Stream this known Jenkins build.
EOF
        ;;
        abort) cat <<EOF
Usage: $script abort (--queue ID | --build NUMBER)

Flags:
  --queue ID            Cancel this Jenkins queue item.
  --build NUMBER        Abort this running Jenkins build.
EOF
        ;;
    esac
}

prepare_paths() {
    mkdir -p "$JENKINS_HOME" "$JENKINS_WEBROOT" "$JENKINS_TMPDIR" \
        "$CI_ROOT/logs" "$CI_ROOT/workspaces"
}

controller_running() {
    tmux has-session -t "$JENKINS_SESSION" 2>/dev/null
}

wait_for_controller() {
    local attempt
    for attempt in {1..30}; do
        # Connection refusals and 503s are expected while Jenkins initializes.
        # Keep individual retry failures quiet; the final error names the
        # controller session to inspect if all attempts fail.
        if curl --fail --silent --max-time 5 \
            --output /dev/null "$JENKINS_URL/login"; then
            return 0
        fi
        sleep 2
    done
    die "Jenkins did not become ready at $JENKINS_URL. Inspect: tmux attach -t $JENKINS_SESSION"
}

start_controller() {
    prepare_paths
    [[ -x "$JAVA_HOME/bin/java" ]] || die "Java 21 was not found at $JAVA_HOME/bin/java"
    [[ -f "$CI_ROOT/jenkins.war" ]] || die "Jenkins WAR is missing: $CI_ROOT/jenkins.war"

    if controller_running; then
        echo "Jenkins tmux session already exists: $JENKINS_SESSION"
    else
        tmux new-session -d -s "$JENKINS_SESSION" "$SCRIPT_PATH" __controller
        echo "Started Jenkins in tmux session: $JENKINS_SESSION"
    fi
    wait_for_controller
    echo "Jenkins is ready at $JENKINS_URL"
}

run_controller() {
    prepare_paths
    [[ -x "$JAVA_HOME/bin/java" ]] || die "Java 21 was not found at $JAVA_HOME/bin/java"
    [[ -f "$CI_ROOT/jenkins.war" ]] || die "Jenkins WAR is missing: $CI_ROOT/jenkins.war"

    export JENKINS_HOME JENKINS_WEBROOT
    export TMPDIR="$JENKINS_TMPDIR"
    export PATH="$JAVA_HOME/bin:$PATH"
    "$JAVA_HOME/bin/java" -Djava.io.tmpdir="$TMPDIR" -jar "$CI_ROOT/jenkins.war" \
        --webroot="$JENKINS_WEBROOT" \
        --httpListenAddress=127.0.0.1 \
        --httpPort="$JENKINS_PORT" \
        2>&1 | tee -a "$CI_ROOT/logs/jenkins.log"
}

cleanup_curl_config() {
    if [[ -n "$CURL_CONFIG" ]]; then
        rm -f "$CURL_CONFIG"
    fi
}

prepare_auth() {
    [[ -f "$JENKINS_CLI_AUTH_FILE" ]] || die "Jenkins API credential file is missing: $JENKINS_CLI_AUTH_FILE"
    [[ -O "$JENKINS_CLI_AUTH_FILE" ]] || die "Jenkins API credential file is not owned by $USER: $JENKINS_CLI_AUTH_FILE"

    local mode credential
    mode="$(stat -c '%a' "$JENKINS_CLI_AUTH_FILE")"
    [[ "$mode" == 600 ]] || die "Jenkins API credential file must have mode 600: $JENKINS_CLI_AUTH_FILE"
    credential="$(<"$JENKINS_CLI_AUTH_FILE")"
    [[ "$credential" =~ ^[^[:space:]:]+:[^[:space:]:]+$ ]] || die "Jenkins API credential file must contain exactly user:api-token"
    [[ "$credential" != *'"'* && "$credential" != *'\\'* ]] || die "Jenkins API credential file contains an unsupported character"

    CURL_CONFIG="$(mktemp "$JENKINS_TMPDIR/jenkins-curl.XXXXXX")"
    chmod 600 "$CURL_CONFIG"
    printf 'user = "%s"\n' "$credential" > "$CURL_CONFIG"
    trap cleanup_curl_config EXIT
}

curl_auth() {
    curl --fail --silent --show-error --globoff --config "$CURL_CONFIG" "$@"
}

curl_auth_quiet() {
    curl --fail --silent --globoff --config "$CURL_CONFIG" "$@"
}

download_cli() {
    if [[ ! -s "$CLI_JAR" ]]; then
        local temporary_jar
        temporary_jar="$(mktemp "$JENKINS_TMPDIR/jenkins-cli.XXXXXX")"
        curl --fail --silent --show-error --output "$temporary_jar" \
            "$JENKINS_URL/jnlpJars/jenkins-cli.jar"
        mv "$temporary_jar" "$CLI_JAR"
    fi
}

require_positive_integer() {
    local name="$1" value="$2"
    [[ "$value" =~ ^[1-9][0-9]*$ ]] || die "$name must be a positive integer, got '$value'"
}

queue_state() {
    local queue_id="$1" payload
    # Jenkins discards a queue item after assigning its build. That expected
    # transition is handled by build_for_queue below, so do not print curl's
    # 404 before attempting the recovery path.
    payload="$(curl_auth_quiet "$JENKINS_URL/queue/item/$queue_id/api/json")" || return 1
    python3 -c '
import json
import sys
item = json.load(sys.stdin)
executable = item.get("executable") or {}
number = executable.get("number", "")
cancelled = "true" if item.get("cancelled", False) else "false"
why = " ".join((item.get("why") or "").split())
print(f"{number}\t{cancelled}\t{why}")
' <<< "$payload"
}

build_for_queue() {
    local queue_id="$1" payload
    payload="$(curl_auth "$JENKINS_URL/job/$JENKINS_JOB/api/json?tree=builds[number,queueId]")" || return 1
    python3 -c '
import json
import sys

queue_id = int(sys.argv[1])
for build in json.load(sys.stdin).get("builds", []):
    if build.get("queueId") == queue_id:
        print(build["number"])
        break
' "$queue_id" <<< "$payload"
}

build_state() {
    local build_number="$1" payload
    payload="$(curl_auth "$JENKINS_URL/job/$JENKINS_JOB/$build_number/api/json")" || return 1
    python3 -c '
import json
import sys
build = json.load(sys.stdin)
print("true" if build.get("building", False) else "false")
print(build.get("result") or "")
' <<< "$payload"
}

job_builds() {
    curl_auth "$JENKINS_URL/job/$JENKINS_JOB/api/json?tree=builds[number,queueId,building,result,timestamp,url,actions[parameters[name,value]]]"
}

queue_items() {
    curl_auth "$JENKINS_URL/queue/api/json?tree=items[id,task[name],why,cancelled,executable[number],actions[parameters[name,value]]]"
}

format_builds() {
    local limit="$1" active_only="$2" payload="$3"
    python3 -c '
import datetime
import json
import sys

limit = int(sys.argv[1])
active_only = sys.argv[2] == "true"

def selectors(actions):
    values = {}
    for action in actions or []:
        for parameter in action.get("parameters") or []:
            name = parameter.get("name")
            if name in ("CANDIDATE_PR", "CANDIDATE_REF", "BASELINE_REF"):
                values[name] = parameter.get("value") or "-"
    return "pr={0} candidate={1} baseline={2}".format(
        values.get("CANDIDATE_PR", "-"),
        values.get("CANDIDATE_REF", "-"),
        values.get("BASELINE_REF", "-"))

count = 0
for build in json.load(sys.stdin).get("builds", []):
    if active_only and not build.get("building", False):
        continue
    if count >= limit:
        break
    state = "RUNNING" if build.get("building", False) else (build.get("result") or "UNKNOWN")
    timestamp = build.get("timestamp")
    when = "-"
    if timestamp is not None:
        when = datetime.datetime.fromtimestamp(
            timestamp / 1000.0, datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    queue_id = build.get("queueId")
    print("BUILD #{number} {state} queue={queue} {when} {selectors} url={url}".format(
        number=build.get("number", "?"), state=state, queue=queue_id if queue_id is not None else "-",
        when=when, selectors=selectors(build.get("actions")), url=build.get("url", "-")))
    count += 1
' "$limit" "$active_only" <<< "$payload"
}

format_queue_items() {
    local payload="$1"
    python3 -c '
import json
import sys

job = sys.argv[1]

def selectors(actions):
    values = {}
    for action in actions or []:
        for parameter in action.get("parameters") or []:
            name = parameter.get("name")
            if name in ("CANDIDATE_PR", "CANDIDATE_REF", "BASELINE_REF"):
                values[name] = parameter.get("value") or "-"
    return "pr={0} candidate={1} baseline={2}".format(
        values.get("CANDIDATE_PR", "-"),
        values.get("CANDIDATE_REF", "-"),
        values.get("BASELINE_REF", "-"))

for item in json.load(sys.stdin).get("items", []):
    if (item.get("task") or {}).get("name") != job:
        continue
    # A started item is displayed as its running build, not twice.
    if item.get("executable"):
        continue
    state = "CANCELLED" if item.get("cancelled", False) else "QUEUED"
    why = " ".join((item.get("why") or "waiting").split())
    print("QUEUE {id} {state} reason={why} {selectors}".format(
        id=item.get("id", "?"), state=state, why=why,
        selectors=selectors(item.get("actions"))))
' "$JENKINS_JOB" <<< "$payload"
}

wait_for_build_number() {
    local queue_id="$1" state number cancelled why previous='' attempt
    while :; do
        if ! state="$(queue_state "$queue_id")"; then
            # Jenkins normally removes a queue item as soon as it starts its
            # build. Build records retain queueId, so recover the assignment
            # from the job rather than requiring the caller to discover it.
            for attempt in {1..6}; do
                if number="$(build_for_queue "$queue_id")" && [[ -n "$number" ]]; then
                    RESOLVED_BUILD_NUMBER="$number"
                    return 0
                fi
                sleep 2
            done
            die "Queue item $queue_id is unavailable and no matching Jenkins build was found"
        fi
        IFS=$'\t' read -r number cancelled why <<< "$state"
        [[ "$cancelled" == false ]] || die "Queue item $queue_id was cancelled"
        if [[ -n "$number" ]]; then
            RESOLVED_BUILD_NUMBER="$number"
            return 0
        fi
        if [[ "$why" != "$previous" ]]; then
            echo "Queue item $queue_id: ${why:-waiting}" >&2
            previous="$why"
        fi
        sleep 5
    done
}

follow_build() {
    local build_number="$1" state building result
    require_positive_integer 'build number' "$build_number"
    download_cli
    echo "Following $JENKINS_JOB #$build_number"
    echo "Build URL: $JENKINS_URL/job/$JENKINS_JOB/$build_number/"

    # The CLI is only a console client. Let an interactive interrupt detach
    # this client cleanly; Jenkins in tmux and any Slurm jobs keep running.
    trap "echo 'Stopped following build #${build_number}'; exit 130" INT TERM HUP

    # -f follows console output without propagating a client interruption to
    # the Jenkins build. The controller and Slurm work remain independent of
    # the SSH terminal.
    if ! "$JAVA_HOME/bin/java" -jar "$CLI_JAR" -s "$JENKINS_URL" -http \
        -auth "@$JENKINS_CLI_AUTH_FILE" console "$JENKINS_JOB" "$build_number" -f; then
        echo "Console follower ended before Jenkins reported a terminal result; checking build status." >&2
    fi

    trap - INT TERM HUP

    while :; do
        mapfile -t state < <(build_state "$build_number")
        building="${state[0]:-}"
        result="${state[1]:-}"
        [[ "$building" == false ]] && break
        sleep 5
    done

    echo "Build #$build_number result: ${result:-UNKNOWN}"
    [[ "$result" == SUCCESS ]]
}

submit_build() {
    local candidate_pr="$1" candidate_ref="$2" baseline_ref="$3"
    local headers queue_url queue_id
    headers="$(mktemp "$JENKINS_TMPDIR/jenkins-headers.XXXXXX")"

    if ! curl_auth --dump-header "$headers" --output /dev/null --request POST \
        --data-urlencode "CANDIDATE_PR=$candidate_pr" \
        --data-urlencode "CANDIDATE_REF=$candidate_ref" \
        --data-urlencode "BASELINE_REF=$baseline_ref" \
        "$JENKINS_URL/job/$JENKINS_JOB/buildWithParameters"; then
        rm -f "$headers"
        die 'Jenkins rejected the build submission'
    fi
    queue_url="$(awk 'BEGIN { IGNORECASE=1 } /^Location:/ { sub(/^[^:]*:[[:space:]]*/, ""); sub(/\r$/, ""); print; exit }' "$headers")"
    rm -f "$headers"
    [[ "$queue_url" =~ /queue/item/([0-9]+)/ ]] || die "Jenkins accepted the request but did not return a queue location"
    queue_id="${BASH_REMATCH[1]}"

    echo "Queued $JENKINS_JOB as queue item $queue_id" >&2
    echo "Queue URL: $queue_url" >&2
    echo "Monitor it with: $SCRIPT_PATH follow --queue $queue_id" >&2
    SUBMITTED_QUEUE_ID="$queue_id"
}

run_command() {
    local candidate_pr='' candidate_ref='' baseline_ref='' follow=false queue_id
    while (($#)); do
        case "$1" in
            --pr)
                (($# >= 2)) || die '--pr requires a number'
                candidate_pr="$2"
                shift 2
                ;;
            --candidate-ref)
                (($# >= 2)) || die '--candidate-ref requires a ref'
                candidate_ref="$2"
                shift 2
                ;;
            --baseline-ref)
                (($# >= 2)) || die '--baseline-ref requires a ref'
                baseline_ref="$2"
                shift 2
                ;;
            --follow)
                follow=true
                shift
                ;;
            *) die "unknown run option: $1" ;;
        esac
    done

    if [[ -n "$candidate_pr" ]]; then
        require_positive_integer '--pr' "$candidate_pr"
        [[ -z "$candidate_ref$baseline_ref" ]] || die '--pr cannot be combined with candidate/baseline refs'
    else
        [[ -n "$candidate_ref" && -n "$baseline_ref" ]] || die 'provide --pr, or both --candidate-ref and --baseline-ref'
    fi

    start_controller
    prepare_auth
    submit_build "$candidate_pr" "$candidate_ref" "$baseline_ref"
    queue_id="$SUBMITTED_QUEUE_ID"
    if [[ "$follow" == true ]]; then
        local build_number
        wait_for_build_number "$queue_id"
        build_number="$RESOLVED_BUILD_NUMBER"
        follow_build "$build_number"
    fi
}

follow_command() {
    [[ $# -eq 2 ]] || die 'usage: follow --queue ID | follow --build NUMBER'
    start_controller
    prepare_auth
    case "$1" in
        --queue)
            require_positive_integer 'queue ID' "$2"
            wait_for_build_number "$2"
            follow_build "$RESOLVED_BUILD_NUMBER"
            ;;
        --build)
            follow_build "$2"
            ;;
        *) die 'usage: follow --queue ID | follow --build NUMBER' ;;
    esac
}

status_command() {
    [[ $# -eq 2 ]] || die 'usage: status --queue ID | status --build NUMBER'
    start_controller
    prepare_auth
    case "$1" in
        --queue)
            require_positive_integer 'queue ID' "$2"
            local state number cancelled why
            if ! state="$(queue_state "$2")"; then
                # Jenkins removes a queue item once it assigns a build. Recover
                # that build from its retained queueId, as follow --queue does.
                number="$(build_for_queue "$2")" || die "Queue item $2 is unavailable"
                [[ -n "$number" ]] || die "Queue item $2 is unavailable"
                echo "Queue item $2 is build #$number"
                status_command --build "$number"
                return
            fi
            IFS=$'\t' read -r number cancelled why <<< "$state"
            if [[ -n "$number" ]]; then
                echo "Queue item $2 is build #$number"
                status_command --build "$number"
            elif [[ "$cancelled" == true ]]; then
                echo "Queue item $2: CANCELLED"
                return 1
            else
                echo "Queue item $2: ${why:-waiting}"
            fi
            ;;
        --build)
            require_positive_integer 'build number' "$2"
            local build
            mapfile -t build < <(build_state "$2")
            if [[ "${build[0]:-}" == true ]]; then
                echo "Build #$2: BUILDING"
            else
                echo "Build #$2: ${build[1]:-UNKNOWN}"
                [[ "${build[1]:-}" == SUCCESS ]]
            fi
            ;;
        *) die 'usage: status --queue ID | status --build NUMBER' ;;
    esac
}

abort_command() {
    [[ $# -eq 2 ]] || die 'usage: abort --queue ID | abort --build NUMBER'
    start_controller
    prepare_auth
    case "$1" in
        --queue)
            require_positive_integer 'queue ID' "$2"
            curl_auth --output /dev/null --request POST --data-urlencode "id=$2" \
                "$JENKINS_URL/queue/cancelItem"
            echo "Requested cancellation of queue item $2"
            ;;
        --build)
            require_positive_integer 'build number' "$2"
            curl_auth --output /dev/null --request POST \
                "$JENKINS_URL/job/$JENKINS_JOB/$2/stop"
            echo "Requested cancellation of $JENKINS_JOB #$2"
            ;;
        *) die 'usage: abort --queue ID | abort --build NUMBER' ;;
    esac
}

active_command() {
    local builds queues output
    [[ $# -eq 0 ]] || die 'usage: active'
    start_controller
    prepare_auth
    builds="$(job_builds)"
    queues="$(queue_items)"
    output="$(format_builds 100 true "$builds")"
    output+=$'\n'"$(format_queue_items "$queues")"
    if [[ -z "${output//$'\n'/}" ]]; then
        echo "No queued or running builds for $JENKINS_JOB"
    else
        printf '%s\n' "$output"
    fi
}

recent_command() {
    local limit=10 builds output
    if [[ $# -gt 0 ]]; then
        [[ $# -eq 2 && "$1" == --limit ]] || die 'usage: recent [--limit NUMBER]'
        limit="$2"
    fi
    require_positive_integer '--limit' "$limit"
    start_controller
    prepare_auth
    builds="$(job_builds)"
    output="$(format_builds "$limit" false "$builds")"
    if [[ -z "$output" ]]; then
        echo "No retained builds for $JENKINS_JOB"
    else
        printf '%s\n' "$output"
    fi
}

main() {
    (($# >= 1)) || { usage; exit 2; }
    if (($# == 2)) && [[ "$2" == -h || "$2" == --help ]]; then
        case "$1" in start|run|active|recent|status|follow|abort) command_usage "$1"; return;; esac
    fi
    case "$1" in
        __controller) shift; (($# == 0)) || die '__controller takes no arguments'; run_controller ;;
        start) shift; (($# == 0)) || die 'start takes no arguments'; start_controller ;;
        run) shift; run_command "$@" ;;
        follow) shift; follow_command "$@" ;;
        status) shift; status_command "$@" ;;
        abort) shift; abort_command "$@" ;;
        active) shift; active_command "$@" ;;
        recent) shift; recent_command "$@" ;;
        -h|--help|help) usage ;;
        *) usage >&2; die "unknown command: $1" ;;
    esac
}

main "$@"
