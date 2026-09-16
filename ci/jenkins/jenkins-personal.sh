#!/usr/bin/env bash
# Submit and inspect the local personal Jenkins Pipeline. Jenkins itself must
# already be running (for example as a Homebrew or systemd service).
set -euo pipefail

JENKINS_URL="${JENKINS_URL:-http://127.0.0.1:8080}"
JENKINS_JOB="${JENKINS_JOB:-gkeyll-ci-personal}"
JENKINS_CLI_AUTH_FILE="${JENKINS_CLI_AUTH_FILE:-}"
JENKINS_CLI_JAR="${JENKINS_CLI_JAR:-${TMPDIR:-/tmp}/gkeyll-jenkins-cli.jar}"
CURL_CONFIG=''
QUEUE_ID=''
RESOLVED_BUILD_NUMBER=''

die() { echo "ERROR: $*" >&2; exit 1; }
usage() {
    cat <<'EOF'
Usage: jenkins-personal.sh <command> [flags]

Commands:
  run --pr NUMBER [--follow]                 Queue a GitHub pull-request build.
  run --candidate-ref REF --baseline-ref REF [--follow]
                                             Queue a branch or commit comparison.
  active                                     List this job's queued and running work.
  recent [--limit NUMBER]                    List retained builds (default: 10).
  status --queue ID                          Show a queued build's current state.
  status --build NUMBER                      Show a known build's current state.
  follow --queue ID                          Wait for and stream a queued build.
  follow --build NUMBER                      Stream a known Jenkins build.
  abort --queue ID                           Cancel a queued Jenkins build.
  abort --build NUMBER                       Abort a running Jenkins build.

The run command returns after Jenkins accepts the request. --follow streams
the build console and returns its final Jenkins result. Press Ctrl-C to stop
following without aborting the Jenkins build.

Jenkins must already be running. Set JENKINS_CLI_AUTH_FILE to a protected file
containing one line: jenkins-user:api-token. JENKINS_URL and JENKINS_JOB
override the loopback URL and gkeyll-ci-personal defaults.
EOF
}
job_path() { local p='/job' n; IFS=/ read -ra n <<< "$JENKINS_JOB"; for x in "${n[@]}"; do p+="/$x/job"; done; printf '%s' "${p%/job}"; }
prepare_auth() {
    [[ -n "$JENKINS_CLI_AUTH_FILE" ]] || die 'Set JENKINS_CLI_AUTH_FILE to a mode-600 user:api-token file'
    [[ -O "$JENKINS_CLI_AUTH_FILE" ]] || die "credential file is not owned by $USER"
    [[ "$(stat -f '%Lp' "$JENKINS_CLI_AUTH_FILE" 2>/dev/null || stat -c '%a' "$JENKINS_CLI_AUTH_FILE")" == 600 ]] || die 'credential file must have mode 600'
    local credential; credential="$(<"$JENKINS_CLI_AUTH_FILE")"
    [[ "$credential" =~ ^[^[:space:]:]+:[^[:space:]:]+$ ]] || die 'credential file must contain user:api-token'
    CURL_CONFIG="$(mktemp "${TMPDIR:-/tmp}/gkeyll-jenkins.XXXXXX")"; chmod 600 "$CURL_CONFIG"
    printf 'user = "%s"\n' "$credential" > "$CURL_CONFIG"
    trap 'rm -f "$CURL_CONFIG"' EXIT
}
curl_auth() { curl --fail --silent --show-error --globoff --config "$CURL_CONFIG" "$@"; }
curl_auth_quiet() { curl --fail --silent --globoff --config "$CURL_CONFIG" "$@"; }
positive() { [[ "$2" =~ ^[1-9][0-9]*$ ]] || die "$1 must be a positive integer"; }
java_command() {
    if [[ -n "${JAVA_HOME:-}" && -x "$JAVA_HOME/bin/java" ]]; then
        printf '%s' "$JAVA_HOME/bin/java"
    else
        command -v java || die 'Set JAVA_HOME or put Java 21 or newer on PATH to stream the Jenkins console'
    fi
}
download_cli() {
    if [[ ! -s "$JENKINS_CLI_JAR" ]]; then
        local directory temporary
        directory="$(dirname "$JENKINS_CLI_JAR")"
        mkdir -p "$directory"
        temporary="$(mktemp "$directory/gkeyll-jenkins-cli.XXXXXX")"
        curl --fail --silent --show-error --output "$temporary" "$JENKINS_URL/jnlpJars/jenkins-cli.jar"
        mv "$temporary" "$JENKINS_CLI_JAR"
    fi
}
state() {
    local payload
    payload="$(curl_auth "${JENKINS_URL}$(job_path)/$1/api/json")" || return 75
    python3 -c 'import json,sys; x=json.load(sys.stdin); print("BUILDING" if x["building"] else (x.get("result") or "UNKNOWN"))' <<< "$payload"
}
queue_state() {
    local queue_id="$1" payload
    payload="$(curl_auth_quiet "${JENKINS_URL}/queue/item/$queue_id/api/json")" || return 1
    python3 -c '
import json,sys
item=json.load(sys.stdin); executable=item.get("executable") or {}
print("{}\t{}\t{}".format(executable.get("number", ""), "true" if item.get("cancelled", False) else "false", " ".join((item.get("why") or "").split())))
' <<< "$payload"
}
build_for_queue() {
    local queue_id="$1" payload
    payload="$(curl_auth "${JENKINS_URL}$(job_path)/api/json?tree=builds[number,queueId]")" || return 1
    python3 -c '
import json,sys
queue_id=int(sys.argv[1])
for build in json.load(sys.stdin).get("builds",[]):
    if build.get("queueId") == queue_id:
        print(build["number"]); break
' "$queue_id" <<< "$payload"
}
wait_for_build_number() {
    local queue_id="$1" record number cancelled why previous='' attempt
    while :; do
        if ! record="$(queue_state "$queue_id")"; then
            for attempt in {1..6}; do
                if number="$(build_for_queue "$queue_id")" && [[ -n "$number" ]]; then
                    RESOLVED_BUILD_NUMBER="$number"
                    return 0
                fi
                sleep 2
            done
            die "Queue item $queue_id is unavailable and no matching Jenkins build was found"
        fi
        IFS=$'\t' read -r number cancelled why <<< "$record"
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
status_build() {
    positive 'build number' "$1"
    local s
    s="$(state "$1")" || die "Build #$1 is not available yet. Retry shortly."
    echo "Build #$1: $s"
    [[ "$s" == SUCCESS ]]
}
status_command() {
    [[ $# == 2 ]] || die 'usage: status --queue ID | status --build NUMBER'
    case "$1" in
        --queue)
            positive 'queue ID' "$2"
            local record number cancelled why
            if ! record="$(queue_state "$2")"; then
                number="$(build_for_queue "$2")" || die "Queue item $2 is unavailable"
                [[ -n "$number" ]] || die "Queue item $2 is unavailable"
                echo "Queue item $2 is build #$number"
                status_build "$number"
                return
            fi
            IFS=$'\t' read -r number cancelled why <<< "$record"
            if [[ -n "$number" ]]; then
                echo "Queue item $2 is build #$number"
                status_build "$number"
            elif [[ "$cancelled" == true ]]; then
                echo "Queue item $2: CANCELLED"
                return 1
            else
                echo "Queue item $2: ${why:-waiting}"
            fi
            ;;
        --build) status_build "$2" ;;
        *) die 'usage: status --queue ID | status --build NUMBER' ;;
    esac
}
abort() {
    [[ $# == 2 ]] || die 'usage: abort --queue ID | abort --build NUMBER'
    case "$1" in
        --queue)
            positive 'queue ID' "$2"
            curl_auth --output /dev/null --request POST --data-urlencode "id=$2" \
                "${JENKINS_URL}/queue/cancelItem"
            echo "Requested cancellation of queue item $2"
            ;;
        --build)
            positive 'build number' "$2"
            curl_auth --output /dev/null --request POST \
                "${JENKINS_URL}$(job_path)/$2/stop"
            echo "Requested cancellation of $JENKINS_JOB #$2"
            ;;
        *) die 'usage: abort --queue ID | abort --build NUMBER' ;;
    esac
}
follow() {
    positive 'build number' "$1"
    local java
    java="$(java_command)"
    download_cli
    echo "Following $JENKINS_JOB #$1"
    echo "Build URL: ${JENKINS_URL}$(job_path)/$1/"
    trap "echo 'Stopped following build #$1'; exit 130" INT TERM HUP
    if ! "$java" -jar "$JENKINS_CLI_JAR" -s "$JENKINS_URL" -http \
        -auth "@$JENKINS_CLI_AUTH_FILE" console "$JENKINS_JOB" "$1" -f; then
        echo 'Console follower ended before Jenkins reported a terminal result; checking build status.' >&2
    fi
    trap - INT TERM HUP
    while :; do
        local s
        if ! s="$(state "$1")"; then
            echo "Build #$1 is not available yet; retrying." >&2
            sleep 2
            continue
        fi
        [[ "$s" == BUILDING ]] || { echo "Build #$1: $s"; [[ "$s" == SUCCESS ]]; return; }
        sleep 5
    done
}
follow_command() {
    [[ $# == 2 ]] || die 'usage: follow --queue ID | follow --build NUMBER'
    case "$1" in
        --queue) positive 'queue ID' "$2"; wait_for_build_number "$2"; follow "$RESOLVED_BUILD_NUMBER" ;;
        --build) follow "$2" ;;
        *) die 'usage: follow --queue ID | follow --build NUMBER' ;;
    esac
}
submit() {
    local pr="$1" candidate="$2" baseline="$3" headers queue
    headers="$(mktemp "${TMPDIR:-/tmp}/gkeyll-jenkins-headers.XXXXXX")"
    curl_auth --dump-header "$headers" --output /dev/null --request POST \
        --data-urlencode "CANDIDATE_PR=$pr" --data-urlencode "CANDIDATE_REF=$candidate" --data-urlencode "BASELINE_REF=$baseline" \
        "${JENKINS_URL}$(job_path)/buildWithParameters" || { rm -f "$headers"; die 'Jenkins rejected the build'; }
    queue="$(awk 'BEGIN{IGNORECASE=1} /^Location:/{sub(/^[^:]*: /,""); sub(/\r$/,""); print; exit}' "$headers")"; rm -f "$headers"
    [[ "$queue" =~ /queue/item/([0-9]+)/ ]] || die 'Jenkins accepted the build but returned no queue ID'
    QUEUE_ID="${BASH_REMATCH[1]}"
    echo "Queued $JENKINS_JOB as queue item $QUEUE_ID"
}
run() {
    local pr='' candidate='' baseline='' want_follow=false
    while (($#)); do case "$1" in --pr) (($#>=2))||die '--pr requires a number'; pr="$2"; shift 2;; --candidate-ref) (($#>=2))||die '--candidate-ref requires a ref'; candidate="$2"; shift 2;; --baseline-ref) (($#>=2))||die '--baseline-ref requires a ref'; baseline="$2"; shift 2;; --follow) want_follow=true; shift;; *) die "unknown run option: $1";; esac; done
    if [[ -n "$pr" ]]; then positive '--pr' "$pr"; [[ -z "$candidate$baseline" ]] || die '--pr cannot be combined with refs'; else [[ -n "$candidate" && -n "$baseline" ]] || die 'provide --pr, or both --candidate-ref and --baseline-ref'; fi
    submit "$pr" "$candidate" "$baseline"
    if [[ "$want_follow" == true ]]; then
        wait_for_build_number "$QUEUE_ID"
        follow "$RESOLVED_BUILD_NUMBER"
    fi
}
job_builds() { curl_auth "${JENKINS_URL}$(job_path)/api/json?tree=builds[number,queueId,building,result,timestamp,url,actions[parameters[name,value]]]"; }
queue_items() { curl_auth "${JENKINS_URL}/queue/api/json?tree=items[id,task[name],why,cancelled,executable[number],actions[parameters[name,value]]]"; }
format_builds() {
    local limit="$1" active_only="$2" payload="$3"
    python3 -c '
import datetime,json,sys
limit=int(sys.argv[1]); active_only=sys.argv[2] == "true"; count=0
def selectors(actions):
    values={}
    for action in actions or []:
        for parameter in action.get("parameters") or []:
            if parameter.get("name") in ("CANDIDATE_PR","CANDIDATE_REF","BASELINE_REF"):
                values[parameter["name"]]=parameter.get("value") or "-"
    return "pr={0} candidate={1} baseline={2}".format(values.get("CANDIDATE_PR","-"),values.get("CANDIDATE_REF","-"),values.get("BASELINE_REF","-"))
for build in json.load(sys.stdin).get("builds",[]):
    if active_only and not build.get("building",False): continue
    if count >= limit: break
    timestamp=build.get("timestamp")
    when="-" if timestamp is None else datetime.datetime.fromtimestamp(timestamp/1000.0,datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    state="RUNNING" if build.get("building",False) else (build.get("result") or "UNKNOWN")
    print("BUILD #{number} {state} queue={queue} {when} {selectors} url={url}".format(number=build.get("number","?"),state=state,queue=build.get("queueId") if build.get("queueId") is not None else "-",when=when,selectors=selectors(build.get("actions")),url=build.get("url","-")))
    count+=1
' "$limit" "$active_only" <<< "$payload"
}
format_queue_items() {
    local payload="$1"
    python3 -c '
import json,sys
job=sys.argv[1]
def selectors(actions):
    values={}
    for action in actions or []:
        for parameter in action.get("parameters") or []:
            if parameter.get("name") in ("CANDIDATE_PR","CANDIDATE_REF","BASELINE_REF"):
                values[parameter["name"]]=parameter.get("value") or "-"
    return "pr={0} candidate={1} baseline={2}".format(values.get("CANDIDATE_PR","-"),values.get("CANDIDATE_REF","-"),values.get("BASELINE_REF","-"))
for item in json.load(sys.stdin).get("items",[]):
    if (item.get("task") or {}).get("name") != job or item.get("executable"): continue
    state="CANCELLED" if item.get("cancelled",False) else "QUEUED"
    print("QUEUE {id} {state} reason={why} {selectors}".format(id=item.get("id","?"),state=state,why=" ".join((item.get("why") or "waiting").split()),selectors=selectors(item.get("actions"))))
' "$JENKINS_JOB" <<< "$payload"
}
active() {
    local builds queues output
    builds="$(job_builds)"; queues="$(queue_items)"
    output="$(format_builds 100 true "$builds")"
    output+=$'\n'"$(format_queue_items "$queues")"
    [[ -n "${output//$'\n'/}" ]] && printf '%s\n' "$output" || echo "No queued or running builds for $JENKINS_JOB"
}
recent() {
    local limit="${1:-10}" output
    positive '--limit' "$limit"
    output="$(format_builds "$limit" false "$(job_builds)")"
    [[ -n "$output" ]] && printf '%s\n' "$output" || echo "No retained builds for $JENKINS_JOB"
}
main() {
    (($#)) || { usage; exit 2; }
    case "$1" in -h|--help|help) usage; return;; esac
    prepare_auth
    case "$1" in run) shift; run "$@";; follow) shift; follow_command "$@";; status) shift; status_command "$@";; abort) shift; abort "$@";; active) shift; [[ $# == 0 ]] || die 'usage: active'; active;; recent) shift; [[ $# == 0 || ( $# == 2 && $1 == --limit ) ]] || die 'usage: recent [--limit NUMBER]'; recent "${2:-10}";; -h|--help|help) usage;; *) usage >&2; die "unknown command: $1";; esac
}
main "$@"
