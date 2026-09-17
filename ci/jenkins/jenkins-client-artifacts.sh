#!/usr/bin/env bash
# Shared build-information and artifact helpers for Jenkins CI clients.

ci_positive_integer() {
    [[ "$2" =~ ^[1-9][0-9]*$ ]] || die "$1 must be a positive integer"
}

ci_build_payload() {
    local build_number="$1"
    curl_auth "$(ci_build_url "$build_number")/api/json?tree=number,queueId,building,result,timestamp,duration,url,actions[parameters[name,value]],artifacts[fileName,relativePath]"
}

ci_artifact_records() {
    local build_number="$1" payload="$2"
    python3 -c '
import json
import sys
import urllib.parse

base = sys.argv[1].rstrip("/")
for artifact in json.load(sys.stdin).get("artifacts", []):
    path = artifact.get("relativePath")
    if isinstance(path, str) and path:
        print("{}\t{}/artifact/{}".format(path, base, urllib.parse.quote(path, safe="/")))
' "$(ci_build_url "$build_number")" <<< "$payload"
}

ci_print_build() {
    local payload="$1"
    python3 -c '
import datetime
import json
import sys

build = json.load(sys.stdin)
values = {}
for action in build.get("actions") or []:
    for parameter in action.get("parameters") or []:
        name = parameter.get("name")
        if name in ("CANDIDATE_PR", "CANDIDATE_REF", "BASELINE_REF"):
            values[name] = parameter.get("value") or "-"
timestamp = build.get("timestamp")
when = "-" if timestamp is None else datetime.datetime.fromtimestamp(
    timestamp / 1000.0, datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
state = "RUNNING" if build.get("building", False) else (build.get("result") or "UNKNOWN")
print("BUILD #{number} {state} queue={queue} {when} pr={pr} candidate={candidate} baseline={baseline} url={url}".format(
    number=build.get("number", "?"), state=state,
    queue=build.get("queueId") if build.get("queueId") is not None else "-",
    when=when, pr=values.get("CANDIDATE_PR", "-"),
    candidate=values.get("CANDIDATE_REF", "-"),
    baseline=values.get("BASELINE_REF", "-"), url=build.get("url", "-")))
' <<< "$payload"
}

ci_artifact_url() {
    local build_number="$1" artifact_path="$2"
    python3 -c 'import sys, urllib.parse; print(sys.argv[1].rstrip("/") + "/artifact/" + urllib.parse.quote(sys.argv[2], safe="/"))' \
        "$(ci_build_url "$build_number")" "$artifact_path"
}

ci_safe_artifact_path() {
    local artifact_path="$1" component
    [[ -n "$artifact_path" && "$artifact_path" != /* && "$artifact_path" != *'//' ]] || return 1
    IFS=/ read -r -a components <<< "$artifact_path"
    for component in "${components[@]}"; do
        [[ -n "$component" && "$component" != . && "$component" != .. ]] || return 1
    done
}

ci_download_artifact() {
    local build_number="$1" artifact_path="$2" destination="$3"
    curl_auth --output "$destination" "$(ci_artifact_url "$build_number" "$artifact_path")"
}

ci_print_failure_summary() {
    local build_number="$1" payload="$2" record artifact_path artifact_url temporary
    while IFS=$'\t' read -r artifact_path artifact_url; do
        [[ "$artifact_path" == ci-failure-summary.txt ]] || continue
        temporary="$(mktemp "${TMPDIR:-/tmp}/gkeyll-ci-failure.XXXXXX")"
        if ci_download_artifact "$build_number" "$artifact_path" "$temporary"; then
            local stage='' message='' result='' key value
            while IFS='=' read -r key value; do
                case "$key" in
                    stage) stage="$value" ;;
                    message) message="$value" ;;
                    result) result="$value" ;;
                esac
            done < "$temporary"
            rm -f "$temporary"
            echo "Failure: ${stage:-unknown}: ${message:-${result:-unknown failure}}"
        else
            rm -f "$temporary"
            echo 'Failure: CI failure summary artifact could not be downloaded.'
        fi
        return
    done < <(ci_artifact_records "$build_number" "$payload")
}

ci_print_regression_failures() {
    local build_number="$1" payload="$2" queryrdb record artifact_path artifact_url temporary layer output found=false
    queryrdb="${GKYL_QUERYRDB:-$SCRIPT_DIR/../../gkylsoft/gkeyll/bin/gkeyll}"
    [[ -x "$queryrdb" ]] || {
        echo "Regression failures: queryrdb is unavailable at $queryrdb"
        return
    }
    while IFS=$'\t' read -r artifact_path artifact_url; do
        [[ "$artifact_path" =~ ^gkylsoft/gkeyll-results/(.*/)?(moments|vlasov|gyrokinetic|pkpm)/regressiondb$ ]] || continue
        layer="${BASH_REMATCH[2]}"
        temporary="$(mktemp "${TMPDIR:-/tmp}/gkeyll-ci-regressiondb.XXXXXX")"
        if ! ci_download_artifact "$build_number" "$artifact_path" "$temporary"; then
            rm -f "$temporary"
            echo "Regression failures: could not download $artifact_path"
            continue
        fi
        if output="$("$queryrdb" queryrdb --db "$temporary" query --id 1 --fail-only 2>/dev/null)"; then
            while IFS= read -r record; do
                if [[ "$record" =~ ^[[:space:]]*[0-9]+[[:space:]]*:[[:space:]]+[^[:space:]]+[[:space:]]+(.+)[[:space:]]+(fail|timeout|compile_fail|no_output|crash)[[:space:]] ]]; then
                    if [[ "$found" == false ]]; then
                        echo 'Regression failures:'
                        found=true
                    fi
                    echo "  $layer/${BASH_REMATCH[1]} [${BASH_REMATCH[2]}]"
                fi
            done <<< "$output"
        else
            echo "Regression failures: queryrdb could not inspect $artifact_path"
        fi
        rm -f "$temporary"
    done < <(ci_artifact_records "$build_number" "$payload")
}

ci_info_command() {
    [[ $# -eq 2 && "$1" == --build ]] || die 'usage: info --build NUMBER'
    ci_positive_integer 'build number' "$2"
    local payload
    payload="$(ci_build_payload "$2")" || die "Build #$2 is not available"
    ci_print_build "$payload"
    ci_print_failure_summary "$2" "$payload"
    ci_print_regression_failures "$2" "$payload"
    echo "Artifacts URL: $(ci_build_url "$2")/artifact/"
    echo 'Artifacts:'
    local artifact_path artifact_url found=false
    while IFS=$'\t' read -r artifact_path artifact_url; do
        echo "  $artifact_path  $artifact_url"
        found=true
    done < <(ci_artifact_records "$2" "$payload")
    [[ "$found" == true ]] || echo '  (none)'
}

ci_artifact_command() {
    local build_number='' mode=list only='' output_dir=''
    while (($#)); do
        case "$1" in
            --build) (($# >= 2)) || die '--build requires a number'; build_number="$2"; shift 2 ;;
            --list) [[ "$mode" == list ]] || die '--list cannot be combined with --fetch'; shift ;;
            --fetch) [[ "$mode" == list ]] || die '--fetch cannot be repeated'; mode=fetch; shift ;;
            --only) (($# >= 2)) || die '--only requires one or more artifact paths'; only="$2"; shift 2 ;;
            --output-dir) (($# >= 2)) || die '--output-dir requires a directory'; output_dir="$2"; shift 2 ;;
            *) die "unknown artifact option: $1" ;;
        esac
    done
    [[ -n "$build_number" ]] || die 'usage: artifact --build NUMBER [--list | --fetch [--only PATH[,PATH...]] [--output-dir DIR]]'
    ci_positive_integer 'build number' "$build_number"
    [[ "$mode" == fetch || -z "$only$output_dir" ]] || die '--only and --output-dir require --fetch'

    local payload
    payload="$(ci_build_payload "$build_number")" || die "Build #$build_number is not available"
    if [[ "$mode" == list ]]; then
        local artifact_path artifact_url found=false
        while IFS=$'\t' read -r artifact_path artifact_url; do
            echo "$artifact_path  $artifact_url"
            found=true
        done < <(ci_artifact_records "$build_number" "$payload")
        [[ "$found" == true ]] || echo 'No archived artifacts.'
        return
    fi

    [[ -n "$output_dir" ]] || output_dir="gkeyll-ci-build-$build_number"
    [[ ! -e "$output_dir" ]] || die "output directory already exists: $output_dir"
    local -a available=() selected=() requested=()
    local artifact_path artifact_url requested_path
    while IFS=$'\t' read -r artifact_path artifact_url; do
        ci_safe_artifact_path "$artifact_path" || die "Jenkins returned unsafe artifact path: $artifact_path"
        available+=("$artifact_path")
    done < <(ci_artifact_records "$build_number" "$payload")
    if [[ -n "$only" ]]; then
        IFS=, read -r -a requested <<< "$only"
        for requested_path in "${requested[@]}"; do
            [[ -n "$requested_path" ]] || die '--only contains an empty artifact path'
            local match=false
            for artifact_path in "${available[@]}"; do
                [[ "$artifact_path" == "$requested_path" ]] || continue
                selected+=("$artifact_path")
                match=true
                break
            done
            [[ "$match" == true ]] || die "artifact is not available in build #$build_number: $requested_path"
        done
    else
        selected=("${available[@]}")
    fi
    ((${#selected[@]})) || die "Build #$build_number has no archived artifacts"

    mkdir -p "$output_dir"
    for artifact_path in "${selected[@]}"; do
        local destination="$output_dir/$artifact_path"
        mkdir -p "$(dirname "$destination")"
        if ! ci_download_artifact "$build_number" "$artifact_path" "$destination"; then
            rm -rf "$output_dir"
            die "could not download artifact: $artifact_path"
        fi
        echo "Downloaded $artifact_path"
    done
    echo "Artifacts downloaded to: $output_dir"
}
