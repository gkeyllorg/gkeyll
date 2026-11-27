#!/bin/bash -l
# Submit site payloads and recover/cancel this build's allocations on every exit.
set -euo pipefail

write_status() {
    local file="$1" job_id="$2" state="$3"
    printf 'job_id=%s\nstate=%s\nexit_code=%s\nelapsed_seconds=%s\n' \
        "$job_id" "$state" "${4:-}" "${5:-}" > "$file.tmp"
    mv -f "$file.tmp" "$file"
}

cleanup_job() {
    local id_file="$1" owner_file="$1.owner" status_file="${1%-id.txt}-status.txt"
    local owner rows id name directory job_id='' cancelled=false attempt empty_queries=0
    [[ -f "$owner_file" ]] || {
        echo "ERROR: No submission ownership record for $id_file; cannot safely cancel it" >&2
        return 1
    }
    owner="$(<"$owner_file")"
    [[ "$owner" =~ ^gkyl-ci-[0-9a-f]{32}$ ]] || {
        echo "ERROR: Invalid submission ownership record: $owner_file" >&2
        return 1
    }
    if [[ -f "$id_file" ]]; then job_id="$(<"$id_file")"; fi
    # Job IDs and status files alone are insufficient: submission may have been
    # interrupted before an ID was saved, or a status write may be incomplete.
    # Match the unique name written BEFORE sbatch, the user and the workspace.
    for attempt in {1..12}; do
        if ! rows="$(squeue --noheader --user="$(id -un)" --name="$owner" --format='%i|%j|%Z')"; then
            echo "ERROR: Cannot verify Slurm cleanup for $owner; preserving ownership records" >&2
            return 1
        fi
        if [[ -z "$rows" ]]; then
            # With no returned ID, allow a brief window for an in-flight submit
            # RPC to become visible after Jenkins has killed its client.
            if [[ -z "$job_id" ]] && (( ++empty_queries < 3 )); then
                sleep 1
                continue
            fi
            if [[ "$cancelled" == true ]]; then
                write_status "$status_file" "$job_id" CANCELLED_BY_JENKINS
            fi
            return 0
        fi
        while IFS='|' read -r id name directory; do
            [[ "$id" =~ ^[0-9]+(_[0-9]+)?$ && "$name" == "$owner" && "$directory" == "$CI_WORKSPACE" ]] || {
                echo "ERROR: Unexpected scheduler ownership response for $owner; refusing cancellation" >&2
                return 1
            }
            if [[ -n "$job_id" && "$job_id" != "$id" ]]; then
                echo "Recovering Slurm job ID from verified ownership ($id_file contained '$job_id')" >&2
            fi
            job_id="$id"
            printf '%s\n' "$job_id" > "$id_file"
            echo "Cancelling unfinished Slurm job $job_id ($owner)"
            # A job may finish between squeue and scancel. The next successful
            # queue query, rather than scancel's exit code, confirms cleanup.
            scancel "$job_id" || true
            cancelled=true
            CLEANUP_CANCELLED=true
        done <<< "$rows"
        sleep 5
    done
    echo "ERROR: Slurm job $job_id is still active after cancellation; preserving its files" >&2
    return 1
}

if [[ "${1:-}" == cleanup ]]; then
    : "${CI_WORKSPACE:?}"
    result=0
    # An interrupted submission can leave only the ownership record. Also find
    # IDs without ownership records and report them instead of silently skipping.
    shopt -s nullglob
    for owner_file in "$CI_WORKSPACE"/slurm-*-job-id.txt.owner; do
        cleanup_job "${owner_file%.owner}" || result=1
    done
    for id_file in "$CI_WORKSPACE"/slurm-*-job-id.txt; do
        [[ -f "$id_file.owner" ]] || { cleanup_job "$id_file" || result=1; }
    done
    exit "$result"
fi
(($# == 0)) || { echo 'Usage: slurm-submit.sh [cleanup]' >&2; exit 2; }
: "${CI_WORKSPACE:?}" "${CI_SLURM_JOB_ID_FILE:?}" "${CI_SLURM_JOB_STATUS_FILE:?}"

on_exit() {
    local result=$?
    trap - EXIT
    CLEANUP_CANCELLED=false
    if [[ -f "$CI_SLURM_JOB_ID_FILE.owner" ]]; then
        cleanup_job "$CI_SLURM_JOB_ID_FILE" || result=1
    fi
    # Accounting must not produce a passing result while its allocation remains
    # active. Cancellation on an otherwise successful exit is still an error.
    if [[ "$CLEANUP_CANCELLED" == true && "$result" == 0 ]]; then result=1; fi
    exit "$result"
}
trap on_exit EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
trap 'exit 129' HUP

# Read the frozen configuration; never source a profile from candidate code.
resource_lines="$(python3 "$CI_TRUSTED_DIR/ci/jenkins/slurm_config.py" allocation "$CI_SLURM_CONFIG" "$CI_SLURM_LANE")"
# Each submission gets a unique scheduler identity, persisted before sbatch can
# accept work. The Pipeline finalizer can recover it even after SIGKILL.
if [[ -f "$CI_SLURM_JOB_ID_FILE.owner" ]]; then cleanup_job "$CI_SLURM_JOB_ID_FILE"; fi
rm -f "$CI_SLURM_JOB_ID_FILE" "$CI_SLURM_JOB_STATUS_FILE"
owner="gkyl-ci-$(python3 -c 'import uuid; print(uuid.uuid4().hex)')"
printf '%s\n' "$owner" > "$CI_SLURM_JOB_ID_FILE.owner"
sbatch_args=(--parsable "--job-name=$owner")
while IFS= read -r arg; do sbatch_args+=("$arg"); done <<< "$resource_lines"
sbatch_args+=("--output=${CI_WORKSPACE}/${CI_SLURM_OUTPUT}" "--chdir=${CI_WORKSPACE}" --export=ALL)
submission="$(sbatch "${sbatch_args[@]}" "$CI_SLURM_PAYLOAD")"
job_id="${submission%%;*}"
if [[ ! "$job_id" =~ ^[0-9]+(_[0-9]+)?$ ]]; then
    echo "ERROR: Could not parse Slurm job ID from: $submission" >&2
    exit 1
fi
printf '%s\n' "$job_id" > "$CI_SLURM_JOB_ID_FILE"
echo "Submitted ${CI_SLURM_JOB_KIND} Slurm job $job_id"

previous_state=''
while :; do
    state="$(squeue --noheader --jobs="$job_id" --format='%T')"
    if [[ -z "$state" ]]; then break; fi
    if [[ "$state" != "$previous_state" ]]; then
        echo "${CI_SLURM_JOB_KIND} Slurm job $job_id state: $state"
        previous_state="$state"
    fi
    sleep 5
done

# Accounting may lag the queue. Only a terminal allocation record is a result.
for attempt in {1..6}; do
    record="$(sacct --noheader --parsable2 --allocations --jobs="$job_id" \
        --format=State,ExitCode,ElapsedRaw)" || record=''
    record="${record%%$'\n'*}"
    IFS='|' read -r final_state exit_code elapsed_seconds <<< "$record"
    case "$final_state" in
        COMPLETED|CANCELLED*|FAILED|TIMEOUT|NODE_FAIL|PREEMPTED|BOOT_FAIL|DEADLINE|OUT_OF_MEMORY|REVOKED)
            write_status "$CI_SLURM_JOB_STATUS_FILE" "$job_id" "$final_state" "$exit_code" "$elapsed_seconds"
            echo "${CI_SLURM_JOB_KIND} Slurm job $job_id finished: state=$final_state exit_code=$exit_code elapsed_seconds=$elapsed_seconds"
            [[ "$final_state" == COMPLETED && "$exit_code" == '0:0' ]] && exit 0
            exit 1
            ;;
    esac
    echo "Waiting for terminal Slurm accounting record for job $job_id ($attempt/6)"
    sleep 5
done
write_status "$CI_SLURM_JOB_STATUS_FILE" "$job_id" UNKNOWN
echo "ERROR: No terminal Slurm accounting record found for job $job_id" >&2
exit 1
