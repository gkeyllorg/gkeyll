#!/usr/bin/env bash
# Shared build operation. A site's build Slurm script may prepare its environment
# and then invoke this helper; the helper never chooses a process launcher.
set -euo pipefail
: "${CI_BUILD_SOURCE:?}" "${CI_BUILD_PREFIX:?}" "${CI_BUILD_ACTION:?}" "${CI_WORKSPACE:?}"
cd "$CI_BUILD_SOURCE"
export PREFIX="$CI_BUILD_PREFIX"
. "$SLURM_CI_ENV_SCRIPT"
metric=''
case "$CI_BUILD_ACTION" in
    deps) bash "$SLURM_CI_MKDEPS_SCRIPT" ;;
    configure)
        bash "$SLURM_CI_CONFIGURE_SCRIPT"
        if grep -Eq '^[[:space:]]*CC[[:space:]]*=[[:space:]]*(nvcc|nvc)([[:space:]]|$)' config.mak; then
            # Apply strict FP to both candidate and baseline GPU builds.
            printf '\nUSE_STRICT_FP=1\nUSE_DETERMINISTIC_REDUCTIONS_GPU=1\n' >> config.mak
        fi
        ;;
    unit|install)
        started="$(date +%s)"
        make -j"$SLURM_CI_BUILD_JOBS" "$CI_BUILD_ACTION"
        if [[ "$CI_BUILD_ACTION" == unit ]]; then metric=unit-build; else metric=install; fi
        ;;
    serial-compile|parallel-compile)
        driver="$CI_BUILD_PREFIX/gkeyll/bin/gkeyll"
        configure_args=(--source-dir "$CI_BUILD_SOURCE" --prefix "$CI_BUILD_PREFIX")
        run_args=(-c)
        if [[ "$CI_BUILD_ACTION" == parallel-compile ]]; then
            configure_args+=(--results-subdir parallel-c-4 --mpiexec "$SLURM_CI_MPIEXEC" --mpi-args "${SLURM_CI_MPI_ARGS:-}")
            run_args+=(--parallel)
        fi
        "$driver" runregression configure "${configure_args[@]}"
        started="$(date +%s)"
        "$driver" runregression run "${run_args[@]}" compile
        metric="${CI_BUILD_ACTION%-compile}-c-compile"
        ;;
    *) echo "Unknown build action: $CI_BUILD_ACTION" >&2; exit 2 ;;
esac
if [[ -n "$metric" ]]; then
    elapsed="$(( $(date +%s) - started ))"
    printf '%s\n' "$elapsed" > "$CI_WORKSPACE/$CI_BUILD_LABEL-$metric-seconds.txt"
    echo "$CI_BUILD_LABEL $CI_BUILD_ACTION: $elapsed seconds"
fi
