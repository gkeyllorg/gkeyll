#!/bin/sh
# Run every test, retaining complete diagnostics and each process's exit code.
set -u
command -v valgrind >/dev/null 2>&1 || { echo 'Valgrind is not installed' >&2; exit 1; }
status=0
for unit do
    result=0
    log="${unit}_val_err"
    if [ -n "${GKYL_VALGRIND_LOG_DIR:-}" ]; then
        log="$GKYL_VALGRIND_LOG_DIR/$(basename "$unit").log"
    fi
    case "$unit" in
        /*) executable=$unit ;;
        *) executable=./$unit ;;
    esac
    valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all \
        --errors-for-leak-kinds=definite,indirect,possible --error-exitcode=99 \
        "$executable" -E > "$log" 2>&1 || result=$?
    printf '%s\n' "$result" > "$log.exit"
    if [ "$result" -eq 0 ]; then
        printf 'PASS valgrind: %s\n' "$unit"
    else
        printf 'FAIL valgrind: %s (exit %s)\n' "$unit" "$result"
        cat "$log"
        status=1
    fi
done
exit "$status"
