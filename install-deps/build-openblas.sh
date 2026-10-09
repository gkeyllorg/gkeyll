#!/bin/bash
set -euo pipefail

source ./build-opts.sh

case "${GKEYLL_USE_VALGRIND:-0}" in
    0|1) ;;
    *) echo 'GKEYLL_USE_VALGRIND must be 0 or 1' >&2; exit 1 ;;
esac
options=(CC="$CC" USE_OPENMP=0 NUM_THREADS=1 NO_FORTRAN=1)
if [ "${GKEYLL_USE_VALGRIND:-0}" = 1 ]; then
    case "$($CC -dumpmachine)" in
        x86_64*|amd64*) options+=(TARGET=GENERIC NO_AVX=1 NO_AVX2=1 NO_AVX512=1) ;;
        aarch64*|arm64*) options+=(TARGET=ARMV8) ;;
        *) echo 'No Valgrind OpenBLAS profile for this compiler target' >&2; exit 1 ;;
    esac
    options+=(DYNAMIC_ARCH=0)
fi

# Default builds retain OpenBLAS CPU detection, including AVX512 where supported.
PREFIX=$GKYLSOFT/OpenBLAS-0.3.30
# Location where dependency sources will be downloaded
DEP_SOURCES=$GKYLSOFT/dep_src/

mkdir -p "$DEP_SOURCES"
cd "$DEP_SOURCES"

if [ "$DOWNLOAD_PKGS" = "yes" ]
then
    echo "Downloading OpenBLAS .."
    # delete old checkout and builds
    rm -rf OpenBLAS-*
    curl -fL https://github.com/xianyi/OpenBLAS/releases/download/v0.3.30/OpenBLAS-0.3.30.tar.gz > OpenBLAS-0.3.30.tar.gz
fi

if [ "$BUILD_PKGS" = "yes" ]
then
    echo "Building OpenBLAS .."
    if [ -f OpenBLAS-0.3.30.tar.gz ]; then
        gunzip -f OpenBLAS-0.3.30.tar.gz
    fi
    # A fresh extraction prevents stale objects when switching build profiles.
    rm -rf OpenBLAS-0.3.30
    tar xvf OpenBLAS-0.3.30.tar
    cd OpenBLAS-0.3.30
    jobs=${NPROC:-}
    if [ -z "$jobs" ]; then
        processors=$(nproc 2>/dev/null || sysctl -n hw.physicalcpu)
        jobs=$((processors > 1 ? processors / 2 : 1))
    fi
    make "${options[@]}" -j "$jobs"
    make "${options[@]}" install PREFIX="$PREFIX" -j "$jobs"

    # soft-link 
    ln -sfn $PREFIX $GKYLSOFT/OpenBLAS
    # remove shared libraries
    rm -rf $GKYLSOFT/OpenBLAS/lib/*.so*
fi
