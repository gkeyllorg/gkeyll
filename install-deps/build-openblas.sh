#!/bin/bash

source ./build-opts.sh

# Edit to suite your system
PREFIX=$GKYLSOFT/OpenBLAS-0.3.30
# Location where dependency sources will be downloaded
DEP_SOURCES=$GKYLSOFT/dep_src/

mkdir -p $DEP_SOURCES
cd $DEP_SOURCES

if [ "$DOWNLOAD_PKGS" = "yes" ]
then
    echo "Downloading OpenBLAS .."
    # delete old checkout and builds
    rm -rf OpenBLAS-*
    curl -L https://github.com/xianyi/OpenBLAS/releases/download/v0.3.30/OpenBLAS-0.3.30.tar.gz > OpenBLAS-0.3.30.tar.gz
fi

if [ "$BUILD_PKGS" = "yes" ]
then
    echo "Building OpenBLAS .."
    gunzip -f OpenBLAS-0.3.30.tar.gz
    tar xvf OpenBLAS-0.3.30.tar
    cd OpenBLAS-0.3.30
    make USE_OPENMP=0 NUM_THREADS=1 NO_FORTRAN=1 -j 32
    make USE_OPENMP=0 NUM_THREADS=1 NO_FORTRAN=1 install PREFIX=$PREFIX -j 32

    # soft-link 
    ln -sfn $PREFIX $GKYLSOFT/OpenBLAS
    # remove shared libraries
    rm -rf $GKYLSOFT/OpenBLAS/lib/*.so*
fi
