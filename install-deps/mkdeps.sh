#!/bin/sh

# Defaults
PREFIX=$HOME/gkylsoft

# default build options
CC=gcc
CXX=g++
FC=gfortran
MPICC=mpicc
MPICXX=mpicxx

# by default, do not build anything
BUILD_OPENBLAS=no
BUILD_SUPERLU=no
BUILD_OPENMPI=no
BUILD_LUAJIT=no
BUILD_TCC=no
BUILD_CUDSS=no
BUILD_EIRENE=no

# by default, download as well as build packages
DOWNLOAD_PKGS=yes
BUILD_PKGS=yes

# ----------------------------------------------------------------------------
# Function definitions
# ----------------------------------------------------------------------------

show_help() {
cat <<EOF

./mkdeps.sh CC=cc CXX=cxx FC=gfortran

Build GkylZero dependencies

CC 
CXX                         C and C++ compilers to use
FC                          Fortran compiler to use (only gfortran is supported)
MPICC                       
MPICXX                      MPI C and C++ compilers to use

-h
--help                      This help.

--download                  [yes] Download packages?
--build                     [yes] Build packages?

--prefix=DIR                Prefix where dependencies should be installed.
                            Default is $HOME/gkylsoft

The following flags specify the libraries to build.

--build-openblas            [no] Should we build OpenBLAS?
--build-superlu             [no] Should we build SuperLU (serial)
--build-openmpi             [no] Should we build OpenMPI?
--build-luajit              [no] Should we build LuaJIT?
--build-cudss               [no] Should we build cuDSS?
--build-tcc                 [no] Should we build tcc?
--build-adas                [no] Should we download ADAS data? (uses python, needs numpy)
--build-eirene              [no] Should we build eirene? (includes installing SOLPS-ITER)

EOF
}

# Helper functions

find_program() {
   prog=`command -v "$1" 2>/dev/null`
   if [ -n "$prog" ]
   then
      dirname "$prog"
   fi
}

die() {
   echo "$*"
   echo
   echo "Dependency builds failed."
   echo
   exit 1
}

# ----------------------------------------------------------------------------
# MAIN PROGRAM
# ----------------------------------------------------------------------------

# Parse options

while [ -n "$1" ]
do
   value="`echo $1 | sed 's/[^=]*.\(.*\)/\1/'`"
   key="`echo $1 | sed 's/=.*//'`"
   if `echo "$value" | grep "~" >/dev/null 2>/dev/null`
   then
      echo
      echo '*WARNING*: the "~" sign is not expanded in flags.'
      echo 'If you mean the home directory, use $HOME instead.'
      echo
   fi
   case "$key" in
   -h)
      show_help
      exit 0
      ;;
   --help)
      show_help
      exit 0
      ;;
   --build)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_PKGS="$value"
      ;;      
   --download)
      [ -n "$value" ] || die "Missing value in flag $key."
      DOWNLOAD_PKGS="$value"
      ;;
   CC)
      [ -n "$value" ] || die "Missing value in flag $key."
      CC="$value"
      ;;
   CXX)
      [ -n "$value" ] || die "Missing value in flag $key."
      CXX="$value"
      ;;
   MPICC)
      [ -n "$value" ] || die "Missing value in flag $key."
      MPICC="$value"
      ;;
   MPICXX)
      [ -n "$value" ] || die "Missing value in flag $key."
      MPICXX="$value"
      ;;
   FC)
      [ -n "$value" ] || die "Missing value in flag $key."
      FC="$value"
      ;;   
   --prefix)
      [ -n "$value" ] || die "Missing value in flag $key."
      PREFIX="$value"
      ;;
   --build-openblas)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_OPENBLAS="$value"
      ;;
   --build-superlu)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_SUPERLU="$value"
      ;;
   --build-openmpi)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_OPENMPI="$value"
      MPICC=$PREFIX/openmpi/bin/mpicc
      MPICXX=$PREFIX/openmpi/bin/mpicxx
      ;;   
   --build-luajit)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_LUAJIT="$value"
      ;;
   --build-tcc)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_TCC="$value"
      ;;
   --build-cudss)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_CUDSS="$value"
      ;;   
   --build-adas)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_ADAS="$value"
      ;;
   --build-eirene)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_EIRENE="$value"
      ;;
   *)
      die "Error: Unknown flag: $1"
      ;;
   esac
   shift
done

CMAKE_SUPERLU_DIST_GPU=OFF
# Set package options
if [ "$BUILD_SUPERLU_DIST_GPU" = "yes" ]
then
    CMAKE_SUPERLU_DIST_GPU=ON
fi

# Write out build options for scripts to use
cat <<EOF1 > build-opts.sh
#!/bin/bash

# Generated automatically! Do not edit

# Download/Build options
DOWNLOAD_PKGS=$DOWNLOAD_PKGS
BUILD_PKGS=$BUILD_PKGS

# Installation directory
GKYLSOFT=$PREFIX
# Various compilers
CC=$CC
CXX=$CXX
MPICC=$MPICC
MPICXX=$MPICXX
FC=gfortran

# Package options
CMAKE_SUPERLU_DIST_GPU=$CMAKE_SUPERLU_DIST_GPU

EOF1

build_openblas() {
    if [ "$BUILD_OPENBLAS" = "yes" ]
    then    
	echo "Building OpenBLAS"
	./build-openblas.sh
    fi
}

build_superlu() {
    if [ "$BUILD_SUPERLU" = "yes" ]
    then    
	echo "Building SUPERLU"
	./build-superlu.sh 
    fi
}

build_openmpi() {
    if [ "$BUILD_OPENMPI" = "yes" ]
    then    
	echo "Building OpenMPI"
	./build-openmpi.sh
    fi
}

build_luajit() {
    if [ "$BUILD_LUAJIT" = "yes" ]
    then    
	echo "Building Luajit"
	./build-luajit.sh
    fi
}

build_tcc() {
    if [ "$BUILD_TCC" = "yes" ]
    then    
	echo "Building TCC"
	./build-tcc.sh
    fi
}

build_cudss() {
    if [ "$BUILD_CUDSS" = "yes" ]
    then    
	echo "Building cuDSS"
	./build-cudss.sh
    fi
}

build_adas() {
    if [ "$BUILD_ADAS" = "yes" ]
    then
	echo "Downloading ADAS data for neutral reactions"
	./download-adas.sh
    fi
}

build_eirene() {
    if [ "$BUILD_EIRENE" = "yes" ]
    then
	echo "Building EIRNE Coupling for neutral coupling"
	./build-eirene.sh
    fi
}

echo "Installations will be in  $PREFIX"

build_openmpi
build_luajit
build_tcc
build_openblas
build_superlu
build_cudss
build_adas
build_eirene
