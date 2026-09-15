. "$(dirname "$0")/module_load.perlmutter-gpu.sh"

: "${PREFIX:=$HOME/gkylsoft}"

cd install-deps
./mkdeps.sh --build-superlu=yes --build-cudss=yes --prefix=$PREFIX --build-luajit=yes MPICC=mpicc  MPICXX=mpicxx --build-adas=yes
