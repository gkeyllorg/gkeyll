. "$(dirname "$0")/module_load.perlmutter-cpu.sh"

: "${PREFIX:=$HOME/gkylsoft}"

./configure CC=cc --prefix=$PREFIX --lapack-lib-name=sci_gnu_mpi --lapack-inc=$CRAY_LIBSCI_PREFIX/include --lapack-lib=$CRAY_LIBSCI_PREFIX/lib --use-mpi=yes --mpi-inc=$CRAY_MPICH_DIR/include --mpi-lib=$CRAY_MPICH_DIR/lib --use-lua=yes;

