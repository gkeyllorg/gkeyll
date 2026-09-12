module unload darshan
module unload craype-accel-nvidia
module unload nccl
module unload cudatoolkit
module unload gpu
module load cray-mpich/9.0.1
module load cray-libsci/25.09.0
module load craype-accel-host

: "${PREFIX:=$HOME/gkylsoft}"

./configure CC=cc --prefix=$PREFIX --lapack-lib-name=sci_gnu_mpi --lapack-inc=$CRAY_LIBSCI_PREFIX/include --lapack-lib=$CRAY_LIBSCI_PREFIX/lib --use-mpi=yes --mpi-inc=$CRAY_MPICH_DIR/include --mpi-lib=$CRAY_MPICH_DIR/lib --use-lua=yes;


