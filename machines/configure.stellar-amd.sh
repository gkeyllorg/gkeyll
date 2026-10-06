. "$(dirname "$0")/module_load.stellar-amd.sh"

: "${PREFIX:=$HOME/gkylsoft}"

./configure CC=nvcc ARCH_FLAGS="-march=native" CUDA_ARCH=80 --prefix=$PREFIX --lapack-inc=/usr/include/openblas/ --lapack-lib=/usr/lib64/ --cudamath-lib=$CUDA_HOME/lib64 --use-mpi=yes --mpi-inc=$MPI_HOME/include --mpi-lib=$MPI_HOME/lib64 --use-nccl=yes --nccl-inc=$NVHPC_ROOT/comm_libs/nccl/include --nccl-lib=$NVHPC_ROOT/comm_libs/nccl/lib --use-cudss=yes --use-lua=yes;
