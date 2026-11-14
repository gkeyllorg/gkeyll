. "$(dirname "$0")/module_load.stellarai-gpu.sh"

: "${PREFIX:=$HOME/gkylsoft}"
: "${ARCH_FLAGS:=-march=emeraldrapids}" # GPU compute nodes are Intel Xeon Platinum 8568Y+.

./configure CC=nvcc ARCH_FLAGS=$ARCH_FLAGS CUDA_ARCH=100 --prefix=$PREFIX --lapack-inc=/usr/include/openblas/ --lapack-lib=/usr/lib64/ --cudamath-lib=$CUDA_HOME/lib64 --use-mpi=yes --mpi-inc=$MPI_HOME/include --mpi-lib=$MPI_HOME/lib64 --use-nccl=yes --nccl-inc=$NVHPC_ROOT/comm_libs/nccl/include --nccl-lib=$NVHPC_ROOT/comm_libs/nccl/lib --use-cudss=yes --use-lua=yes;

