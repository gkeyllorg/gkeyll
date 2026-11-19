module load cuda
module load nvidia_math
module load nccl
: "${PREFIX:=$HOME/gkylsoft}"
./configure CC=nvcc ARCH_FLAGS="-march=native" CUDA_ARCH=90 --prefix=$PREFIX --lapack-inc=$PREFIX/OpenBLAS/include --lapack-lib=$PREFIX/OpenBLAS/lib/libopenblas.a --superlu-inc=$PREFIX/superlu/include --superlu-lib=$PREFIX/superlu/lib/libsuperlu.a --cudamath-libdir=$TACC_NVIDIA_MATH_LIB --use-mpi=yes --mpi-inc=$TACC_MPI_DIR/include --mpi-lib=$TACC_MPI_DIR/lib --use-nccl=yes --nccl-inc=$TACC_NCCL_DIR/include --nccl-lib=$TACC_NCCL_DIR/lib --use-lua=yes;


