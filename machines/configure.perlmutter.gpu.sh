. "$(dirname "$0")/module_load.perlmutter-gpu.sh"

: "${PREFIX:=$HOME/gkylsoft}"

./configure CC=nvcc ARCH_FLAGS="-march=native" CUDA_ARCH=80 --prefix=$PREFIX --lapack-lib-name=sci_gnu --lapack-inc=$CRAY_LIBSCI_PREFIX/include --lapack-lib=$CRAY_LIBSCI_PREFIX/lib --cudamath-lib=/opt/nvidia/hpc_sdk/Linux_x86_64/25.9/math_libs/13.0/lib64 --use-mpi=yes --mpi-inc=$CRAY_MPICH_DIR/include --mpi-lib=$CRAY_MPICH_DIR/lib --use-nccl=yes --nccl-inc=$NCCL_DIR/include --nccl-lib=$NCCL_DIR/lib --use-lua=yes --use-cudss=yes;
