module load intel/2022.2.0
module load openmpi/intel-2022.0/4.1.8
module load openblas/0.3.x

: "${PREFIX:=$HOME/gkylsoft}"
: "${ARCH_FLAGS:=-march=cascadelake}" # The Intel compute nodes are Cascade Lake.

./configure CC=cc ARCH_FLAGS=$ARCH_FLAGS --prefix=$PREFIX --lapack-inc=/usr/include/openblas/ --lapack-lib=/usr/lib64/ --use-mpi=yes --mpi-inc=$MPI_HOME/include --mpi-lib=$MPI_HOME/lib64 --use-lua=yes
