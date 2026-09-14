module load intel-classic/2023.2.3
module load openmpi/intel-classic-2023.2.3/5.0.10
module load openblas/0.3.x

: "${PREFIX:=/scratch/gpfs/CSD/manaurer/gkeyll/code/gkeyll_cpu0/gkylsoft}"

./configure CC=cc --prefix=$PREFIX --lapack-inc=/usr/include/openblas/ --lapack-lib=/usr/lib64/ --use-mpi=yes --mpi-inc=$MPI_HOME/include --mpi-lib=$MPI_HOME/lib64 --use-lua=yes

