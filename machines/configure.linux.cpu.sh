: "${PREFIX:=$HOME/gkylsoft}"
: "${MPI_HOME:=$HOME/gkylsoft/openmpi}"
./configure CC=clang --prefix=$PREFIX --use-lua=yes --use-mpi=yes --mpi-inc=$MPI_HOME/include
