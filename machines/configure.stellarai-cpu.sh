. "$(dirname "$0")/module_load.stellarai-cpu.sh"

: "${PREFIX:=$HOME/gkylsoft}"

./configure CC=cc --prefix=$PREFIX --lapack-inc=/usr/include/openblas/ --lapack-lib=/usr/lib64/ --use-mpi=yes --mpi-inc=$MPI_HOME/include --mpi-lib=$MPI_HOME/lib64 --use-lua=yes
