: "${PREFIX:=$HOME/gkylsoft}"
: "${MPI_HOME:=$HOME/gkylsoft/openmpi}"
# -fno-slp-vectorize works around a clang-15 SLP-vectorizer crash (segfault in
# the optimizer) on Gkeyll's large generated p3 kernels at -O3 -march=native.
./configure CC=clang ARCH_FLAGS="-march=native -fno-slp-vectorize" --prefix=$PREFIX --use-lua=yes --use-mpi=yes --mpi-inc=$MPI_HOME/include
