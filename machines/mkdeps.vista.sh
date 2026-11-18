module load cuda
module load nvidia_math
module load nccl
cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
./mkdeps.sh --build-openblas=yes --build-superlu=yes --prefix=$PREFIX --build-luajit=yes
