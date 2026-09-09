module unload darshan
module unload craype-accel-nvidia
module unload nccl
module unload cudatoolkit
module unload gpu
module load cray-mpich/9.0.1
module load cray-libsci/25.09.0
module load craype-accel-host

: "${PREFIX:=$HOME/gkylsoft}"

cd install-deps
./mkdeps.sh --build-superlu=yes --build-luajit=yes --prefix=$PREFIX
