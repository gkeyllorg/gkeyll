module load intel/2022.2.0
module load openmpi/intel-2022.0/4.1.8
module load openblas/0.3.x
module load anaconda3/2026.7

: "${PREFIX:=$HOME/gkylsoft}"

cd install-deps
./mkdeps.sh --build-superlu=yes --prefix=$PREFIX --build-luajit=yes --build-adas=yes
