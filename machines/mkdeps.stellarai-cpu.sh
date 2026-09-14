module load intel-classic/2023.2.3
module load openmpi/intel-classic-2023.2.3/5.0.10
module load openblas/0.3.x
module load anaconda3/2026.7

: "${PREFIX:=/scratch/gpfs/CSD/manaurer/gkeyll/code/gkeyll_cpu0/gkylsoft}"

cd install-deps
./mkdeps.sh --build-superlu=yes --prefix=$PREFIX --build-luajit=yes --build-adas=yes
