. "$(dirname "$0")/module_load.stellarai-cpu.sh"

: "${PREFIX:=$HOME/gkylsoft}"

# NumPy is needed to download ADAS; this script does not launch MPI.
module load anaconda3/2026.7

cd install-deps
./mkdeps.sh --build-superlu=yes --prefix=$PREFIX --build-luajit=yes --build-adas=yes
