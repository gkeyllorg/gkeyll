. "$(dirname "$0")/module_load.della-gpu.sh"
cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
./mkdeps.sh --build-openblas=yes --build-superlu=yes --build-luajit=yes --prefix=$PREFIX
