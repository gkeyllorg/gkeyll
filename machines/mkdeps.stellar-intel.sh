. "$(dirname "$0")/module_load.stellar-intel.sh"

: "${PREFIX:=$HOME/gkylsoft}"

cd install-deps
./mkdeps.sh --build-superlu=yes --prefix=$PREFIX --build-luajit=yes --build-adas=yes
