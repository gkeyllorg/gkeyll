. "$(dirname "$0")/module_load.perlmutter-cpu.sh"

: "${PREFIX:=$HOME/gkylsoft}"

cd install-deps
./mkdeps.sh --build-superlu=yes --build-luajit=yes --prefix=$PREFIX
