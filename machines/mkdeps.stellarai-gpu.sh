. "$(dirname "$0")/module_load.stellarai-gpu.sh"

: "${PREFIX:=$HOME/gkylsoft}"

cd install-deps
./mkdeps.sh --build-superlu=yes --prefix=$PREFIX --build-cudss=yes --build-luajit=yes --build-adas=yes
