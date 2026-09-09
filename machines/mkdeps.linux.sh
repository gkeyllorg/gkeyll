cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
./mkdeps.sh --build-openblas=yes --build-superlu=yes --build-openmpi=yes --build-luajit=yes --prefix=$PREFIX
