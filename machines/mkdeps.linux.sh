cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
./mkdeps.sh --prefix=$PREFIX --build-adas=yes --build-openblas=yes --build-superlu=yes --build-openmpi=yes --build-luajit=yes --build-adas=yes CC=clang CXX=clang++ 