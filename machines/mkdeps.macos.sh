cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
export MACOSX_DEPLOYMENT_TARGET=14.7
./mkdeps.sh --build-openblas=no --build-superlu=yes --build-luajit=yes --prefix=$PREFIX --build-adas=yes
