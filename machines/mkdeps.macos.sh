cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
# Target the macOS version of the build host: the dependencies are built locally.
export MACOSX_DEPLOYMENT_TARGET=$(sw_vers -productVersion | cut -d. -f1,2)

./mkdeps.sh --build-openblas=no --build-superlu=yes --build-openmpi=yes --build-luajit=yes --build-adas=yes --prefix=$PREFIX
