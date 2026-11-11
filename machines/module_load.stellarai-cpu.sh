# Source this file to load modules needed by Gkeyll builds and simulations
# on Princeton Stellar AI's CPU partition.

# Some module files source /usr/share/Modules/init/bash, which expects
# PS1 even in a non-interactive shell.
: "${PS1:=}"

# CI batch jobs start without login profiles or inherited module state.
if ! type module >/dev/null 2>&1; then
  . /usr/share/Modules/init/bash
fi

module purge
module load intel-classic/2023.2.3
module load openmpi/intel-classic-2023.2.3/5.0.10
module load openblas/0.3.x
# Keep Anaconda MPI out of the build and simulation runtime.
