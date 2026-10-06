# Source this file to load modules needed by Gkeyll builds and simulations
# on Princeton Stellar AI's CPU partition.

# Some module files source /usr/share/Modules/init/bash, which expects
# PS1 even in a non-interactive shell.
: "${PS1:=}"

module purge
module load intel-classic/2023.2.3
module load openmpi/intel-classic-2023.2.3/5.0.10
module load openblas/0.3.x
module load anaconda3/2026.7 # To download ADAS.
