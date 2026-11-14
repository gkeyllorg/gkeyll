# Source this file load modules needed by Gkeyll builds and sims
# on Princeton Stellar's Intel partition.

# Some module files source /usr/share/Modules/init/bash, which expects
# PS1 even in a non-interactive shell. 
: "${PS1:=}"

module purge
module load intel/2022.2.0
module load openmpi/intel-2022.0/4.1.8
module load openblas/0.3.x
module load anaconda3/2026.7
