# Source this file load modules needed by Gkeyll builds and sims
# on Princeton StellarAI (GPU builds/sims).

# Some module files source /usr/share/Modules/init/bash, which expects
# PS1 even in a non-interactive shell. 
: "${PS1:=}"

module purge
module load nvhpc/26.5 # For NCCL.
module load cudatoolkit/13.3
module load openmpi/nvhpc-26.5/5.0.10
module load openblas/0.3.x
module load anaconda3/2026.7 # To download ADAS.
