#!/bin/bash -l

#.Declare a name for this job, preferably with 16 or fewer characters.
#SBATCH -J <Job Name>

#.Request the queue. See stellar docs online.
#SBATCH --qos pppl-short

#.Number of nodes to request (Stellar has 96 cores and 2 GPUs per node).
#SBATCH -N 1

#.Total number of cores/tasks/MPI processes.
#SBATCH --tasks-per-node=2

#.Request wall time
#SBATCH -t 00:30:00

#.Specify GPUs per node (stellar-amd has 2):
#SBATCH --gres=gpu:2

#.Mail is sent to you when the job starts and when it terminates or aborts.
#SBATCH --mail-user=<your email>
#SBATCH --mail-type=END,FAIL,REQUEUE

#.Specify name format of output file.
#SBATCH -o slurm-%j.out

#.Load modules (this must match those in the machines/configure script).
. ./machines/module_load.stellar-amd.sh

#.Run the rt_gk_sheath_2x2v_p1 executable using 1 GPU along x (-c 1) and 2
#.GPUs along the field line (-d 2). See './rt_gk_sheath_2x2v_p1 -h' for
#.more details/options on decomposition. It also assumes the executable is
#.in the present directory. If it isn't, change `./` to point to the
#.directory containing the executable.

echo "mpirun -np 2 ./rt_gk_sheath_2x2v_p1 -g -M -c 1 -d 2"
mpirun -np 2 ./rt_gk_sheath_2x2v_p1 -g -M -c 1 -d 2

exit 0
