#!/bin/bash -l

#.This jobscript is for full nodes (using all the GPUs/node).

#.Declare a name for this job, preferably 16 or fewer characters.
#SBATCH -J <Job Name>

#.Enter the account to charge.
#SBATCH -A <Account Number>

#.Specify a queue.
#SBATCH -q regular

#.Number of nodes to request (Perlmutter 4 GPUs per node).
#SBATCH -N 2

#.Number of MPI processes (and we match 1 process to 1 GPU).
#SBATCH --ntasks 8

#.Specify GPU needs:
#SBATCH --constraint gpu
#SBATCH --gpus 8

#.Request wall time
#SBATCH -t 00:30:00

#.Mail is sent to you when the job starts and when it terminates or aborts.
#SBATCH --mail-user=<your email>
#SBATCH --mail-type=END,FAIL,REQUEUE

#.Load modules (this must match those in the machines/configure script).
. ./machines/module_load.perlmutter-gpu.sh
module unload python

#.Run the rt_gk_sheath_2x2v_p1 executable using 1 GPU along x (-c 1) and 8
#.GPUs along the field line (-d 8). See './rt_gk_sheath_2x2v_p1 -h' for
#.more details/options on decomposition. It also assumes the executable is
#.in the present directory. If it isn't, change `./` to point to the
#.directory containing the executable.

echo "srun -u -n 8 --gpus 8 ./rt_gk_sheath_2x2v_p1 -g -M -c 1 -d 8"
srun -u -n 8 --gpus 8 ./rt_gk_sheath_2x2v_p1 -g -M -c 1 -d 8



