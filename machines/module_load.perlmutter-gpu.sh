#!/bin/bash
# Shared Perlmutter GPU build and runtime environment.
module load PrgEnv-gnu/8.6.0
module load craype-accel-nvidia80
module load cray-mpich/9.0.1
module load cudatoolkit/13.0
module load nccl/2.29.2-cu13
module load cray-libsci/25.09.0
module load python

# Gkeyll uses NCCL for GPU communication rather than CUDA-aware MPI.
export MPICH_GPU_SUPPORT_ENABLED=0
export DVS_MAXNODES=24
export MPICH_MPIIO_DVS_MAXNODES=24
export NCCL_NET_GDR_LEVEL=PHB
export NCCL_NET='AWS Libfabric'
export NCCL_CROSS_NIC=1
export FI_CXI_DISABLE_HOST_REGISTER=1
export FI_CXI_RDZV_GET_MIN=0
export FI_CXI_RDZV_THRESHOLD=0
export FI_CXI_RDZV_EAGER_SIZE=0
