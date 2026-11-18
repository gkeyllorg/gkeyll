# StellarAI CUDA/NCCL CI profile. Set GKEYLL_CI_ROOT and JAVA_HOME externally.
# Options, defaults and file types: ci/jenkins/slurm-clusters/README.md
# Setup and usage: ci/jenkins/README.slurm_cluster.md

SLURM_CI_CONFIG_VERSION=1
SLURM_CI_DISPLAY_NAME="StellarAI GPU"
SLURM_CI_ENV_SCRIPT=machines/module_load.stellarai-gpu.sh
SLURM_CI_MKDEPS_SCRIPT=machines/mkdeps.stellarai-gpu.sh
SLURM_CI_CONFIGURE_SCRIPT=machines/configure.stellarai-gpu.sh
SLURM_CI_PARTITION=gpu
SLURM_CI_QOS_REQUIRED=true
# Set SLURM_CI_QOS and, if required by your project, SLURM_CI_ACCOUNT externally.
# Match the site's mpirun example; the driver supplies each test's rank count.
SLURM_CI_MPIEXEC=mpirun
# Run one serial regression at a time on the allocated GPU.
SLURM_CI_SERIAL_JOBS=1
SLURM_CI_UNIT_SCRIPT=ci/jenkins/slurm-clusters/slurm-unit-tests.stellarai_gpu.sh
SLURM_CI_UNIT_CPUS_PER_TASK=4
SLURM_CI_UNIT_GRES=gpu:1
SLURM_CI_SERIAL_SCRIPT=ci/jenkins/slurm-clusters/slurm-serial-regression-tests.stellarai_gpu.sh
SLURM_CI_SERIAL_CPUS_PER_TASK=4
SLURM_CI_SERIAL_GRES=gpu:1
SLURM_CI_PARALLEL_SCRIPT=ci/jenkins/slurm-clusters/slurm-parallel-regression-tests.stellarai_gpu.sh
SLURM_CI_PARALLEL_GRES=gpu:4
SLURM_CI_PARALLEL_MEM=64G
# Defaults: one node; unit/serial one task; parallel four tasks, one CPU each;
# 00:30:00 unit and 04:00:00 regression wall time. Keep parallel on one node:
# the site scripts use UCX_TLS=sm,cuda_copy,self for intra-node communication.
