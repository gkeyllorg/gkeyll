# Reviewed Perlmutter CUDA/NCCL CI profile. Set GKEYLL_CI_ROOT and JAVA_HOME externally.
# Options, defaults and file types: ci/jenkins/slurm-clusters/README.md
# Setup and usage: ci/jenkins/README.slurm_cluster.md

SLURM_CI_CONFIG_VERSION=1
SLURM_CI_DISPLAY_NAME="Perlmutter GPU"
SLURM_CI_ENV_SCRIPT=machines/module_load.perlmutter-gpu.sh
SLURM_CI_MKDEPS_SCRIPT=machines/mkdeps.perlmutter.gpu.sh
SLURM_CI_CONFIGURE_SCRIPT=machines/configure.perlmutter.gpu.sh
SLURM_CI_ACCOUNT_REQUIRED=true
SLURM_CI_QOS=shared
SLURM_CI_CONSTRAINT=gpu
SLURM_CI_UNIT_SCRIPT=ci/jenkins/slurm-clusters/slurm-unit-tests.perlmutter_gpu.sh
SLURM_CI_UNIT_CPUS_PER_TASK=32
SLURM_CI_UNIT_GPUS_PER_TASK=1
SLURM_CI_SERIAL_SCRIPT=ci/jenkins/slurm-clusters/slurm-serial-regression-tests.perlmutter_gpu.sh
SLURM_CI_SERIAL_CPUS_PER_TASK=32
SLURM_CI_SERIAL_GPUS_PER_TASK=1
SLURM_CI_PARALLEL_SCRIPT=ci/jenkins/slurm-clusters/slurm-parallel-regression-tests.perlmutter_gpu.sh
SLURM_CI_PARALLEL_CPUS_PER_TASK=8
SLURM_CI_PARALLEL_GPUS_PER_TASK=1
# Defaults: one node; unit/serial one task; parallel four tasks;
# 00:30:00 unit and 04:00:00 regression wall time. Scripts own srun/binding.
