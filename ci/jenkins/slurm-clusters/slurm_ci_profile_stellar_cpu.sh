# Reviewed Stellar CPU CI profile. Operator environment overrides these defaults.
# Options, defaults and file types: ci/jenkins/slurm-clusters/README.md
# Setup and usage: ci/jenkins/README.slurm_cluster.md

SLURM_CI_CONFIG_VERSION=1
SLURM_CI_DISPLAY_NAME="Stellar CPU"
SLURM_CI_ROOT="/scratch/gpfs/${USER:?}/gkeyll_ci"
SLURM_CI_JAVA_HOME=/usr/lib/jvm/java-21-openjdk-21.0.13.0.9-1.1.el8.x86_64
SLURM_CI_ENV_SCRIPT=machines/module_load.stellar-intel.sh
SLURM_CI_MKDEPS_SCRIPT=machines/mkdeps.stellar-intel.sh
SLURM_CI_CONFIGURE_SCRIPT=machines/configure.stellar-intel.sh
SLURM_CI_QOS_REQUIRED=true
# Set SLURM_CI_QOS and, if required by your project, SLURM_CI_ACCOUNT externally.
SLURM_CI_UNIT_SCRIPT=ci/jenkins/slurm-clusters/slurm-unit-tests.stellar_cpu.sh
SLURM_CI_UNIT_CPUS_PER_TASK=4
SLURM_CI_SERIAL_SCRIPT=ci/jenkins/slurm-clusters/slurm-serial-regression-tests.stellar_cpu.sh
SLURM_CI_SERIAL_CPUS_PER_TASK=8
SLURM_CI_PARALLEL_SCRIPT=ci/jenkins/slurm-clusters/slurm-parallel-regression-tests.stellar_cpu.sh
# Defaults: one node; unit/serial one task; parallel four tasks, one CPU each;
# no GPU request; 00:30:00 unit and 04:00:00 regression wall time.
