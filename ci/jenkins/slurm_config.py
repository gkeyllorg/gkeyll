#!/usr/bin/env python3
"""Load reviewed Bash profiles. Only validated data crosses into Jenkins or argv."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]*\Z")
LANES = ('UNIT', 'SERIAL', 'PARALLEL', 'BUILD')
DEFAULTS = dict(CONFIG_VERSION='1', DISPLAY_NAME='', JOB_NAME='', CACHE_NAMESPACE='',
    STATUS_CONTEXT='', ROOT='', JAVA_HOME='', PORT='8080', SESSION='gkeyll_ci',
    ENV_SCRIPT='', MKDEPS_SCRIPT='', CONFIGURE_SCRIPT='', BUILD_LOCATION='agent',
    BUILD_JOBS='3', SERIAL_JOBS='4', SERIAL_TEST_TIMEOUT='900',
    ACCOUNT='', QOS='', PARTITION='', CONSTRAINT='', ACCOUNT_REQUIRED='false',
    QOS_REQUIRED='false', MPIEXEC='srun', MPI_ARGS='')
for lane in LANES:
    DEFAULTS.update({f'{lane}_{key}': value for key, value in dict(SCRIPT='', NODES='1',
        TASKS='4' if lane == 'PARALLEL' else '1', CPUS_PER_TASK='1', GPUS_PER_TASK='0',
        TIME='00:30:00' if lane == 'UNIT' else '04:00:00', MEM='', GRES='',
        ACCOUNT='', QOS='', PARTITION='', CONSTRAINT='').items()})
ARRAYS = ['SBATCH_EXTRA_ARGS'] + [f'{lane}_EXTRA_ARGS' for lane in LANES]
PREFIX = 'SLURM_CI_'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def load(profile, root=HERE, environ=None):
    environ = dict(os.environ if environ is None else environ)
    require(ID.fullmatch(profile) and profile not in ('personal', 'team', 'slurm'), 'invalid Slurm profile ID')
    path = root / 'slurm-clusters' / f'slurm_ci_profile_{profile}.sh'
    require(path.is_file(), f'profile does not exist: {path}')
    defaults = dict(DEFAULTS, DISPLAY_NAME=profile, JOB_NAME=f'gkeyll-ci-{profile}',
                    CACHE_NAMESPACE=profile, STATUS_CONTEXT=f'continuous-integration/jenkins/{profile}')
    # Generate references from our fixed variable vocabulary, never from profile text.
    script = ['set -eo pipefail']
    script += [f'{PREFIX}{key}={shlex.quote(value)}' for key, value in defaults.items()]
    script += [f'{PREFIX}{key}=()' for key in ARRAYS]
    script += ['source "$1" >&2', 'compgen -A variable SLURM_CI_ | while IFS= read -r key; do printf "%s\\0" "$key"; done', 'printf "\\0"']
    for key in defaults:
        script.append(f'printf "%s\\0" "${{{PREFIX}{key}}}"')
    for key in ARRAYS:
        script += [f'printf "%s\\0" "${{#{PREFIX}{key}[@]}}"',
                   f'if (( ${{#{PREFIX}{key}[@]}} )); then printf "%s\\0" "${{{PREFIX}{key}[@]}}"; fi']
    shell_env = {k: v for k, v in environ.items() if not k.startswith(PREFIX)}
    output = subprocess.check_output(['bash', '-c', '\n'.join(script), 'profile', str(path)],
                                     env=shell_env, cwd=root.parent.parent)
    fields = iter(output.decode().split('\0'))
    known = {PREFIX + k for k in [*defaults, *ARRAYS]}
    for key in fields:
        if not key:
            break
        require(key in known, f'unknown profile variable: {key}')
    config = {key: next(fields) for key in defaults}
    config.update({key: [next(fields) for _ in range(int(next(fields)))] for key in ARRAYS})
    for key in defaults:
        if PREFIX + key in environ:
            config[key] = environ[PREFIX + key]
    for key, name in dict(ROOT='GKEYLL_CI_ROOT', JAVA_HOME='JAVA_HOME', PORT='JENKINS_PORT',
                          SESSION='JENKINS_SESSION', JOB_NAME='JENKINS_JOB').items():
        if name in environ:
            config[key] = environ[name]
    require(config['CONFIG_VERSION'] == '1', 'unsupported profile version')
    require(ID.fullmatch(config['CACHE_NAMESPACE']), 'invalid cache namespace')
    for key in ('JOB_NAME', 'SESSION'):
        require(ID.fullmatch(config[key]), f'{key} must be a simple name')
    require(re.fullmatch(r'[A-Za-z0-9_./-]+', config['STATUS_CONTEXT']), 'invalid status context')
    require(config['BUILD_LOCATION'] in ('agent', 'slurm'), 'BUILD_LOCATION must be agent or slurm')
    for key in ('ACCOUNT_REQUIRED', 'QOS_REQUIRED'):
        require(config[key] in ('true', 'false'), f'{key} must be true or false')
    for key, value in config.items():
        if isinstance(value, str):
            require(not any(c in value for c in '\n\r\0'), f'{key} contains a control character')
    for key in ('BUILD_JOBS', 'SERIAL_JOBS', 'SERIAL_TEST_TIMEOUT', 'PORT'):
        require(re.fullmatch(r'[1-9][0-9]*', config[key]), f'{key} must be positive')
    require(int(config['PORT']) <= 65535, 'invalid port')
    for key in ('ROOT', 'JAVA_HOME'):
        value = config[key]
        require(not value or (value.startswith('/') and value != '/' and '..' not in Path(value).parts
                and not any(c in value for c in "'\"`$\\")), f'{key} must be a safe absolute directory')
    script_keys = ['ENV_SCRIPT', 'MKDEPS_SCRIPT', 'CONFIGURE_SCRIPT'] + [f'{lane}_SCRIPT' for lane in LANES]
    for key in script_keys:
        value = config[key]
        if key == 'BUILD_SCRIPT' and config['BUILD_LOCATION'] == 'agent' and not value:
            continue
        require(re.fullmatch(r'[A-Za-z0-9_./-]+', value) and not value.startswith('/')
                and '..' not in Path(value).parts, f'{key} must be a repository-relative script')
        script_path = root.parent.parent / value
        require(script_path.is_file(), f'missing {key}: {script_path}')
        subprocess.run(['bash', '-n', str(script_path)], check=True)
    for lane in LANES:
        for key in ('NODES', 'TASKS', 'CPUS_PER_TASK'):
            require(re.fullmatch(r'[1-9][0-9]*', config[f'{lane}_{key}']), f'{lane}_{key} must be positive')
        require(re.fullmatch(r'0|[1-9][0-9]*', config[f'{lane}_GPUS_PER_TASK']), f'{lane}: invalid GPU count')
        require(not (config[f'{lane}_GRES'] and int(config[f'{lane}_GPUS_PER_TASK'])), f'{lane}: GRES conflicts with GPUS_PER_TASK')
        require(re.fullmatch(r'(?:[0-9]+-)?[0-9]+:[0-5][0-9]:[0-5][0-9]', config[f'{lane}_TIME']), f'{lane}: invalid wall time')
    # Only extra resource options with unambiguous long names and explicit values.
    allowed = {'--reservation', '--licenses', '--hint', '--cpu-bind', '--gpu-bind',
               '--gpus-per-node', '--ntasks-per-node', '--exclusive', '--mem-per-cpu', '--gres-flags'}
    for key in ARRAYS:
        for arg in config[key]:
            option = arg.split('=', 1)[0]
            require(option in allowed and (option == '--exclusive' or '=' in arg)
                    and not any(c in arg for c in '\n\r\0'), f'{key}: unsupported extra resource argument {arg!r}')
    for lane in LANES:
        options = [arg.split('=', 1)[0] for arg in config['SBATCH_EXTRA_ARGS'] + config[f'{lane}_EXTRA_ARGS']]
        require(len(options) == len(set(options)), f'{lane}: duplicate extra resource option')
        require(not ('--gpus-per-node' in options and
                    (config[f'{lane}_GRES'] or int(config[f'{lane}_GPUS_PER_TASK']))),
                f'{lane}: gpus-per-node conflicts with GRES or GPUS_PER_TASK')
        require(not ('--mem-per-cpu' in options and config[f'{lane}_MEM']),
                f'{lane}: mem-per-cpu conflicts with MEM')
    config['PROFILE'] = profile
    return config


def allocation(config, lane):
    lane = lane.upper()
    require(lane in LANES, 'invalid lane')
    args = [f'--nodes={config[lane+"_NODES"]}', f'--ntasks={config[lane+"_TASKS"]}',
            f'--cpus-per-task={config[lane+"_CPUS_PER_TASK"]}', f'--time={config[lane+"_TIME"]}']
    for key in ('ACCOUNT', 'QOS', 'PARTITION', 'CONSTRAINT', 'MEM', 'GRES'):
        value = config.get(f'{lane}_{key}') or config.get(key, '')
        if value:
            args.append(f'--{key.lower()}={value}')
    if int(config[f'{lane}_GPUS_PER_TASK']):
        args.append(f'--gpus-per-task={config[lane+"_GPUS_PER_TASK"]}')
    return args + config['SBATCH_EXTRA_ARGS'] + config[f'{lane}_EXTRA_ARGS']


def validate_deployment(config):
    for key in ('ACCOUNT', 'QOS'):
        if config[key + '_REQUIRED'] == 'true':
            for lane in LANES:
                if lane == 'BUILD' and config['BUILD_LOCATION'] == 'agent':
                    continue
                require(config[key] or config[lane + '_' + key], f'{key} is required for {lane}')


def fingerprint(config, revision):
    excluded = {'DISPLAY_NAME', 'JOB_NAME', 'STATUS_CONTEXT', 'ROOT', 'JAVA_HOME', 'PORT', 'SESSION'}
    data = {k: v for k, v in config.items() if k not in excluded}
    return hashlib.sha256(json.dumps([revision, data], sort_keys=True).encode()).hexdigest()


def main():
    # Client flags belong to the client, including --help.
    if len(sys.argv) >= 3 and sys.argv[1] == 'client':
        config = load(sys.argv[2])
        env = dict(os.environ, GKEYLL_CI_CLUSTER=sys.argv[2])
        for name, key in dict(GKEYLL_CI_ROOT='ROOT', JAVA_HOME='JAVA_HOME', JENKINS_JOB='JOB_NAME',
                              JENKINS_PORT='PORT', JENKINS_SESSION='SESSION').items():
            env.setdefault(name, config[key])
        os.execve(str(HERE / 'jenkins-slurm_cluster.sh'),
                  ['jenkins-slurm_cluster.sh', '__configured', *sys.argv[3:]], env)
    if len(sys.argv) == 4 and sys.argv[1] == 'allocation':
        config = json.loads(Path(sys.argv[2]).read_text())
        print('\n'.join(allocation(config, sys.argv[3])))
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('resolve', 'validate', 'client'))
    parser.add_argument('profile')
    parser.add_argument('--revision', default='local-validation')
    parser.add_argument('--output')
    parser.add_argument('--env-output')
    args, client_args = parser.parse_known_args()
    require(args.command == 'client' or not client_args, f'unexpected arguments: {client_args}')
    config = load(args.profile)
    validate_deployment(config)
    config['CONFIG_HASH'] = fingerprint(config, args.revision)
    config['ALLOCATIONS'] = {lane: allocation(config, lane) for lane in LANES}
    if args.env_output:
        Path(args.env_output).write_text(''.join(f'{key}={value}\n' for key, value in config.items() if isinstance(value, str)))
    output = json.dumps(config, indent=2) + '\n'
    if args.output:
        Path(args.output).write_text(output)
    else:
        print(output, end='')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        sys.exit(f'slurm-config: {error}')
