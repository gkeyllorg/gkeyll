#!/usr/bin/env python3
"""Copy configured dependencies into one CI run and write its build config."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile

GROUPS = (
    ('', ('LAPACK_INC_DIR', 'LAPACK_LIB_DIR')),
    ('', ('SUPERLU_INC_DIR', 'SUPERLU_LIB_DIR')),
    ('USE_MPI', ('CONF_MPI_INC_DIR', 'CONF_MPI_LIB_DIR')),
    ('USE_LUA', ('CONF_LUA_INC_DIR', 'CONF_LUA_LIB_DIR')),
    ('USE_NCCL', ('CONF_NCCL_INC_DIR', 'CONF_NCCL_LIB_DIR')),
    ('USE_CUDSS', ('CONF_CUDSS_INC_DIR', 'CONF_CUDSS_LIB_DIR')),
)
VARIABLES = ('PREFIX', 'BUILD_APP', 'CC', 'USE_LAPACK_LITE', 'CUDAMATH_LIB_DIR') + tuple(dict.fromkeys(
    name for flag, paths in GROUPS for name in (flag, *paths) if name))


def evaluate_config(source):
    """Let Make evaluate expressions; do not approximate its assignment syntax."""
    with tempfile.TemporaryDirectory(prefix='gkeyll-config-values-') as directory:
        root = Path(directory)
        makefile = root / 'probe.mk'
        escaped = str(source).replace(' ', '\\ ').replace('#', '\\#')
        makefile.write_text('include ' + escaped + '\n.PHONY: ci-values\nci-values:\n'
                           + ''.join('\t$(file >' + str(root / name) + ',$(' + name + '))\n'
                                     for name in VARIABLES) + '\t@:\n')
        environment = dict(os.environ, MAKEFLAGS='')
        environment.pop('PREFIX', None)
        environment.pop('INSTALL_PREFIX', None)
        subprocess.run(['make', '-s', '-f', str(makefile), 'ci-values'],
                       check=True, env=environment, stdout=subprocess.PIPE)
        return {name: (root / name).read_text().strip() for name in VARIABLES}


def system_path(path):
    # OS libraries and platform SDKs remain part of the machine's toolchain.
    return any(path == root or root in path.parents for root in
               map(Path, ('/usr', '/lib', '/lib64', '/System')))


def dependency_root(path):
    for index in range(len(path.parts) - 1, 0, -1):
        if path.parts[index] in ('include', 'lib', 'lib64', 'bin'):
            return Path(*path.parts[:index])
    return path


def remap(path, mappings):
    for old, new in sorted(mappings.items(), key=lambda pair: len(pair[0]), reverse=True):
        if path == old or path.startswith(old + '/'):
            return new + path[len(old):]
    return path


def prepare_dependencies(source, destination, mpiexec):
    source_text = source.read_text()
    if not re.search(r'(?m)^(?:override\s+)?PREFIX\s*[:?]?=', source_text):
        raise ValueError('Prebuilt config.mak must define PREFIX')
    digest = hashlib.sha256(source_text.encode()).hexdigest()
    manifest_file = destination / 'manifest.json'
    if manifest_file.is_file():
        manifest = json.loads(manifest_file.read_text())
        if manifest['source_sha256'] != digest:
            raise ValueError('Prebuilt config changed during this run; start a new CI run')
        if 'adas_dir' not in manifest:
            raise ValueError('Dependency snapshot predates ADAS copying; start a new CI run')
        return manifest
    if destination.exists():
        raise ValueError(f'Incomplete dependency copy already exists: {destination}')
    values = evaluate_config(source)
    original_prefix = Path(values['PREFIX'])
    if original_prefix == Path('/') or not original_prefix.is_absolute() or not original_prefix.is_dir():
        raise ValueError(f'Prebuilt PREFIX is not an absolute directory: {original_prefix}')
    adas_source = original_prefix / 'gkeyll/share/adas'
    if values['BUILD_APP'] in ('gyrokinetic', 'pkpm') and not any(
            path.is_file() for path in adas_source.glob('*.npy')):
        raise ValueError(f'Missing ADAS .npy data in {adas_source}; '
                         'prepare the prebuilt dependencies with --build-adas=yes')
    paths = {}
    for flag, names in GROUPS:
        if flag and values[flag] != '1':
            continue
        if names[0] == 'LAPACK_INC_DIR' and values['USE_LAPACK_LITE'] == '1':
            continue
        paths.update((name, values[name]) for name in names if values[name])
    if values['CC'] == 'nvcc' and values['CUDAMATH_LIB_DIR']:
        paths['CUDAMATH_LIB_DIR'] = values['CUDAMATH_LIB_DIR']
    roots = set()
    for name, value in paths.items():
        path = Path(value)
        if not path.is_absolute() or not path.is_dir():
            raise ValueError(f'{name} is not an absolute dependency directory: {value}')
        if not system_path(path):
            roots.add(dependency_root(path))
    if mpiexec and Path(mpiexec).is_absolute() and not system_path(Path(mpiexec)):
        if not Path(mpiexec).is_file():
            raise ValueError(f'MPI launcher does not exist: {mpiexec}')
        roots.add(dependency_root(Path(mpiexec)))

    mappings = {} if system_path(original_prefix) else {
        str(original_prefix): str(destination / 'gkylsoft')}
    copied = {}
    canonical_copies = {}
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='.dependencies-', dir=destination.parent))
    try:
        for root in sorted(roots, key=lambda path: (len(path.parts), str(path))):
            if any(parent == root or parent in root.parents for parent in copied):
                continue
            canonical = root.resolve()
            if canonical in canonical_copies:
                mappings[str(root)] = str(canonical_copies[canonical])
                continue
            if destination == canonical or canonical in destination.parents:
                raise ValueError(f'Dependency destination is inside its source: {root}')
            if original_prefix == root or original_prefix in root.parents:
                relative = Path('gkylsoft') / root.relative_to(original_prefix)
            else:
                key = hashlib.sha256(str(root).encode()).hexdigest()[:12]
                relative = Path('external') / key / root.name
            target = destination / relative
            # Dereference links so neither versioned install aliases nor library
            # symlinks retain references to the machine owner's installation.
            def ignore(directory, names):
                if Path(directory) == original_prefix:
                    return {'gkeyll', 'gkeyll-results'} & set(names)
                return set()
            shutil.copytree(root, stage / relative, symlinks=False, ignore=ignore)
            copied[root] = target
            canonical_copies[canonical] = target
            mappings[str(root)] = str(target)
            mappings[str(canonical)] = str(target)
        # ADAS is a runtime dependency inside the otherwise excluded Gkeyll
        # installation. Snapshot it separately, including for system prefixes.
        adas_dir = ''
        if adas_source.is_dir():
            relative = Path('data/adas')
            shutil.copytree(adas_source, stage / relative, symlinks=False)
            adas_dir = str(destination / relative)
        paths = {name: remap(value.rstrip('/'), mappings) for name, value in paths.items()}
        libdirs = list(dict.fromkeys(value for name, value in paths.items() if 'LIB_DIR' in name))
        mpi_home = ''
        if 'CONF_MPI_LIB_DIR' in paths:
            mpi_home = str(dependency_root(Path(paths['CONF_MPI_LIB_DIR'])))
            if not mpiexec and (Path(mpi_home) / 'bin/mpiexec').is_file():
                mpiexec = str(Path(mpi_home) / 'bin/mpiexec')
            elif not mpiexec:
                mpiexec = str(dependency_root(Path(values['CONF_MPI_LIB_DIR'])) / 'bin/mpiexec')
        manifest = dict(source_config=str(source), source_sha256=digest,
                        mappings=mappings, paths=paths, library_dirs=libdirs,
                        mpi_home=mpi_home, mpiexec=remap(mpiexec, mappings),
                        adas_dir=adas_dir)
        (stage / 'source-config.mak').write_text(source_text)
        (stage / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
        environment = runtime_environment(manifest, destination)
        exports = ''.join(
            'export ' + key + '=' + shlex.quote(value) + '\n'
            for key, value in environment.items() if value and key != 'LD_LIBRARY_PATH')
        exports += ('export LD_LIBRARY_PATH=' + shlex.quote(':'.join(libdirs))
                    + '${LD_LIBRARY_PATH:+:"$LD_LIBRARY_PATH"}\n')
        (stage / 'env.sh').write_text(exports)
        stage.rename(destination)
        return manifest
    finally:
        if stage.exists():
            shutil.rmtree(stage)


def runtime_environment(manifest, destination):
    previous = (remap(path, manifest['mappings'])
                for path in os.environ.get('LD_LIBRARY_PATH', '').split(':') if path)
    libdirs = list(dict.fromkeys([*manifest['library_dirs'], *previous]))
    return dict(LD_LIBRARY_PATH=':'.join(libdirs), MPI_HOME=manifest['mpi_home'],
                OPAL_PREFIX=manifest['mpi_home'], MPIEXEC=manifest['mpiexec'],
                GKEYLL_CI_DEPENDENCY_ENV=str(destination / 'env.sh'))


def write_config(source, prefix, destination, output, mpiexec=''):
    manifest = prepare_dependencies(source, destination, mpiexec)
    if manifest['adas_dir']:
        # Each build owns its share directory: Make installs revision-specific
        # radiation fits here, so do not link it to the shared dependency copy.
        shutil.copytree(manifest['adas_dir'], prefix / 'gkeyll/share/adas',
                        symlinks=False, dirs_exist_ok=True)
    config = (destination / 'source-config.mak').read_text()
    config = re.sub(r'(?m)^(override\s+)?PREFIX(?=\s*[:?]?=)',
                    lambda match: (match[1] or '') + 'GKEYLL_CI_DEPENDENCY_PREFIX', config)
    config = config.replace('$(PREFIX)', '$(GKEYLL_CI_DEPENDENCY_PREFIX)')
    config = config.replace('${PREFIX}', '${GKEYLL_CI_DEPENDENCY_PREFIX}')
    alternatives = '|'.join(re.escape(old) for old in
                            sorted(manifest['mappings'], key=len, reverse=True))
    if alternatives:
        config = re.sub('(' + alternatives + r')(?=/|[\s\'"#;)]|$)',
                        lambda match: manifest['mappings'][match[1]], config)
    config += '\n# Private dependency snapshot for this CI run\n'
    for name, value in manifest['paths'].items():
        config += f'override {name} := {value}\n'
    # Keep the literal assignment for historical regression-tool readers.
    config += (f'\n# Jenkins run installation\nPREFIX={prefix}\n'
               f'override PREFIX := {prefix}\noverride INSTALL_PREFIX := {prefix}\n')
    # Legacy configs can pin GKYL_SHARE_DIR to the original installation.
    # Generic dependency remapping points that at an excluded gkeyll tree;
    # compile against this build's data copy and revision-specific fits instead.
    config += f'override GKYL_SHARE_DIR := {prefix}/gkeyll/share\n'
    output.write_text(config)
    return runtime_environment(manifest, destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--prefix', type=Path, required=True)
    parser.add_argument('--dependencies', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path('config.mak'))
    parser.add_argument('--mpiexec', default='')
    parser.add_argument('--format', choices=('json', 'env'), default='json',
                        help='Runtime environment output (env uses literal KEY=value lines)')
    args = parser.parse_args()
    for path in (args.config, args.prefix, args.dependencies):
        if not path.is_absolute():
            parser.error(f'Expected an absolute path: {path}')
    try:
        result = write_config(args.config, args.prefix, args.dependencies, args.output, args.mpiexec)
        if args.format == 'env':
            # Jenkins can read these using sandbox-approved string operations.
            # Reject line breaks so one value cannot turn into multiple entries.
            if any('\n' in value or '\r' in value for value in result.values()):
                raise ValueError('Runtime environment values must not contain line breaks')
            output = '\n'.join(f'{key}={value}' for key, value in result.items())
        else:
            output = json.dumps(result)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f'prebuilt-config: {error}\n')
    print(output)


if __name__ == '__main__':
    main()
