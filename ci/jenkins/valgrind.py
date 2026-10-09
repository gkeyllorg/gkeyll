#!/usr/bin/env python3
"""Prepare CPU builds and retain complete serial unit-test results for Jenkins."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import shlex
import subprocess
import sys


FLAGS = '-march=x86-64 -mtune=generic'
APPS = ('core', 'moments', 'vlasov', 'gyrokinetic', 'pkpm')


def enabled():
    value = os.environ.get('GKEYLL_USE_VALGRIND', '0')
    if value not in ('0', '1'):
        raise ValueError('GKEYLL_USE_VALGRIND must be 0 or 1')
    return value == '1'


def config_values():
    names = ('CC', 'USE_NCCL', 'USE_CUDSS', 'BUILD_DIR', 'BUILD_APP')
    makefile = 'include config.mak\nBUILD_DIR ?= build\n.PHONY: ci-values\nci-values:\n'
    makefile += ''.join("\t@printf '%s\\n' '{}=$({})'\n".format(n, n) for n in names)
    result = subprocess.run(['make', '-s', '--no-print-directory', '-f', '-', 'ci-values'],
                            input=makefile, text=True, capture_output=True, check=True)
    return dict(line.split('=', 1) for line in result.stdout.splitlines())


def configure():
    values = config_values()
    if (any(Path(word).name in ('nvcc', 'nvc', 'nvc++') for word in values['CC'].split())
            or values['USE_NCCL'] == '1' or values['USE_CUDSS'] == '1'):
        raise ValueError('Valgrind requires a CPU config.mak (no CUDA, NCCL, or cuDSS)')
    target = subprocess.check_output(shlex.split(values['CC']) + ['-dumpmachine'], text=True).strip()
    if target.startswith(('x86_64', 'amd64')):
        flags = FLAGS
    elif target.startswith(('aarch64', 'arm64')):
        flags = '-march=armv8-a -mtune=generic'
    else:
        raise ValueError(f'No Valgrind build profile for compiler target {target}')
    # Append after configuration, including for source revisions predating this
    # feature. Never modify the machine owner's prebuilt configuration.
    with Path('config.mak').open('a') as stream:
        stream.write('\n# Jenkins Valgrind CPU configuration\n'
                     'GKEYLL_USE_VALGRIND = 1\n'
                     f'override ARCH_FLAGS := {flags}\n'
                     f'CFLAGS := $(filter-out -m%,$(CFLAGS)) {flags} -g -gdwarf-4\n'
                     'SQL_CFLAGS ?= -fPIC -Wno-implicit-int-float-conversion\n'
                     f'SQL_CFLAGS := $(filter-out -m%,$(SQL_CFLAGS)) {flags} -g -gdwarf-4\n')


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    summary = {'status': 'running', 'tests': [], 'make_exit': None}
    summary_file = output / 'summary.json'

    def save():
        summary_file.write_text(json.dumps(summary, indent=2) + '\n')

    save()
    if not shutil.which('valgrind'):
        summary.update(status='failed', error='Valgrind is not installed on this CPU worker')
        save()
        return 1
    values = config_values()
    app = values['BUILD_APP'] or 'pkpm'
    if app == 'all':
        app = 'pkpm'
    if app not in APPS:
        raise ValueError(f'Unknown BUILD_APP: {app}')
    apps = APPS[:APPS.index(app) + 1]
    expected = [(layer, source.stem) for layer in apps
                for source in sorted(Path(layer, 'unit').glob('ctest_*.c'))]
    summary['layers'] = list(apps)
    try:
        # Use the same runner as core-valcheck, from the trusted CI checkout.
        # Historical candidate refs may have a valcheck recipe that always
        # succeeds and discards assertions, so do not depend on that recipe.
        summary['make_exit'] = subprocess.call(['make', app + '-unit'])
        if summary['make_exit'] == 0:
            runner = Path(__file__).with_name('valcheck.sh')
            if not runner.is_file():
                runner = Path(__file__).resolve().parents[2] / 'core/minus/valcheck.sh'
            for layer in apps:
                log_dir = output / layer
                log_dir.mkdir()
                result = subprocess.call([
                    'sh', str(runner),
                    *(str(Path(values['BUILD_DIR']) / layer / 'unit' / name)
                      for test_layer, name in expected if test_layer == layer)],
                    env=dict(os.environ, GKYL_VALGRIND_LOG_DIR=str(log_dir.resolve())))
                if result:
                    summary['make_exit'] = result
    finally:
        # Also retain partial logs when the command fails or is interrupted.
        for layer, name in expected:
            log = output / layer / (name + '.log')
            exit_file = output / layer / (name + '.log.exit')
            text = log.read_text(errors='replace') if log.exists() else ''
            counts = re.findall(r'ERROR SUMMARY:\s+([\d,]+) errors', text)
            errors = sum(int(n.replace(',', '')) for n in counts) if counts else None
            # Fail closed on missing/truncated output, assertions, and crashes.
            exit_code = exit_file.read_text().strip() if exit_file.exists() else None
            passed = errors == 0 and exit_code == '0'
            summary['tests'].append(dict(name=layer + '/' + name, errors=errors, exit=exit_code,
                                         status='passed' if passed else 'failed',
                                         log=layer + '/' + name + '.log' if log.exists() else None))
        summary['status'] = ('passed' if expected and summary['make_exit'] == 0 and
                             all(t['status'] == 'passed' for t in summary['tests']) else 'failed')
        save()
    return 0 if summary['status'] == 'passed' else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('configure', 'run'))
    parser.add_argument('--output', type=Path, default=Path('ci-valgrind/candidate'))
    args = parser.parse_args()
    if not enabled():
        return 0
    if args.action == 'configure':
        configure()
        return 0
    return run(args.output.resolve())


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        sys.exit(str(error))
