"""Exercise the runner's comparisons with real G0.Zero readers and tiny files.

Set GKEYLL to a built executable and make its shared libraries discoverable.
"""
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
GKEYLL = os.environ.get('GKEYLL') or shutil.which('gkeyll')
SOURCE = (ROOT / 'gkeyll/lua/Tool/runregression.lua').read_text()
COMPARISONS = SOURCE[SOURCE.index('local function shortPath('):
                     SOURCE.index('-- ---- GPU comparison helpers')]


def pack(value):
    """Encode the small subset of MessagePack used by these fixtures."""
    if isinstance(value, str):
        data = value.encode()
        return b'\xd9' + bytes([len(data)]) + data
    if isinstance(value, int):
        return b'\xd3' + struct.pack('>q', value)
    if isinstance(value, float):
        return b'\xcb' + struct.pack('>d', value)
    if isinstance(value, list):
        return bytes([0x90 + len(value)]) + b''.join(map(pack, value))
    return bytes([0x80 + len(value)]) + b''.join(
        pack(key) + pack(item) for key, item in value.items())


def meta_file(meta, file_type=5):
    payload = pack(meta)
    return b'gkyl0' + struct.pack('=QQQ', 1, file_type, len(payload)) + payload


def field_file(value):
    # Type 1: double, one dimension, one cell, [0, 1], one component.
    return (b'gkyl0' + struct.pack('=QQQQQQddQQd',
            1, 1, 0, 2, 1, 1, 0., 1., 8, 1, value))


def parallel_field_file(tot_cells=6):
    # Two uneven MPI blocks, with a 2D grid and two values per cell.
    header = (b'gkyl0' + struct.pack('=QQQ', 1, 3, 0)
              + struct.pack('=QQQQddddQQQ', 2, 2, 3, 2,
                            0., 0., 1., 1., 16, tot_cells, 2))
    first = struct.pack('=QQQQQ8d', 1, 1, 2, 2, 4, *range(8))
    second = struct.pack('=QQQQQ4d', 3, 1, 3, 2, 2, *range(8, 12))
    return header + first + second


@unittest.skipUnless(GKEYLL, 'Set GKEYLL to run the comparison fixtures')
class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix='gkeyll-comparison-')
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.accepted = self.root / 'moments/creg-accepted/cpu_serial/rt_fixture'
        self.run = self.root / 'run'
        self.accepted.mkdir(parents=True)
        self.run.mkdir()
        self.meta = dict(time=0., frame=0, topo_file='rt_fixture_btopo.gkyl',
                         app_name='rt_fixture')
        self.files = {
            'rt_fixture-euler_0.gkyl': meta_file(self.meta),
            'rt_fixture_b0-euler_0.gkyl': field_file(1.),
            'rt_fixture_btopo.gkyl': meta_file(
                dict(ndim=1, num_blocks=1, connections=[0, 0, 5, 0, 0, 5]), 4),
        }
        for name, data in self.files.items():
            (self.accepted / name).write_bytes(data)
            (self.run / name).write_bytes(data)

    def check(self, passed=True, detail=''):
        script = self.root / 'compare.lua'
        script.write_text('''
local luaRoot = %s
package.path = luaRoot .. "/?.lua;" .. luaRoot .. "/Lib/?.lua;" .. package.path
local lfs = require "lfs"
local configVals = {results_dir=%s}
local runMode = "cpu_serial"
local layerCounts = {moments={passed=0, failed=0}}
local function log(msg) end
local function verboseLog(msg) end
local function basename(path) return path:match("([^/]+)$") end
local function stripext(path) return path:gsub("%%.[^.]+$", "") end
''' % (json.dumps(str(ROOT / 'gkeyll/lua')), json.dumps(str(self.root)))
            + COMPARISONS + '''
local status, detail = check_action(
   {name="moments/creg/rt_fixture", layer="moments", src="rt_fixture.c"},
   %s, "c")
assert(status == %d, detail)
assert(detail:find(%s, 1, true), detail)
''' % (json.dumps(str(self.run)), 1 if passed else 0, json.dumps(detail)))
        result = subprocess.run([GKEYLL, '-S', str(script)], cwd=self.root,
                                capture_output=True, text=True, timeout=20)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn('*** ERROR:', result.stderr)

    def test_identical_outputs(self):
        self.check()

    def test_metadata_map_order(self):
        (self.run / 'rt_fixture-euler_0.gkyl').write_bytes(
            meta_file(dict(reversed(list(self.meta.items())))))
        self.check()

    def test_changed_metadata(self):
        for key, value in dict(time=1., frame=1, topo_file='other.gkyl',
                               app_name='other', extra=1).items():
            with self.subTest(key=key):
                (self.run / 'rt_fixture-euler_0.gkyl').write_bytes(
                    meta_file(dict(self.meta, **{key: value})))
                self.check(False, 'multiblock metadata mismatch: ' + key)

    def test_malformed_metadata(self):
        valid = self.files['rt_fixture-euler_0.gkyl']
        for data in (valid[:-1], valid + b'x', meta_file({}),
                     meta_file(dict(self.meta, time=float('nan'))),
                     b'gkyl0' + struct.pack('=QQQ', 1, 5, 1) + b'\xc1'):
            with self.subTest(data=data):
                # Identical corruption must fail too.
                for directory in (self.accepted, self.run):
                    (directory / 'rt_fixture-euler_0.gkyl').write_bytes(data)
                self.check(False, 'multiblock metadata read failed')

    def test_changed_block_data(self):
        for value in (0., 2.):
            with self.subTest(value=value):
                (self.run / 'rt_fixture_b0-euler_0.gkyl').write_bytes(field_file(value))
                self.check(False, 'rt_fixture_b0-euler_0.gkyl  [DIFF]  max_abs=1')

    def test_block_data_within_tolerance(self):
        (self.run / 'rt_fixture_b0-euler_0.gkyl').write_bytes(field_file(1. + 1e-13))
        self.check()

    def test_parallel_arrays(self):
        for directory in (self.accepted, self.run):
            (directory / 'rt_fixture_b0-euler_0.gkyl').write_bytes(parallel_field_file())
        self.check()

    def test_parallel_array_read_failures_identify_side(self):
        name = 'rt_fixture_b0-euler_0.gkyl'
        valid = parallel_field_file()
        # The historical BGK writer advertised only configuration-space cells.
        for corrupt in (parallel_field_file(tot_cells=3), valid[:-8]):
            for bad_baseline, bad_candidate, detail in (
                    (True, False, 'baseline array read failed (candidate readable)'),
                    (False, True, 'candidate array read failed'),
                    (True, True, 'candidate array read failed')):
                with self.subTest(baseline=bad_baseline, candidate=bad_candidate):
                    (self.accepted / name).write_bytes(corrupt if bad_baseline else valid)
                    (self.run / name).write_bytes(corrupt if bad_candidate else valid)
                    self.check(False, detail)

    def test_changed_topology(self):
        (self.run / 'rt_fixture_btopo.gkyl').write_bytes(meta_file(
            dict(ndim=1, num_blocks=1, connections=[0, 0, 1, 0, 0, 5]), 4))
        self.check(False, 'topology mismatch')

    def test_missing_outputs(self):
        for name, data in self.files.items():
            with self.subTest(name=name):
                (self.run / name).unlink()
                self.check(False, name + '  [MISSING]')
                (self.run / name).write_bytes(data)


if __name__ == '__main__':
    unittest.main()
