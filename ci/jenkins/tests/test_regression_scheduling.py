"""Exercise the runner's Lua helpers with tiny jobs, without simulations.

Set LUAJIT to the LuaJIT executable when it is not on PATH.
"""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
LUAJIT = os.environ.get('LUAJIT') or shutil.which('luajit')
SOURCE = (ROOT / 'gkeyll/lua/Tool/runregression.lua').read_text()


def fragment(start, end):
    return SOURCE[SOURCE.index(start):SOURCE.index(end, SOURCE.index(start))]


@unittest.skipUnless(LUAJIT, 'Set LUAJIT to run the scheduling fixtures')
class SchedulingTests(unittest.TestCase):
    def test_mpi_worker_budget(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lua = r'''
local root = arg[1]
local function hasCmd(name) return os.execute("command -v " .. name .. " >/dev/null") == 0 end
local function log(msg) end
'''
            lua += fragment('local function shellQuote(', 'local function hasGkylOutput(')
            lua += fragment('local function wrapWithTimeout(', '-- ---- Test classification predicates')
            lua += 'local executeBatch\n'
            lua += fragment('executeBatch = function(', 'local function runLuaTest(')
            lua += r'''
local function run(ranks, budget, durations, commands)
   local items = {}
   for idx, count in ipairs(ranks) do
      local dir = root .. "/job" .. idx
      assert(os.execute("mkdir -p " .. shellQuote(dir)) == 0)
      items[idx] = {test={name="job" .. idx}, runDir=dir, parallelRanks=count,
         cmd=wrapWithTimeout(commands and commands[idx] or "sleep " .. durations[idx], 1, dir)}
   end
   local completions = {}
   local results = executeBatch(items, budget, function(item, result, completed, total)
      assert(completed == #completions + 1 and total == #items)
      table.insert(completions, item.test.name)
   end, true)
   local starts, ends, events = {}, {}, {}
   for idx, item in ipairs(items) do
      local file = assert(io.open(item.runDir .. "/_parallel_out.txt"))
      local raw = file:read("*a"); file:close()
      starts[idx] = assert(tonumber(raw:match("__START__:([%d.]+)")))
      ends[idx] = assert(tonumber(raw:match("__END__:([%d.]+)")))
      local slots = math.min(item.parallelRanks, budget)
      table.insert(events, {time=starts[idx], slots=slots})
      table.insert(events, {time=ends[idx], slots=-slots})
   end
   table.sort(events, function(a, b)
      return a.time < b.time or (a.time == b.time and a.slots < b.slots)
   end)
   local active = 0
   for _, event in ipairs(events) do
      active = active + event.slots
      assert(active <= budget, "MPI rank budget exceeded")
   end
   assert(active == 0)
   return starts, ends, completions, results
end
-- Four-rank tests overlap and refill without waiting for the slow first test.
local starts, ends, completions = run({4, 4, 4, 4}, 8, {0.8, 0.1, 0.1, 0.1})
assert(completions[1] == "job2")
assert(starts[2] < ends[1] and starts[3] < ends[1] and starts[4] < ends[1])
assert(starts[3] >= ends[2] and starts[4] >= ends[3])
-- Different decompositions consume their actual rank counts.
starts, ends = run({4, 2, 2, 4}, 8, {0.8, 0.1, 0.1, 0.1})
assert(starts[3] < ends[2], "smaller decompositions did not share capacity")
assert(starts[4] < ends[1] and starts[4] >= math.max(ends[2], ends[3]))
-- The default budget, and budgets below one test's size, run it alone.
for _, budget in ipairs({1, 3, 4, 7}) do
   starts, ends = run({4, 4}, budget, {0.1, 0.1})
   assert(starts[2] >= ends[1])
end
-- Failed and timed-out collectives release every reserved slot.
local results
starts, ends, completions, results = run({4, 4, 4}, 4, {}, {"false", "sleep 2", "true"})
assert(results[1].exitCode == 1 and results[2].timedOut and results[3].exitCode == 0)
-- A budget larger than the entire queue must also drain cleanly.
run({2, 4}, 20, {0.1, 0.1})
assert(#executeBatch({}, 8) == 0)
'''
            script = root / 'fixture.lua'
            script.write_text(lua)
            result = subprocess.run([LUAJIT, str(script), directory],
                                    capture_output=True, text=True, timeout=20)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_workers_and_compilation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'prefix/gkeyll/share').mkdir(parents=True)
            # Real make jobs, one deliberate failure and one tiny executable.
            (root / 'prefix/gkeyll/share/Makefile').write_text(
                'good:\n\tprintf "#!/bin/sh\\nexit 0\\n" > good\n\tchmod +x good\n'
                'bad:\n\t@echo deliberate compile failure\n\t@false\n')
            for name in ('good', 'bad'):
                (root / (name + '.c')).touch()
            lua = r'''
local root = arg[1]
local configVals = {source_dir=root, results_dir=root, prefix=root .. "/prefix"}
local runMode, runID = "cpu_serial", "fixture"
local messages = {}
local function log(msg) table.insert(messages, msg) end
local function verboseLog(msg) end
local function hasCmd(name) return os.execute("command -v " .. name .. " >/dev/null") == 0 end
local function mkdir(path) assert(os.execute("mkdir -p '" .. path .. "'") == 0) end
local lfs = {attributes=function(path)
   local f = io.open(path); if f then f:close(); return true end
end}
'''
            lua += fragment('local function shellQuote(', 'local function hasGkylOutput(')
            lua += '''
local function basename(path) return path:match("([^/]+)$") end
local function dirname(path) return path:match("(.+)/") end
local function stripext(path) return path:gsub("%.[^.]+$", "") end
'''
            lua += fragment('local function wrapWithTimeout(', '-- ---- Test classification predicates')
            lua += fragment('local function prepareCompileCommand(', '-- prepareLuaRun(')
            lua += fragment('local executeBatch\n', '-- Verify a complete prior compile')
            lua += fragment('executeBatch = function(', 'local function runLuaTest(')
            lua += r'''
local items = {}
for idx, duration in ipairs({0.8, 0.1, 0.1, 0.1}) do
   local dir = root .. "/job" .. idx
   mkdir(dir)
   items[idx] = {test={name="job" .. idx}, runDir=dir, cmd=wrapWithTimeout("sleep " .. duration, 0, dir)}
end
local completions = {}
local results = executeBatch(items, 2, function(item, result, completed, total)
   assert(completed == #completions + 1 and total == 4)
   table.insert(completions, item.test.name)
   if completed == 1 then
      assert(item.test.name == "job2", "results were not streamed in completion order")
      local f = assert(io.open(items[1].runDir .. "/_parallel_out.txt"))
      local raw = f:read("*a"); f:close()
      assert(not raw:find("__END__:", 1, true), "callback waited for slow job")
   end
end, true)
assert(#completions == 4)
local progress = table.concat(messages)
assert(progress:find("[C] [Worker 1] running job1", 1, true))
assert(progress:find("[C] [Worker 2] running job2", 1, true))
assert(progress:find("[C] [Worker 2] running job3", 1, true))
assert(progress:find("[C] [Worker 2] running job4", 1, true))
assert(not progress:find("Batch", 1, true))
for idx, result in ipairs(results) do
   assert(result.exitCode == 0 and result.runtm > 0)
end
local function interval(idx)
   local f = assert(io.open(items[idx].runDir .. "/_parallel_out.txt"))
   local raw = f:read("*a"); f:close()
   return tonumber(raw:match("__START__:([%d.]+)")), tonumber(raw:match("__END__:([%d.]+)"))
end
local start1, end1 = interval(1)
local start2, end2 = interval(2)
local start3, end3 = interval(3)
local start4, end4 = interval(4)
assert(start3 >= end2 and start3 < end1, "worker slot was not refilled")
assert(start4 >= end3, "worker limit exceeded")
items[1].cmd = "false; echo __EXIT__:$?"
items[2].cmd = wrapWithTimeout("sleep 2", 1, items[2].runDir)
results = executeBatch({items[1], items[2]}, 2)
assert(results[1].exitCode == 1)
assert(results[2].timedOut)
items[1].cmd = "exit 0" -- Missing exit marker must fail closed.
messages = {}
items[1].test.testType = "lua"
assert(executeBatch({items[1]}, 1, nil, true)[1].exitCode ~= 0)
progress = table.concat(messages)
assert(progress == "[Lua] running job1\n")
assert(items[1].workerLabel == "")
local cTests = {}
for _, name in ipairs({"good", "bad"}) do
   table.insert(cTests, {name=name, layer="fixture", src=root .. "/" .. name .. ".c"})
end
messages = {}
local ok, compiled = compile_c_regressions(cTests, 2)
assert(not ok and not compiled.good.compileFailed and compiled.bad.compileFailed)
assert(compiled.bad.compileLog:find("deliberate compile failure", 1, true))
progress = table.concat(messages)
assert(progress:find("[C] Compiled good\n", 1, true))
assert(progress:find("[C] COMPILE FAILED bad\n", 1, true))
assert(not progress:find("Batch", 1, true))
assert(not progress:find("running", 1, true))
'''
            script = root / 'fixture.lua'
            script.write_text(lua)
            result = subprocess.run([LUAJIT, str(script), directory],
                                    capture_output=True, text=True, timeout=20)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
