-- Gkyl ------------------------------------------------------------------------
--
-- Runs regression tests and compares with accepted results.
-- Each layer (moments -> vlasov -> gyrokinetics -> pkpm) 
-- has its own luareg/ (Lua regression tests) and creg/
-- (C regression test sources).
--
-- Lua tests:  discovered from source_dir/<layer>/luareg/rt_*.lua.
--             Run via 'gkeyll <file>' from a per-test scratch directory.
--
-- C tests:    discovered from source_dir/<layer>/creg/rt_*.c.
--             Each test is compiled on-the-fly using the installed
--             share/Makefile (processed from Makefile_for_ext_C_input at
--             'make install' time), then run from its scratch directory.
--             No separate 'make regression' step is required.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local lfs = require "lfs"
local uuid = require "Lib.UUID"

if GKYL_HAVE_SQLITE3 == false then
   -- can't run without SQLITE3
   print("Sorry, runregression needs Sqlite3. This executable was built without it.")
   return 1
end

local Logger = require "Lib.Logger"
local Time = require "Lib.Time"
local argparse = require "Lib.argparse"
local date = require "xsys.date"
local lume = require "Lib.lume"
local sql = require "sqlite3"

-- GKYL_OUT_PREFIX is required by some Lua modules; set to a safe default here
-- (each test run redirects output to its own scratch directory).
GKYL_OUT_PREFIX = lfs.currentdir() .. "/" .. "runregression"

local log = Logger { logToFile = true }
local verboseLog = function (msg) end -- default: no verbose output
local verboseLogger = function (msg) log(msg) end

-- ---- Layer definitions ------------------------------------------------------
-- The four physics layers that have regression tests. 'src' is the name of
-- the subdirectory under the gkeyll source root and build root.
local LAYERS = {
   { name = "moments",     src = "moments"     },
   { name = "vlasov",      src = "vlasov"      },
   { name = "gyrokinetic", src = "gyrokinetic" },
   { name = "pkpm",        src = "pkpm"        },
}

-- Fast lookup: set of valid layer names.
local VALID_LAYER_NAMES = {}
for _, L in ipairs(LAYERS) do VALID_LAYER_NAMES[L.name] = true end

-- ---- GPU build detection ----------------------------------------------------
-- Layers with production GPU (CUDA) code paths.  The moments layer has no GPU
-- support and is intentionally excluded.
local GPU_LAYERS = { vlasov = true, gyrokinetic = true, pkpm = true }

-- Global flag: is the current build GPU-enabled (CC=nvcc)?
-- Set once by detectGpuBuild() during loadConfigure().
local GPU_BUILD = false

-- Detect a GPU build by reading the installed config.mak.
-- Returns the new value of GPU_BUILD for convenience.
local function detectGpuBuild(prefix)
   local configMak = prefix .. "/gkeyll/share/config.mak"
   local f = io.open(configMak, "r")
   if f then
      for line in f:lines() do
         if line:match("^CC%s*=%s*nvcc") then
            GPU_BUILD = true
            break
         end
      end
      f:close()
   end
   return GPU_BUILD
end

-- ---- Per-test timeout command detection -------------------------------------
-- GNU 'timeout' is not part of macOS's standard utilities. It ships as
-- 'gtimeout' via 'brew install coreutils', or as 'timeout' on Linux and on
-- macOS with coreutils linked without the 'g' prefix.
-- We detect whichever is available at startup. If neither is found we fall
-- back to a Perl-based fork+alarm implementation (Perl is always present on
-- macOS and Linux).
local TIMEOUT_CMD
do
   local function hasCmd(name)
      return os.execute(
         string.format("which %s > /dev/null 2>&1", name)) == 0
   end
   if     hasCmd("timeout")  then TIMEOUT_CMD = "timeout"
   elseif hasCmd("gtimeout") then TIMEOUT_CMD = "gtimeout"
   else                           TIMEOUT_CMD = nil  -- use Perl fallback
   end
end

-- ---- Physical CPU count detection ------------------------------------------
-- Used by --jobs 0 to auto-detect how many tests to run concurrently.
-- Tries macOS sysctl (physical cores only), then Linux nproc, defaults to 1.
local function physicalCpuCount()
   local f = io.popen("sysctl -n hw.physicalcpu 2>/dev/null", "r")
   if f then
      local n = tonumber(f:read("*a"))
      f:close()
      if n and n > 0 then return n end
   end
   local f2 = io.popen("nproc 2>/dev/null", "r")
   if f2 then
      local n2 = tonumber(f2:read("*a"))
      f2:close()
      if n2 and n2 > 0 then return n2 end
   end
   return 1
end

-- ---- Layer pre-processing of command arguments ------------------------------
-- The argparse library does not natively support positional sub-sub-command
-- names that collide with free text (e.g. "run moments check"). We therefore
-- scan GKYL_COMMANDS_L before argparse sees it, extract any token that is a
-- known layer name, and store it in 'detectedLayer'. The remaining tokens are
-- passed to argparse as normal, so "run moments check" becomes "run check"
-- from argparse's perspective.
local detectedLayer = nil
local filteredCmds  = {}
for i = 1, #GKYL_COMMANDS_L do
   local arg = GKYL_COMMANDS_L[i]
   if VALID_LAYER_NAMES[arg] then
      -- First layer name wins; extras (if any) are silently dropped.
      if detectedLayer == nil then
         detectedLayer = arg
      end
   else
      filteredCmds[#filteredCmds + 1] = arg
   end
end

-- ---- Global counters --------------------------------------------------------
-- Per-layer pass/fail counts, keyed by layer name.
-- GPU counters track the GPU variant independently (only populated for
-- GPU-capable layers on a GPU build).
local layerCounts = {}
for _, L in ipairs(LAYERS) do
   layerCounts[L.name] = {
      total = 0, passed = 0, failed = 0,
      gpu_total = 0, gpu_passed = 0, gpu_failed = 0,
   }
end

-- ---- Configuration ----------------------------------------------------------
local isConfiguring = false
local configVals    = nil

-- Path of the configuration file written by 'configure' and read by 'run'.
-- Preferred location: <prefix>/gkeyll-results/runregression.config.lua, derived
-- from config.mak (same logic as prefix auto-detection in config_action).
-- Falls back to ~/runregression.config.lua for backwards compatibility.
local function computeConfFile()
   local gkeyllDir = GKYL_EXEC_PATH and GKYL_EXEC_PATH:match("^(.+)/bin$")
   if gkeyllDir then
      local mf = io.open(gkeyllDir .. "/share/config.mak", "r")
      if mf then
         for line in mf:lines() do
            local p = line:match("^PREFIX%s*=%s*(.+)%s*$")
            if p then
               mf:close()
               return p .. "/gkeyll-results/runregression.config.lua",
                      p .. "/gkeyll-results"
            end
         end
         mf:close()
      end
   end
   return os.getenv("HOME") .. "/runregression.config.lua", nil
end
local confFile, confFileResultsDir = computeConfFile()

-- Per-layer SQLite database connections and prepared statements.
-- Populated lazily by getLayerDB().
local layerDBs = {}

-- Per-layer ignore lists and MOAT lists, keyed by layer name.
-- Loaded during loadConfigure().
local ignoreTests = {}  -- ignoreTests[layerName] = {lua={...}, c={...}}
local moatTests   = {}  -- moatTests[layerName]   = {path, ...}

-- Unique ID and timestamp for this invocation.
local runID  = uuid()
local runDate = date(false):fmt("${iso}")

-- ---- Database schema --------------------------------------------------------
-- Two tables per layer database:
--
-- RegressionMeta stores one row per 'run' invocation:
--   guid               unique ID shared with RegressionData rows
--   tstamp             ISO-8601 timestamp of the run
--   GKYL_EXEC          path to the gkeyll executable used
--   GKYL_GIT_CHANGESET git commit hash
--   GKYL_BUILD_DATE    build date string
--   ntotal / npass / nfail  aggregate counts for this run (CPU)
--   gpu_build          1 if this was a GPU build, 0 otherwise
--   ngpu_pass / ngpu_fail  aggregate GPU-variant counts
--
-- RegressionData stores one row per individual test:
--   guid        links to a RegressionMeta row
--   name        test name (relative path / binary name)
--   test_type   'lua' or 'c'  (new in hierarchical design)
--   status      -4=compile_fail  -3=timeout  -2=create  -1=skip  0=fail  1=pass
--   runtime     wall-clock seconds (CPU variant)
--   runlog      captured stdout+stderr from the test
--   gpu_status  GPU variant: 1=pass, 0=fail, -1=skip, -3=timeout, -5=crash
--   gpu_runtime wall-clock seconds (GPU variant; 0 if not run)
--   cpu_gpu_diff  1=match, 0=differ, -1=not-applicable

local SCHEMA_SQL = [[
  drop table if exists RegressionMeta;
  create table RegressionMeta (
    guid               text,
    tstamp             text,
    GKYL_EXEC          text,
    GKYL_GIT_CHANGESET text,
    GKYL_BUILD_DATE    text,
    ntotal             integer,
    npass              integer,
    nfail              integer,
    gpu_build          integer,
    ngpu_pass          integer,
    ngpu_fail          integer
  );

  drop table if exists RegressionData;
  create table RegressionData (
    guid         text,
    name         text,
    test_type    text,
    status       integer,
    runtime      real,
    runlog       text,
    gpu_status   integer,
    gpu_runtime  real,
    cpu_gpu_diff integer
  );
]]

-- ---- Helper utilities -------------------------------------------------------

-- Recursively creates all directories in 'path'.
local function mkdir(path)
   local sep, pStr = package.config:sub(1,1), "/"
   for dir in path:gmatch("[^" .. sep .. "]+") do
      pStr = pStr .. dir .. sep
      lfs.mkdir(pStr)
   end
end

-- Splits a comma-separated string into a list of words.
local function splitList(listStr)
   local words = {}
   for w in listStr:gmatch('[^,%s]+') do
      table.insert(words, w)
   end
   return words
end

-- Returns the basename of a path (everything after the last '/').
local function basename(path)
   return path:match("([^/]+)$") or path
end

-- Returns the directory component of a path (e.g. "/a/b/c.lua" -> "/a/b").
local function dirname(path)
   return path:match("^(.+)/[^/]+$") or "."
end

-- Strips the extension from a filename (e.g. "rt_foo.lua" -> "rt_foo").
local function stripext(fn)
   return fn:match("^(.+)%..+$") or fn
end

-- Walks a directory tree recursively. Returns a coroutine that yields
-- (dir, filename, attr) for every entry.
local function dirtree(dir)
   assert(dir and dir ~= "", "directory parameter is missing or empty")
   if string.sub(dir, -1) == "/" then
      dir = string.sub(dir, 1, -2)
   end
   local function yieldtree(d)
      for fn in lfs.dir(d) do
         if fn ~= "." and fn ~= ".." then
            local fullNm = d .. "/" .. fn
            local attr = lfs.attributes(fullNm)
            if attr then
               coroutine.yield(d, fn, attr)
               if attr.mode == "directory" then
                  yieldtree(fullNm)
               end
            end
         end
      end
   end
   return coroutine.wrap(function() yieldtree(dir) end)
end

-- Wraps a shell command string with a per-test timeout.
-- 'innerCmd' must already contain any needed 2>&1 redirection.
-- Writes a small helper script to runDir to avoid quoting problems.
-- Appends '; echo __EXIT__:$?' so the caller can parse the exit code.
-- Returns exit code 124 (the POSIX convention used by GNU timeout) when the
-- test exceeds timeoutSecs seconds.
local function wrapWithTimeout(innerCmd, timeoutSecs, runDir)
   if not timeoutSecs or timeoutSecs <= 0 then
      return innerCmd .. "; echo __EXIT__:$?"
   end

   -- Write innerCmd to a tiny shell script so we don't have to escape it.
   local scriptPath = runDir .. "/_rr_cmd.sh"
   local sf = io.open(scriptPath, "w")
   sf:write("#!/bin/sh\n")
   sf:write(innerCmd .. "\n")
   sf:close()
   os.execute(string.format("chmod +x '%s'", scriptPath))

   if TIMEOUT_CMD then
      -- timeout / gtimeout available: straightforward.
      return string.format("%s %d '%s'; echo __EXIT__:$?",
         TIMEOUT_CMD, timeoutSecs, scriptPath)
   else
      -- Perl-based fork+alarm fallback (always available on macOS and Linux).
      -- Key: child calls setpgrp so that gkeyll (a grandchild launched by sh)
      -- is in the same process group. We then kill the entire group with
      -- kill(SIGNAL, -$pid), which reaches gkeyll even though it is not a
      -- direct child of the Perl process.
      local perlPath = runDir .. "/_rr_timeout.pl"
      local pf = io.open(perlPath, "w")
      pf:write(string.format([[#!/usr/bin/env perl
use strict;
my $secs   = %d;
my $script = '%s';
my $pid = fork // die "fork: $!\n";
if (!$pid) {
    # Put the child in its own process group so kill -PGID reaches all
    # grandchildren (e.g. the gkeyll process launched by the sh wrapper).
    setpgrp(0, 0);
    exec($script) or die "exec: $!\n";
}
local $SIG{ALRM} = sub {
    # Kill the entire process group: -$pid means PGID == $pid.
    kill 'TERM', -$pid;
    sleep 2;
    kill 'KILL', -$pid;
    exit 124;
};
alarm($secs);
waitpid($pid, 0);
exit($? >> 8);
]], timeoutSecs, scriptPath:gsub("'", "\\'")))
      pf:close()
      return string.format("perl '%s'; echo __EXIT__:$?", perlPath)
   end
end

-- ---- Test classification predicates ----------------------------------------

local function isLuaRegressionTest(fn)
   return string.find(fn, "/rt_[^/]+%.lua$") ~= nil
end

local function isLuaUnitTest(fn)
   return string.find(fn, "/test_[^/]-%.lua$") ~= nil
end

local function isCxxUnitTest(fn)
   if string.find(fn, "/ctest_[^/]+$") ~= nil then
      local attr = lfs.attributes(fn)
      if attr and attr.mode == "file" and string.sub(attr.permissions, 3, 3) == "x" then
         return true
      end
   end
   return false
end

-- ---- Database helpers -------------------------------------------------------

-- Returns (and caches) the DB connection and prepared statements for a layer.
-- The database lives at: results_dir/<layer>/regressiondb
local function getLayerDB(layerName)
   if layerDBs[layerName] then return layerDBs[layerName] end

   local dbPath = configVals.results_dir .. "/" .. layerName .. "/regressiondb"
   local conn   = sql.open(dbPath)

   -- insertRegressionData binds: guid, name, test_type, status, runtime, runlog,
   --                             gpu_status, gpu_runtime, cpu_gpu_diff
   local insertData = conn:prepare [[
     insert into RegressionData values (?, ?, ?, ?, ?, ?, ?, ?, ?)
   ]]

   -- insertRegressionMeta binds: guid, tstamp, exec, changeset, builddate,
   --                             ntotal, npass, nfail, gpu_build, ngpu_pass, ngpu_fail
   local insertMeta = conn:prepare [[
     insert into RegressionMeta values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
   ]]

   layerDBs[layerName] = {
      conn       = conn,
      insertData = insertData,
      insertMeta = insertMeta,
   }
   return layerDBs[layerName]
end

-- Insert a single test result row into the appropriate layer's database.
-- gpuStatus / gpuRuntime / cpuGpuDiff default to -1/0/-1 when not provided
-- (non-GPU layers or non-GPU builds).
local function insertRegressionData(layerName, guid, name, testType,
      status, runtm, runlog, gpuStatus, gpuRuntime, cpuGpuDiff)
   gpuStatus  = gpuStatus  or -1
   gpuRuntime = gpuRuntime or 0
   cpuGpuDiff = cpuGpuDiff or -1
   local db = getLayerDB(layerName)
   db.insertData:reset():bind(
      guid, name, testType, status, runtm, runlog,
      gpuStatus, gpuRuntime, cpuGpuDiff
   ):step()
end

-- Insert the run-level summary row for a layer.
-- gpuBuild / ngpuPass / ngpuFail default to 0 when not provided.
local function insertRegressionMeta(layerName, guid, tm,
      ntotal, npass, nfail, gpuBuild, ngpuPass, ngpuFail)
   gpuBuild = gpuBuild or 0
   ngpuPass = ngpuPass or 0
   ngpuFail = ngpuFail or 0
   local db = getLayerDB(layerName)
   db.insertMeta:reset():bind(
      guid, tm,
      GKYL_EXEC,
      GKYL_GIT_CHANGESET,
      GKYL_BUILD_DATE,
      ntotal, npass, nfail,
      gpuBuild, ngpuPass, ngpuFail
   ):step()
end

-- ---- Configure command ------------------------------------------------------
-- Creates the directory structure, per-layer databases, and configuration file.
-- Explanation of the directory layout:
--   gkeyll-results/
--     <layer>/
--       regressiondb         SQLite database for this layer
--       luareg-runs/         per-test scratch dirs for Lua tests
--       luareg-accepted/     accepted (baseline) outputs for Lua tests
--       creg-runs/           per-test scratch dirs for C tests
--       creg-accepted/       accepted (baseline) outputs for C tests
--
-- Separating 'runs' from 'accepted' makes it possible to inspect the most
-- recent test run at any time without disturbing the accepted baseline.
-- Including the input file (.lua or .c) in the run directory allows the
-- developer to diff the input if outputs change unexpectedly.
local function configure(prefix, mpiExec, sourceDir, args)

   -- Validate prefix directory.
   local prefixAttr = lfs.attributes(prefix)
   if not prefixAttr or prefixAttr.mode ~= "directory" then
      assert(false, string.format("Prefix '%s' is not a directory!", prefix))
   end
   log(string.format("Will write accepted results to %s/gkeyll-results ...\n", prefix))

   -- Validate source directory.
   local srcAttr = lfs.attributes(sourceDir)
   if not srcAttr or srcAttr.mode ~= "directory" then
      assert(false, string.format("Source directory '%s' is not a directory!", sourceDir))
   end

   -- Create the top-level results directory and per-layer subdirectories.
   local resultsDir = string.format("%s/gkeyll-results", prefix)
   mkdir(resultsDir)

   for _, layer in ipairs(LAYERS) do
      local ld = resultsDir .. "/" .. layer.name
      mkdir(ld)
      mkdir(ld .. "/luareg-runs")
      mkdir(ld .. "/luareg-accepted")
      mkdir(ld .. "/creg-runs")
      mkdir(ld .. "/creg-accepted")

      -- Create (or re-create) the SQLite database for this layer.
      -- The --drop-tables flag forces recreation even if the DB already exists.
      local dbPath  = ld .. "/regressiondb"
      local dbExist = lfs.attributes(dbPath)

      if args.drop_tables or dbExist == nil then
         log(string.format("Creating DB for layer '%s' at %s\n", layer.name, dbPath))
         local conn = sql.open(dbPath)
         -- DROP + CREATE ensures the schema is always consistent.
         -- test_type is the key new column: distinguishes Lua ('lua') from C ('c')
         -- regression tests within a single layer database.
         conn:exec(SCHEMA_SQL)
         conn:close()
      else
         log(string.format(
            "DB for layer '%s' already exists at %s (use --drop-tables to reset)\n",
            layer.name, dbPath))
      end
   end

   -- Write the configuration file. This file is read by all subsequent
   -- 'run', 'check', 'list', etc. invocations.
   -- 'prefix' is stored so that runCTest can locate share/Makefile for
   -- compiling C regression tests on-the-fly.
   local fn = io.open(confFile, "w")
   fn:write("return {\n")
   fn:write(string.format("  mpiExec     = \"%s\",\n", mpiExec))
   fn:write(string.format("  prefix      = \"%s\",\n", prefix))
   fn:write(string.format("  results_dir = \"%s\",\n", resultsDir))
   fn:write(string.format("  source_dir  = \"%s\",\n", sourceDir))
   fn:write("}\n")
   fn:close()
   log(string.format("Configuration written to %s\n", confFile))
end

-- ---- Load configuration -----------------------------------------------------
-- Reads the config file, validates it, opens global state (ignore lists, etc.).
local function loadConfigure(args)
   local f = loadfile(confFile)
   if not f then
      -- Try the legacy home-directory location as a fallback.
      local legacyPath = os.getenv("HOME") .. "/runregression.config.lua"
      f = loadfile(legacyPath)
      if f then
         log(string.format(
            "NOTE: config loaded from legacy location %s.\n"
            .. "Re-run 'runregression configure' to migrate it to gkeyll-results/.\n",
            legacyPath))
      else
         log("Regression tests not configured! Run 'runregression configure' first.\n")
         os.exit(1)
      end
   end
   configVals = f()

   -- Sanity-check that required keys are present.
   if not configVals.source_dir then
      log("Configuration is missing source_dir. "
         .. "Please re-run 'runregression configure' with --source-dir.\n")
      os.exit(1)
   end

   -- Backward compatibility: older config files lack 'prefix'. Derive it from
   -- GKYL_EXEC_PATH (e.g. ~/gkylsoft/gkeyll/bin -> ~/gkylsoft).
   if not configVals.prefix then
      local execParent = GKYL_EXEC_PATH:match("^(.+)/gkeyll/bin/?$")
      configVals.prefix = execParent or (os.getenv("HOME") .. "/gkylsoft")
      log(string.format(
         "NOTE: config lacks 'prefix'; derived as '%s'. "
         .. "Re-run 'configure' to make this permanent.\n", configVals.prefix))
   end

   if args.verbose then
      verboseLog = verboseLogger
   end

   -- Load per-layer ignore and MOAT lists.
   -- ignore_lua_tests.lua (in luareg/): return { tests = {...}, gpu = {...} }
   -- ignore_c_tests.lua   (in creg/):   return { tests = {...}, gpu = {...} }
   -- moat.lua (in luareg/): return {"rt_test1", "rt_test2", ...}
   -- Missing files are silently treated as empty lists.
   for _, layer in ipairs(LAYERS) do
      local srcBase    = configVals.source_dir .. "/" .. layer.src
      local luaIgnFile = srcBase .. "/luareg/ignore_lua_tests.lua"
      local cIgnFile   = srcBase .. "/creg/ignore_c_tests.lua"
      local luaMoatFile = srcBase .. "/luareg/moat_lua.lua"
      local cMoatFile   = srcBase .. "/creg/moat_c.lua"

      local luaIgn = { tests = {}, gpu = {} }
      local cIgn   = { tests = {}, gpu = {} }

      local gLua = loadfile(luaIgnFile)
      if gLua then
         local ok, loaded = pcall(gLua)
         if ok and type(loaded) == "table" then luaIgn = loaded end
      end

      local gC = loadfile(cIgnFile)
      if gC then
         local ok, loaded = pcall(gC)
         if ok and type(loaded) == "table" then cIgn = loaded end
      end

      ignoreTests[layer.name] = {
         lua     = luaIgn.tests or {},
         c       = cIgn.tests   or {},
         gpu_lua = luaIgn.gpu   or {},
         gpu_c   = cIgn.gpu     or {},
      }

      moatTests[layer.name] = { lua = {}, c = {} }
      local gLuaMoat = loadfile(luaMoatFile)
      if gLuaMoat then
         local ok, t = pcall(gLuaMoat)
         if ok and type(t) == "table" then moatTests[layer.name].lua = t end
      end
      local gCMoat = loadfile(cMoatFile)
      if gCMoat then
         local ok, t = pcall(gCMoat)
         if ok and type(t) == "table" then moatTests[layer.name].c = t end
      end

      if args.all then
         -- --all bypasses the ignore list for the test types that will actually run.
         -- Combined with --lua-only, only the Lua ignore list is cleared (C skipped anyway).
         -- Combined with --c-only, only the C ignore list is cleared (Lua skipped anyway).
         if not args.c_only then
            ignoreTests[layer.name].lua     = {}
            ignoreTests[layer.name].gpu_lua = {}
         end
         if not args.lua_only then
            ignoreTests[layer.name].c     = {}
            ignoreTests[layer.name].gpu_c = {}
         end
      end
   end

   -- Detect GPU build from the installed config.mak (CC=nvcc).
   detectGpuBuild(configVals.prefix)
   if GPU_BUILD then
      log("GPU build detected (CC=nvcc in config.mak)\n")
   else
      log("CPU-only build\n")
   end
end

-- ---- Test discovery ---------------------------------------------------------
-- Returns two lists of test descriptors:
--   luaTests: { {file, layer, name}, ... }   -- Lua .lua test files
--   cTests:   { {bin, src, layer, name}, ... } -- compiled C test binaries
--
-- 'activeLayers' is either nil (→ all layers) or a single layer name.
-- 'appPrefix' optionally filters Lua tests to those whose name starts with
-- "rt_<appPrefix>_" (mirrors the old --app option).
local function list_tests(activeLayers, args)
   local luaTests, cTests = {}, {}

   -- Track which --run-only entries matched at least one existing test file
   -- (in any scanned layer, Lua or C) so we can report unknown test names
   -- instead of silently running nothing.
   local runOnlyEntries, runOnlyFound = {}, {}
   if args.run_only then
      runOnlyEntries = splitList(args.run_only)
   end

   -- Determine which layers to scan.
   local layersToScan = {}
   if activeLayers then
      for _, L in ipairs(LAYERS) do
         if L.name == activeLayers then
            table.insert(layersToScan, L)
            break
         end
      end
   else
      layersToScan = LAYERS
   end

   for _, layer in ipairs(layersToScan) do
      -- ---- Lua tests --------------------------------------------------------
      local luaregDir = configVals.source_dir .. "/" .. layer.src .. "/luareg"
      local dirAttr   = lfs.attributes(luaregDir)
      if dirAttr and dirAttr.mode == "directory" then
         local function addLuaTest(fn)
            -- When an absolute path is given (e.g. via --run-only), only accept
            -- files that actually live under this layer's luareg directory.
            -- This prevents the per-layer loop from adding the same file to
            -- every layer in layersToScan.
            if string.sub(fn, 1, 1) == "/" then
               local luaregPath = "/" .. layer.src .. "/luareg/"
               if not string.find(fn, luaregPath, 1, true) then return end
            end
            if not isLuaRegressionTest(fn) then return end
            -- Skip tests in ignore list.  The lua table stores basenames
            -- (no path, no extension) so the file works on any machine.
            local ignLua = ignoreTests[layer.name] and ignoreTests[layer.name].lua or {}
            if lume.find(ignLua, stripext(basename(fn))) then
               if args.run_only then
                  log(string.format(
                     "NOTE: skipping '%s/luareg/%s': test is in the ignore list.\n",
                     layer.name, stripext(basename(fn))))
               end
               return
            end
            -- Derive a short display name from the absolute path for DB storage.
            local nm = layer.name .. "/luareg/" .. basename(fn)
            table.insert(luaTests, {
               file  = fn,
               layer = layer.name,
               name  = nm,
            })
         end

         if args.moat then
            -- Only run the MOAT Lua tests for this layer (basenames from moat_lua.lua).
            for _, t in ipairs(moatTests[layer.name].lua or {}) do
               local candidate = luaregDir .. "/" .. t .. ".lua"
               if lfs.attributes(candidate) then addLuaTest(candidate) end
            end
         elseif args.run_only then
            for _, ro in ipairs(runOnlyEntries) do
               local a = lfs.attributes(ro)
               if a then
                  runOnlyFound[ro] = true
                  if a.mode == "file" then
                     addLuaTest(ro)
                  elseif a.mode == "directory" then
                     for dir, fn, _ in dirtree(ro) do addLuaTest(dir .. "/" .. fn) end
                  end
               else
                  -- Bare test name (e.g. "rt_euler_sodshock"): search this layer's luareg/.
                  local candidate = luaregDir .. "/" .. ro .. ".lua"
                  if lfs.attributes(candidate) then
                     runOnlyFound[ro] = true
                     addLuaTest(candidate)
                  end
               end
            end
         else
            for dir, fn, _ in dirtree(luaregDir) do
               addLuaTest(dir .. "/" .. fn)
            end
         end
      else
         verboseLog(string.format(
            "Lua regression directory '%s' not found, skipping.\n", luaregDir))
      end

      -- ---- C tests ----------------------------------------------------------
      -- C regression sources live in source_dir/<layer>/creg/.
      -- Each rt_*.c file is compiled on-the-fly when the test runs; no
      -- pre-built binaries are required.
      local cregSrcDir = configVals.source_dir .. "/" .. layer.src .. "/creg"
      local cDirAttr   = lfs.attributes(cregSrcDir)
      if cDirAttr and cDirAttr.mode == "directory" then
         local function addCTest(fn)
            -- Layer-affinity guard for absolute paths (e.g. from --run-only).
            if string.sub(fn, 1, 1) == "/" then
               local cregPath = "/" .. layer.src .. "/creg/"
               if not string.find(fn, cregPath, 1, true) then return end
            end
            -- Must match rt_*.c pattern.
            if not string.match(fn, "/rt_[^/]+%.c$") then return end
            local testname = stripext(basename(fn))
            -- Skip tests in ignore list (C ignores are stored as basenames without .c).
            local ignC = ignoreTests[layer.name] and ignoreTests[layer.name].c or {}
            if lume.find(ignC, testname) then
               if args.run_only then
                  log(string.format(
                     "NOTE: skipping '%s/creg/%s': test is in the ignore list.\n",
                     layer.name, testname))
               end
               return
            end
            table.insert(cTests, {
               src      = fn,
               layer    = layer.name,
               layer_src = layer.src,   -- source subdir name (e.g. "gyrokinetic")
               name     = layer.name .. "/creg/" .. testname,
            })
         end

         if args.moat then
            -- Only run the MOAT C tests for this layer (basenames from moat_c.lua).
            for _, t in ipairs(moatTests[layer.name].c or {}) do
               local candidate = cregSrcDir .. "/" .. t .. ".c"
               if lfs.attributes(candidate) then addCTest(candidate) end
            end
         elseif args.run_only then
            for _, ro in ipairs(runOnlyEntries) do
               local a = lfs.attributes(ro)
               if a then
                  runOnlyFound[ro] = true
                  if a.mode == "file" then
                     addCTest(ro)
                  elseif a.mode == "directory" then
                     for dir, fn, attr in dirtree(ro) do
                        if attr.mode == "file" then addCTest(dir .. "/" .. fn) end
                     end
                  end
               else
                  -- Bare test name (e.g. "rt_10m_sodshock"): search this layer's creg/.
                  local candidate = cregSrcDir .. "/" .. ro .. ".c"
                  if lfs.attributes(candidate) then
                     runOnlyFound[ro] = true
                     addCTest(candidate)
                  end
               end
            end
         else
            for fn in lfs.dir(cregSrcDir) do
               if string.match(fn, "^rt_.+%.c$") then
                  addCTest(cregSrcDir .. "/" .. fn)
               end
            end
         end
      else
         verboseLog(string.format(
            "C regression directory '%s' not found, skipping.\n", cregSrcDir))
      end
   end

   -- Fail with a clear message if any --run-only entry matched no test in the
   -- scanned layers (e.g. 'rt_gk_sheath_3x2v' given when the actual tests are
   -- 'rt_gk_sheath_3x2v_p1' and 'rt_gk_sheath_3x2v_p1_cons').
   local missing = {}
   for _, ro in ipairs(runOnlyEntries) do
      if not runOnlyFound[ro] then table.insert(missing, ro) end
   end
   if #missing > 0 then
      local scannedNames = {}
      for _, L in ipairs(layersToScan) do
         table.insert(scannedNames, L.name)
      end
      for _, ro in ipairs(missing) do
         log(string.format(
            "ERROR: no regression test named '%s' found in layer(s): %s.\n",
            ro, table.concat(scannedNames, ", ")))
         -- Suggest tests whose names contain the requested name.
         local suggestions, seen = {}, {}
         for _, layer in ipairs(layersToScan) do
            for _, sub in ipairs({ "luareg", "creg" }) do
               local dirPath = configVals.source_dir .. "/" .. layer.src .. "/" .. sub
               local attr = lfs.attributes(dirPath)
               if attr and attr.mode == "directory" then
                  for fn in lfs.dir(dirPath) do
                     if string.match(fn, "^rt_.+%.lua$") or string.match(fn, "^rt_.+%.c$") then
                        local nm = stripext(fn)
                        local key = layer.name .. "/" .. sub .. "/" .. nm
                        if not seen[key] and string.find(nm, ro, 1, true) then
                           seen[key] = true
                           table.insert(suggestions, key)
                        end
                     end
                  end
               end
            end
         end
         if #suggestions > 0 then
            table.sort(suggestions)
            log("Did you mean one of these?\n")
            for _, s in ipairs(suggestions) do log("   " .. s .. "\n") end
         end
      end
      log("Use 'gkeyll runregression list' to see all available tests.\n")
      os.exit(1)
   end

   return luaTests, cTests
end

-- ---- Test runners -----------------------------------------------------------

-- Runs a single Lua regression test.
-- Workflow:
--   1. Create a per-test scratch directory under luareg-runs/<name>/
--   2. Copy the .lua input file there (so the developer can diff it later)
--   3. Remove stale .gkyl output files from the scratch directory
--   4. Run gkeyll <abs_path_to_lua> from within the scratch directory.
--      Without the -d flag, GKYL_OUT_PREFIX = basename(lua_file_without_ext),
--      so all .gkyl output lands in the CWD = scratch directory.
-- Returns (runtm, runlog, runDir).
-- timeoutSecs: 0 or nil = no limit; > 0 = kill test after that many seconds.
-- mode: "cpu" = force CPU via -G on a GPU build; "gpu" = use GPU (default on
--       GPU build, no extra flag needed); nil = legacy behaviour (no mode flag).
-- Returns: runtm, runlog, runDir, timedOut (boolean).
-- ---- Parallel execution infrastructure -------------------------------------
-- These three functions underpin both the serial (--jobs 1) and parallel
-- (--jobs N) execution paths.
--
-- prepareLuaRun / prepareCRun: extract all setup work (mkdir, copy, compile)
--   from runLuaTest / runCTest without calling io.popen for the test itself.
--   They return a table describing the ready-to-run command.
--
-- executeBatch: given a list of such tables, runs them all concurrently using
--   shell background jobs and collects the results.

-- Compiles a single C regression test in its scratch directory.
-- Copies test.src and the installed share/Makefile into scratchDir, then
-- runs 'make <testname>' there.
-- Returns: ok (boolean), compileLog (string).
local function compileCTest(test, scratchDir)
   local testname      = stripext(basename(test.src))
   local shareMakefile = configVals.prefix .. "/gkeyll/share/Makefile"

   if not lfs.attributes(shareMakefile) then
      return false, string.format(
         "share/Makefile not found at '%s'.\n"
         .. "Ensure 'make install' has been run for the current build.\n",
         shareMakefile)
   end

   -- Copy the C source and the installed Makefile into the scratch dir.
   os.execute(string.format("cp -f '%s' '%s/'", test.src, scratchDir))
   os.execute(string.format("cp -f '%s' '%s/'", shareMakefile, scratchDir))
   -- Copy rt_arg_parse.h from the same creg/ directory so the local -I.
   -- takes precedence over the installed header.  This ensures any updates to
   -- rt_arg_parse.h in the source tree (e.g. new fields) are picked up without
   -- requiring a full 'make install'.
   local argParseH = dirname(test.src) .. "/rt_arg_parse.h"
   if lfs.attributes(argParseH) then
      os.execute(string.format("cp -f '%s' '%s/'", argParseH, scratchDir))
   end

   -- Compile: 'make <testname>' from the scratch directory.
   local compileCmd = string.format(
      "cd '%s' && make '%s' 2>&1; echo __COMPILE_EXIT__:$?", scratchDir, testname)
   local proc    = io.popen(compileCmd, "r")
   local rawOut  = proc:read("*a")
   proc:close()
   local exitCode   = tonumber(rawOut:match("__COMPILE_EXIT__:(%d+)%s*$")) or 1
   local compileLog = rawOut:gsub("\n?__COMPILE_EXIT__:%d+%s*$", "")
   return (exitCode == 0), compileLog
end

-- prepareLuaRun(test, timeoutSecs, mode) → table
-- Sets up the scratch dir and builds the run command, but does NOT execute it.
-- Returns {cmd, runDir, test} normally, or {mpiSkip=true, runDir, test} for
-- MPI-only tests (numProc set) that cannot be run via this framework.
local function prepareLuaRun(test, timeoutSecs, mode)
   local testBasename = stripext(basename(test.file))
   local runDir = configVals.results_dir .. "/" .. test.layer
      .. "/luareg-runs/" .. testBasename

   mkdir(runDir)
   os.execute(string.format("cp -f '%s' '%s/'", test.file, runDir))
   os.execute(string.format("rm -f '%s'/*.gkyl 2>/dev/null", runDir))

   local opts = {}
   local fh = io.open(test.file, "r")
   if fh then
      local line1 = fh:read()
      if line1 and string.find(line1, "--!") then
         local r = string.sub(line1, 4, -1)
         opts = loadstring("return " .. r)()
      end
      fh:close()
   end

   if opts.numProc then
      return { mpiSkip = true, runDir = runDir, test = test }
   end

   local gkylExec = GKYL_EXEC_PATH .. "/gkeyll"
   local modeFlag = ""
   if mode == "cpu" and GPU_BUILD then modeFlag = " -G" end
   local innerCmd = string.format(
      "cd '%s' && '%s' '%s'%s 2>&1", runDir, gkylExec, test.file, modeFlag)
   local cmd = wrapWithTimeout(innerCmd, timeoutSecs or 0, runDir)

   return { cmd = cmd, runDir = runDir, test = test }
end

-- prepareCRun(test, timeoutSecs, mode, skipCompile) → table
-- Creates the scratch dir, compiles the C test (unless skipCompile), creates
-- the layer source symlink so tests can find data files by relative path, 
-- and builds the run command.
-- Returns {compileFailed=true, ...} on compile failure (no cmd field).
-- Returns {compileFailed=false, cmd, runDir, test, compileLog, compileSecs} on success.
local function prepareCRun(test, timeoutSecs, mode, skipCompile)
   local testname = stripext(basename(test.src))
   local runDir   = configVals.results_dir .. "/" .. test.layer
      .. "/creg-runs/" .. testname

   mkdir(runDir)
   os.execute(string.format("rm -f '%s'/*.gkyl 2>/dev/null", runDir))

   local compileLog  = ""
   local compileSecs = 0
   if not skipCompile then
      log(string.format("\n[C]   Compiling %s ...\n", test.name))
      local tmComp   = Time.clock()
      local compileOk
      compileOk, compileLog = compileCTest(test, runDir)
      compileSecs = Time.clock() - tmComp

      if not compileOk then
         log(string.format("... COMPILE FAILED in %g sec\n", compileSecs))
         verboseLog(compileLog)
         return {
            compileFailed = true,
            runDir        = runDir,
            test          = test,
            compileLog    = compileLog,
            compileSecs   = compileSecs,
         }
      end
      verboseLog(compileLog)
   end

   -- Symlink the layer source dir so tests can find data files by relative path.
   -- Always recreate the link so it tracks the currently configured source_dir
   -- (a stale link from a previous configure would point at the old tree).
   -- Remove the old link first: 'ln -sf' onto a symlink-to-directory would
   -- otherwise create the new link *inside* the old target.
   -- Also uses only rm -f to avoid accidentally deleting a real directory
   -- if for whatever reason a directory of the same name (vlasov or gyrokinetic)
   -- exists at the location we are trying to make the symlink. 
   if test.layer_src then
      local layerSrcPath = configVals.source_dir .. "/" .. test.layer_src
      local symlinkPath  = runDir .. "/" .. test.layer_src
      if lfs.attributes(layerSrcPath, "mode") == "directory" then
         os.execute(string.format("rm -f '%s' && ln -s '%s' '%s' 2>/dev/null",
            symlinkPath, layerSrcPath, symlinkPath))
      end
   end

   local binPath  = runDir .. "/" .. testname
   local gpuFlag  = (mode == "gpu") and " -g" or ""
   local innerCmd = string.format("cd '%s' && '%s'%s 2>&1", runDir, binPath, gpuFlag)
   local cmd      = wrapWithTimeout(innerCmd, timeoutSecs or 0, runDir)

   return {
      compileFailed = false,
      cmd           = cmd,
      runDir        = runDir,
      test          = test,
      compileLog    = compileLog,
      compileSecs   = compileSecs,
   }
end

-- executeBatch(items) → list of {runtm, runlog, timedOut}
-- Runs all items concurrently via shell background jobs.  Each item must have
-- {cmd, runDir}.  The function blocks until every job in the batch finishes.
--
-- Per-item timing uses __START__:epoch / __END__:epoch markers written around
-- each command's execution; exit status comes from wrapWithTimeout's existing
-- __EXIT__:N marker.  All output (stdout + stderr) goes to
-- runDir/_parallel_out.txt.  The coordinator script is placed in results_dir
-- and rewritten each call; per-item scripts go in their own runDir.
local function executeBatch(items)
   if #items == 0 then return {} end

   -- Step 1: write a per-item wrapper script so we never have to embed
   -- arbitrary command strings inside the coordinator (avoids quoting issues).
   for _, item in ipairs(items) do
      local sf = io.open(item.runDir .. "/_rr_batch_item.sh", "w")
      sf:write("#!/bin/sh\n")
      sf:write("echo __START__:$(date +%s)\n")
      -- item.cmd already ends with '; echo __EXIT__:$?' from wrapWithTimeout.
      sf:write(item.cmd .. "\n")
      sf:write("echo __END__:$(date +%s)\n")
      sf:close()
   end

   -- Step 2: write the batch coordinator script.
   -- Use the first item's runDir (not results_dir) so concurrent runregression
   -- processes (e.g. different layers running in parallel) never share a file.
   local coordPath = items[1].runDir .. "/_rr_batch_coordinator.sh"
   local cf = io.open(coordPath, "w")
   cf:write("#!/bin/sh\n")
   for _, item in ipairs(items) do
      local itemScript = item.runDir .. "/_rr_batch_item.sh"
      local outFile    = item.runDir .. "/_parallel_out.txt"
      cf:write(string.format("sh '%s' > '%s' 2>&1 &\n", itemScript, outFile))
   end
   cf:write("wait\n")
   cf:close()

   -- Step 3: run the coordinator (blocks until all background jobs finish).
   os.execute(string.format("sh '%s'", coordPath))

   -- Step 4: collect results.
   local results = {}
   for _, item in ipairs(items) do
      local rf = io.open(item.runDir .. "/_parallel_out.txt", "r")
      local raw = rf and rf:read("*a") or ""
      if rf then rf:close() end

      local startEpoch = tonumber(raw:match("__START__:(%d+)"))
      local endEpoch   = tonumber(raw:match("__END__:(%d+)"))
      local runtm = (startEpoch and endEpoch) and (endEpoch - startEpoch) or 0

      -- Strip timing markers first so they don't interfere with EXIT parsing.
      local stripped = raw
         :gsub("\n?__START__:%d+\n?", "\n")
         :gsub("\n?__END__:%d+\n?",   "\n")
      local exitCode = tonumber(stripped:match("__EXIT__:(%d+)%s*$")) or 0
      local runlog   = stripped:gsub("\n?__EXIT__:%d+%s*$", "")

      table.insert(results, {
         runtm    = runtm,
         runlog   = runlog,
         timedOut = (exitCode == 124),
      })
   end

   return results
end

local function runLuaTest(test, timeoutSecs, mode)
   local modeTag = (mode == "gpu") and "[GPU-Lua]" or "[Lua]"
   log(string.format("\n%s Running %s ...\n", modeTag, test.name))

   local prep = prepareLuaRun(test, timeoutSecs, mode)
   if prep.mpiSkip then
      log(string.format("**** NOT RUNNING PARALLEL TEST %s\n", test.name))
      return 0, "", prep.runDir, false
   end

   local results = executeBatch({ prep })
   local r = results[1]

   if r.timedOut then
      log(string.format("... TIMED OUT after %g sec\n", r.runtm))
   else
      log(string.format("... completed in %g sec\n", r.runtm))
   end
   verboseLog(r.runlog)

   return r.runtm, r.runlog, prep.runDir, r.timedOut
end

-- Runs a single C regression test (thin wrapper over prepareCRun + executeBatch).
-- prepareCRun handles: scratch dir creation, stale-file removal, compilation, symlink 
-- setup so tests can find data files by relative path, and run-command construction.
-- mode: "gpu" = append '-g' to the binary invocation; nil/omitted = CPU.
-- skipCompile: true = binary already exists (GPU re-run); skip compile phase.
-- keepBinary: true = do not delete the binary after running (caller handles it).
-- Returns: runtm, runlog, runDir, timedOut (bool), compileFailed (bool).
local function runCTest(test, timeoutSecs, mode, skipCompile, keepBinary)
   -- For GPU re-runs (skipCompile=true), prepareCRun does no logging, so we
   -- announce the test here.  For normal CPU runs, prepareCRun logs
   -- "\n[C]   Compiling X ..." which serves as the announcement.
   if skipCompile then
      local modeTag = (mode == "gpu") and "[GPU-C]" or "[C]"
      log(string.format("\n%s   Running %s ...\n", modeTag, test.name))
   end

   local prep = prepareCRun(test, timeoutSecs, mode, skipCompile)
   if prep.compileFailed then
      return prep.compileSecs, "COMPILE FAILED:\n" .. prep.compileLog,
             prep.runDir, false, true
   end

   local results = executeBatch({ prep })
   local r       = results[1]
   local testname = stripext(basename(test.src))
   local runDir   = prep.runDir

   -- Remove compiled binary (keep .c source, Makefile, .gkyl outputs).
   -- When skipCompile or keepBinary is true the caller handles cleanup.
   if not skipCompile and not keepBinary then
      os.execute(string.format("rm -f '%s/%s' '%s/%s.d' 2>/dev/null",
         runDir, testname, runDir, testname))
      os.execute(string.format("rm -rf '%s/%s.dSYM' 2>/dev/null",
         runDir, testname))
   end

   if r.timedOut then
      log(string.format("... TIMED OUT after %g sec\n", r.runtm))
   else
      log(string.format("... completed in %g sec\n", r.runtm))
   end
   verboseLog(r.runlog)

   return r.runtm, (prep.compileLog or "") .. "\n" .. r.runlog,
          runDir, r.timedOut, false
end

-- ---- File comparison --------------------------------------------------------
-- Returns a shortened version of an absolute path for display purposes.
-- Strips the results_dir prefix so paths display as e.g.
-- "moments/creg-accepted/rt_10m_sodshock/rt_10m_sodshock-0.gkyl".
local function shortPath(p)
   local rd = configVals.results_dir .. "/"
   if p:sub(1, #rd) == rd then return p:sub(#rd + 1) end
   return p
end

-- Compares two .gkyl files (field data or dynvector) element-by-element.
-- absTol / relTol: optional absolute and relative tolerance thresholds.
--   Defaults to 1e-12 for both (strict CPU comparison).
--   Pass looser values (e.g. 1e-7) for GPU-vs-accepted comparisons.
-- File comparison uses the G0.Zero Lua-C API registered by
-- gkyl_zero_lw_openlibs at startup (core/apps/zero_lw.c). 
local function compareFiles(f1, f2, absTol, relTol)
   absTol = absTol or 1e-12
   relTol = relTol or 1e-12
   verboseLog(string.format("  Comparing %s\n             %s ...\n", shortPath(f1), shortPath(f2)))
   if not lfs.attributes(f1) or not lfs.attributes(f2) then
      verboseLog(string.format(
         "    ... files %s and/or %s do not exist!\n", shortPath(f1), shortPath(f2)))
      return false
   end

   -- Wrap the C-level comparison calls in pcall.  Some .gkyl files can
   -- trigger Lua-level errors in the G0.Zero functions; pcall prevents
   -- those from aborting the entire regression run.  Note: C-level
   -- assertion failures (abort/SIGABRT) cannot be caught by pcall —
   -- tests that trigger those must be added to the ignore list.
   -- pcall returns (true, result, detail) on success or (false, errmsg) on throw.
   local ok, result, detail = pcall(function()
      local f1type = G0.Zero.gkylFileType(f1)
      local f2type = G0.Zero.gkylFileType(f2)

      if f1type ~= f2type then
         verboseLog(string.format(
            "    ... type mismatch: %s vs %s\n", f1type, f2type))
         return false, string.format("type mismatch: %s vs %s", f1type, f2type)
      end

      if f1type == "dynvector" then
         local diff = G0.Zero.dynvecDiff(f1, f2)
         if not diff.is_compatible then
            verboseLog("    ... dynvector files not compatible (size/type mismatch or read failure)\n")
            return false, "dynvec not compatible"
         end
         -- Combined tolerance: fail only when BOTH absolute and relative
         -- thresholds are exceeded.  This avoids false failures for near-zero
         -- values (tiny abs but large rel%) and for large values (large abs but
         -- tiny rel%).  When max_abs_diff is 0, the abs check is false so the
         -- condition short-circuits (handles the 0/0 → DBL_MAX rel case).
         if diff.max_abs_diff > absTol and diff.max_rel_diff > relTol then
            verboseLog(string.format(
               "    ... dynvec max abs diff %g (tol %g), max rel diff %g (tol %g)\n",
               diff.max_abs_diff, absTol, diff.max_rel_diff, relTol))
            return false, string.format("dynvec max_abs=%.3g max_rel=%.3g",
               diff.max_abs_diff, diff.max_rel_diff)
         end
         if diff.tm_max_abs_diff > 1e-10 then
            verboseLog(string.format(
               "    ... dynvec timestamp max abs diff %g (informational)\n", diff.tm_max_abs_diff))
         end
         return true
      end

      if f1type == "block-topology" then
         local equal = G0.Zero.blockTopoCmp(f1, f2)
         if not equal then
            verboseLog("    ... block topology mismatch\n")
            return false, "topology mismatch"
         end
         return true
      end

      -- arrayNewFromFile returns (nil, nil) on failure rather than throwing.
      local g1, a1 = G0.Zero.arrayNewFromFile(f1)
      local g2, a2 = G0.Zero.arrayNewFromFile(f2)
      if not g1 or not g2 then
         verboseLog(string.format(
            "    ... skipping %s (unsupported file format)\n", shortPath(f1)))
         return true
      end

      if not G0.Zero.rectGridCmp(g1, g2) then
         return false, "grid mismatch"
      end

      local nghost = { 0, 0, 0, 0, 0, 0, 0 }
      local r1, er1 = G0.Zero.createGridRanges(g1, nghost)
      local r2, er2 = G0.Zero.createGridRanges(g2, nghost)

      local diff = G0.Zero.arrayDiff(a1, a2, r1)

      if not diff.is_compatible then return false, "incompatible arrays" end
      -- Combined tolerance: fail only when BOTH absolute and relative thresholds
      -- are exceeded.  Near-zero values naturally have large relative differences
      -- (e.g. 1e-15 vs -1e-15 → rel=200%) but negligible absolute differences;
      -- large values can have large absolute differences but tiny relative
      -- differences.  Failing on either alone produces false positives.
      -- When max_abs_diff is 0, the condition short-circuits safely (handles the
      -- 0/0 → DBL_MAX rel case from gkyl_array_diff).
      if diff.max_abs_diff > absTol and diff.max_rel_diff > relTol then
         verboseLog(string.format(
            "    ... max abs diff %g (tol %g), max rel diff %g (tol %g)\n",
            diff.max_abs_diff, absTol, diff.max_rel_diff, relTol))
         return false, string.format("max_abs=%.3g max_rel=%.3g",
            diff.max_abs_diff, diff.max_rel_diff)
      end

      return true
   end)

   if not ok then
      -- result holds the error message when pcall catches a throw.
      verboseLog(string.format(
         "    ... comparison CRASHED: %s\n", tostring(result)))
      log(string.format(
         "WARNING: comparison crashed for %s (C-level error: %s)\n",
         basename(f1), tostring(result)))
      return false, "crash: " .. tostring(result)
   end
   return result, detail
end

-- ---- Post-run actions: create and check -------------------------------------
-- Each action takes (test, runDir, testType) where testType is 'lua' or 'c'.
-- 'create' copies the run's .gkyl output to the accepted directory, making
-- it the new baseline. Status -2 is recorded in the database.
-- 'check' compares the current run's output against the accepted baseline.
-- Status 1 = pass, 0 = fail.

-- Returns the accepted results directory for the given test and type.
local function acceptedDir(test, testType)
   -- testType is 'lua' or 'c'; subdirectory names are 'luareg-accepted' / 'creg-accepted'.
   local sub = (testType == "lua") and "luareg-accepted" or "creg-accepted"
   local nm  = (testType == "lua")
      and stripext(basename(test.file))
      or  stripext(basename(test.src))
   return configVals.results_dir .. "/" .. test.layer .. "/" .. sub .. "/" .. nm
end

local function create_action(test, runDir, testType)
   local aDir = acceptedDir(test, testType)
   log(string.format("... saving accepted results to %s ...\n", aDir))
   mkdir(aDir)
   -- Copy all .gkyl output files from the scratch directory to the accepted dir.
   os.execute(string.format("cp -f '%s'/*.gkyl '%s/' 2>/dev/null", runDir, aDir))
   return -2
end

local function check_action(test, runDir, testType, absTol, relTol)
   local aDir = acceptedDir(test, testType)
   -- For Lua tests the output prefix matches the input basename (e.g. rt_euler_sodshock).
   -- For C tests the app's .name field sets the prefix, which we have standardised to
   -- match the source basename (e.g. rt_10m_sodshock), so all .gkyl files in the run
   -- directory belong to this test and we compare them all.
   local testPrefix = (testType == "lua") and stripext(basename(test.file)) or nil

   local passed, count = true, 0
   local failedFiles = {}
   -- Walk the accepted directory (ground truth) and verify each file is present
   -- in the run directory and matches within tolerance. Walking accepted (not run)
   -- means a test that crashes before writing its last frame is correctly flagged:
   -- the missing run file is detected rather than silently skipped.
   if lfs.attributes(aDir) then
      for fn in lfs.dir(aDir) do
         local isGkyl  = (string.sub(fn, -5) == ".gkyl")
         local inScope = isGkyl and (
            testPrefix == nil or string.find(fn, "^" .. testPrefix) ~= nil)
         if inScope then
            count = count + 1
            local accFile = aDir   .. "/" .. fn
            local runFile = runDir .. "/" .. fn
            if not lfs.attributes(runFile) then
               verboseLog(string.format("  MISSING run file: %s\n", fn))
               table.insert(failedFiles, fn .. "  [MISSING]")
               passed = false
            else
               local ok, detail = compareFiles(accFile, runFile, absTol, relTol)
               if not ok then
                  local entry = fn .. "  [DIFF]"
                  if detail and detail ~= "" then entry = entry .. "  " .. detail end
                  table.insert(failedFiles, entry)
               end
               passed = passed and ok
            end
         end
      end
   end

   -- count == 0 means either no accepted files exist (create not yet run) or the
   -- accepted directory itself is absent. Either way the check cannot pass.
   if count == 0 then passed = false end

   -- Build the comparison log (stored in the DB so queryrdb can surface it later).
   local checkLog = ""
   if passed then
      layerCounts[test.layer].passed = layerCounts[test.layer].passed + 1
      log("... passed.\n")
   else
      layerCounts[test.layer].failed = layerCounts[test.layer].failed + 1
      log(string.format("... %s FAILED!\n", test.name))
      if #failedFiles > 0 then
         -- Always log failing file names with diff magnitudes (not gated on verbose).
         log(string.format("  Failing files (%d):\n", #failedFiles))
         for _, ff in ipairs(failedFiles) do
            log(string.format("    %s\n", ff))
         end
         log("  Legend: [DIFF] values exceed tolerance  [MISSING] file not produced by run\n")
         -- Build checkLog for DB storage so queryrdb --test N surfaces magnitudes.
         checkLog = "--- Comparison failures ---\n"
            .. table.concat(failedFiles, "\n")
         -- Write a machine-parseable failures file to the run directory.
         local ffPath = runDir .. "/_rr_failures.txt"
         local ffFile = io.open(ffPath, "w")
         if ffFile then
            ffFile:write(string.format("test: %s\n", test.name))
            ffFile:write(string.format("accepted: %s\n", aDir))
            ffFile:write(string.format("run: %s\n", runDir))
            ffFile:write(string.format("failures: %d\n", #failedFiles))
            for _, ff in ipairs(failedFiles) do ffFile:write(ff .. "\n") end
            ffFile:close()
         end
      end
   end
   return passed and 1 or 0, checkLog
end

-- ---- GPU comparison helpers -------------------------------------------------
-- These mirror the file-walking logic in check_action but do NOT update
-- layerCounts (GPU counts are tracked separately in run_action).

-- Compare GPU run output against accepted baselines using the given tolerance.
-- Returns: 1 = pass, 0 = fail.
local function gpuCheck_vs_accepted(test, runDir, testType, gpuTol)
   local aDir = acceptedDir(test, testType)
   local testPrefix = (testType == "lua") and stripext(basename(test.file)) or nil
   local passed, count = true, 0
   if lfs.attributes(aDir) then
      for fn in lfs.dir(aDir) do
         local isGkyl  = (string.sub(fn, -5) == ".gkyl")
         local inScope = isGkyl and (
            testPrefix == nil or string.find(fn, "^" .. testPrefix) ~= nil)
         if inScope then
            count = count + 1
            local accFile = aDir   .. "/" .. fn
            local runFile = runDir .. "/" .. fn
            if not lfs.attributes(runFile) then
               passed = false
            else
               local ok = compareFiles(accFile, runFile, gpuTol, gpuTol)
               passed = passed and ok
            end
         end
      end
   end
   if count == 0 then passed = false end
   return passed and 1 or 0
end

-- Compare CPU output (in cpuDir) against GPU output (in gpuDir).
-- Returns: 1 = outputs match within tolerance, 0 = outputs differ.
local function compareCpuGpu(test, cpuDir, gpuDir, testType, gpuTol)
   local testPrefix = (testType == "lua") and stripext(basename(test.file)) or nil
   local passed, count = true, 0
   for fn in lfs.dir(gpuDir) do
      local isGkyl   = (string.sub(fn, -5) == ".gkyl")
      local inScope  = isGkyl and (
         testPrefix == nil or string.find(fn, "^" .. testPrefix) ~= nil)
      if inScope then
         count = count + 1
         local cpuFile = cpuDir .. "/" .. fn
         local gpuFile = gpuDir .. "/" .. fn
         if lfs.attributes(cpuFile) then
            local ok = compareFiles(cpuFile, gpuFile, gpuTol, gpuTol)
            passed = passed and ok
         else
            -- CPU file missing: cannot compare.
            passed = false
         end
      end
   end
   if count == 0 then passed = false end
   return passed and 1 or 0
end

-- ---- Unit test helpers ----------------------------
local function runLuaUnitTest(test)
   log(string.format("\nRunning unit test %s ...\n", test))
   local gkylExec = GKYL_EXEC_PATH .. "/gkeyll"
   local f = io.popen(string.format("'%s' '%s'", gkylExec, test), "r")
   local outPut = f:read("*a")
   f:close()
   if string.find(outPut, "PASSED") then
      log("... passed.\n")
   else
      if string.find(outPut, "without CUDA") then
         log(string.format("... %s skipped (no CUDA).\n", test))
      else
         log(string.format("... %s FAILED!\n", test))
      end
   end
end

local function runCxxUnitTest(test)
   log(string.format("\nRunning C++ unit test %s ...\n", test))
   local f = io.popen(string.format("'%s'", test), "r")
   local outPut = f:read("*a")
   f:close()
   if string.find(outPut, "FAILED") then
      log(string.format("... %s FAILED!\n", test))
   else
      log("... passed.\n")
   end
end

local function list_unit_tests(args)
   local luaUnitTests, cxxUnitTests = {}, {}
   -- Unit tests still live in the source tree (parallel structure to regression).
   local layersToScan = {}
   if detectedLayer then
      for _, L in ipairs(LAYERS) do
         if L.name == detectedLayer then
            table.insert(layersToScan, L); break
         end
      end
   else
      layersToScan = LAYERS
   end
   for _, layer in ipairs(layersToScan) do
      local unitDir = configVals.source_dir .. "/" .. layer.src .. "/unit"
      if lfs.attributes(unitDir) then
         for dir, fn, _ in dirtree(unitDir) do
            local fullNm = dir .. "/" .. fn
            if isLuaUnitTest(fullNm) then
               table.insert(luaUnitTests, fullNm)
            elseif isCxxUnitTest(fullNm) then
               table.insert(cxxUnitTests, fullNm)
            end
         end
      end
   end
   return luaUnitTests, cxxUnitTests
end

-- ---- Ignore-list updater ----------------------------------------------------
-- Called after a run with --timeout when some tests exceeded the limit.
-- Merges timed-out names into the layer's per-suite ignore files and writes
-- them back so subsequent runs automatically skip those tests.
--
-- Lua timeouts  → luareg/ignore_lua_tests.lua  (tests / gpu keys)
-- C timeouts    → creg/ignore_c_tests.lua       (tests / gpu keys)
--
-- Re-reads each file from disk before writing: manual edits made after this
-- process started are preserved.
local function updateIgnoreTests(layer, newLuaNames, newCNames,
      newGpuLuaNames, newGpuCNames)
   newGpuLuaNames = newGpuLuaNames or {}
   newGpuCNames   = newGpuCNames   or {}

   local srcBase    = configVals.source_dir .. "/" .. layer.src
   local luaIgnFile = srcBase .. "/luareg/ignore_lua_tests.lua"
   local cIgnFile   = srcBase .. "/creg/ignore_c_tests.lua"

   -- Adds entries from newList into existingList (no duplicates).
   -- Returns count of entries actually added.
   local function mergeList(existingList, newList)
      local set = {}
      for _, v in ipairs(existingList) do set[v] = true end
      local added = 0
      for _, v in ipairs(newList) do
         if not set[v] then
            set[v] = true
            table.insert(existingList, v)
            added = added + 1
         end
      end
      return added
   end

   -- Write { tests = {...}, gpu = {...} } to path with a header comment.
   local function writeIgnoreFile(path, header, tbl)
      local f = io.open(path, "w")
      f:write(header)
      f:write("return {\n")
      f:write("   tests = {\n")
      for _, v in ipairs(tbl.tests) do f:write(string.format("      %q,\n", v)) end
      f:write("   },\n")
      f:write("   gpu = {\n")
      for _, v in ipairs(tbl.gpu)   do f:write(string.format("      %q,\n", v)) end
      f:write("   },\n")
      f:write("}\n")
      f:close()
   end

   -- ---- Lua ignore file ----
   if #newLuaNames > 0 or #newGpuLuaNames > 0 then
      local existing = { tests = {}, gpu = {} }
      local gi = loadfile(luaIgnFile)
      if gi then
         local ok, loaded = pcall(gi)
         if ok and type(loaded) == "table" then existing = loaded end
      end
      existing.tests = existing.tests or {}
      existing.gpu   = existing.gpu   or {}

      local addedTests = mergeList(existing.tests, newLuaNames)
      local addedGpu   = mergeList(existing.gpu,   newGpuLuaNames)
      if addedTests > 0 or addedGpu > 0 then
         writeIgnoreFile(luaIgnFile,
            "-- Tests skipped by the Lua regression suite.\n"
            .. "-- Remove an entry manually to re-enable the test.\n"
            .. "-- gpu: Lua tests whose GPU variant timed out"
            .. " (CPU variant still runs).\n",
            existing)
         log(string.format(
            "[ignore] Updated %s (+%d timed-out, +%d GPU timed-out)\n",
            luaIgnFile, addedTests, addedGpu))
      end
   end

   -- ---- C ignore file ----
   if #newCNames > 0 or #newGpuCNames > 0 then
      local existing = { tests = {}, gpu = {} }
      local gi = loadfile(cIgnFile)
      if gi then
         local ok, loaded = pcall(gi)
         if ok and type(loaded) == "table" then existing = loaded end
      end
      existing.tests = existing.tests or {}
      existing.gpu   = existing.gpu   or {}

      local addedTests = mergeList(existing.tests, newCNames)
      local addedGpu   = mergeList(existing.gpu,   newGpuCNames)
      if addedTests > 0 or addedGpu > 0 then
         writeIgnoreFile(cIgnFile,
            "-- Tests skipped by the C regression suite.\n"
            .. "-- Remove an entry manually to re-enable the test.\n"
            .. "-- gpu: C tests whose GPU variant timed out"
            .. " (CPU variant still runs).\n",
            existing)
         log(string.format(
            "[ignore] Updated %s (+%d timed-out, +%d GPU timed-out)\n",
            cIgnFile, addedTests, addedGpu))
      end
   end
end

-- ---- Command action functions -----------------------------------------------

-- 'configure' command: set up the regression system for this machine.
local function config_action(args, name)
   isConfiguring = true

   local prefix = args.config_prefix
   if not prefix then
      -- Auto-detect from the installed config.mak.
      -- GKYL_EXEC_PATH is the bin dir (e.g. ~/gkylsoft/gkeyll/bin);
      -- config.mak lives one level up in share/.
      local gkeyllDir = GKYL_EXEC_PATH:match("^(.+)/bin$")
      if gkeyllDir then
         local mf = io.open(gkeyllDir .. "/share/config.mak", "r")
         if mf then
            for line in mf:lines() do
               local p = line:match("^PREFIX%s*=%s*(.+)%s*$")
               if p then prefix = p; break end
            end
            mf:close()
         end
      end
      prefix = prefix or (os.getenv("HOME") .. "/gkylsoft")
      log(string.format("Auto-detected prefix from config.mak: %s\n", prefix))
   end

   local mpiexec = args.config_mpiexec
      or (os.getenv("HOME") .. "/gkylsoft/openmpi/bin/mpiexec")
   local sourceDir = args.config_source_dir
   if not sourceDir then
      log("ERROR: --source-dir is required.\n"
         .. "Example: gkeyll runregression configure"
         .. " --source-dir /path/to/gkeyll\n")
      os.exit(1)
   end

   configure(prefix, mpiexec, sourceDir, args)
end

-- 'list' command: print all regression tests that would be run.
local function list_action(args, name)
   loadConfigure(args)
   local luaTests, cTests = list_tests(detectedLayer, args)
   if not args.c_only then
      for _, t in ipairs(luaTests) do print("[lua] " .. t.name) end
   end
   if not args.lua_only then
      for _, t in ipairs(cTests)   do print("[c]   " .. t.name) end
   end
end

-- 'run' command: execute regression tests and optionally create or check results.
-- On a GPU build (CC=nvcc), GPU-capable layers (vlasov, gyrokinetic, pkpm) run
-- each test twice: once in CPU mode, once in GPU mode.  Both are compared against
-- the same accepted baselines, and then CPU vs GPU output is compared to detect
-- GPU-specific divergence.  The 'create' action always forces CPU mode so that
-- accepted baselines are deterministic.
local function run_action(args, name)
   loadConfigure(args)

   local luaTests, cTests = list_tests(detectedLayer, args)
   local gpuTol = args.gpu_tol or 1e-7

   -- Per-test timeout in seconds (0 = unlimited).
   local timeoutSecs = args.timeout or 0
   if timeoutSecs > 0 and not TIMEOUT_CMD then
      log(string.format(
         "NOTE: 'timeout'/'gtimeout' not found; using Perl-based timeout.\n"))
   end

   -- Determine the post-run action. Without a sub-command the tests are run
   -- but results are neither saved nor compared (status = -1, "skip").
   local postRun = function(test, runDir, testType) return -1 end
   if args.create then
      postRun = create_action
   elseif args.check then
      postRun = check_action
   end

   -- Track timed-out tests per layer so we can update ignoretests.lua.
   local timedOutByLayer = {}
   local gpuTimedOutByLayer = {}
   for _, L in ipairs(LAYERS) do
      timedOutByLayer[L.name]    = { lua = {}, c = {} }
      gpuTimedOutByLayer[L.name] = { lua = {}, c = {} }
   end

   -- Helper: should this test get a GPU run?
   local function shouldDoGpu(test, testType)
      if not GPU_BUILD then return false end
      if not GPU_LAYERS[test.layer] then return false end
      if args.no_gpu then return false end
      -- Check GPU-specific ignore list.
      local ignKey = "gpu_" .. testType  -- "gpu_lua" or "gpu_c"
      local ignList = ignoreTests[test.layer] and ignoreTests[test.layer][ignKey] or {}
      local testBasename = (testType == "lua")
         and stripext(basename(test.file))
         or  stripext(basename(test.src))
      if lume.find(ignList, testBasename) then return false end
      return true
   end

   -- Helper: check whether any .gkyl output files exist in a directory.
   local function hasGkylOutput(dir)
      for fn in lfs.dir(dir) do
         if string.sub(fn, -5) == ".gkyl" then return true end
      end
      return false
   end

   -- Helper: run the GPU variant of a test after the CPU run.
   -- Saves CPU .gkyl files to _cpu_output/, clears .gkyl, runs GPU, compares.
   -- Returns: gpuStatus, gpuRuntime, cpuGpuDiff.
   local function runGpuVariant(test, runDir, testType, runFn, runArgs)
      local cpuOutputDir = runDir .. "/_cpu_output"
      mkdir(cpuOutputDir)
      os.execute(string.format("cp -f '%s'/*.gkyl '%s/' 2>/dev/null", runDir, cpuOutputDir))
      os.execute(string.format("rm -f '%s'/*.gkyl 2>/dev/null", runDir))

      -- Run the GPU variant.
      local gpuTm, gpuLog, _, gpuTimedOut = runFn(unpack(runArgs))
      local gpuStatus, gpuRuntime, cpuGpuDiff = -1, gpuTm, -1

      if gpuTimedOut then
         gpuStatus = -3
         log(string.format("... GPU variant TIMED OUT after %g sec\n", gpuTm))
      elseif not hasGkylOutput(runDir) then
         gpuStatus = -5  -- crash: no output produced
         log(string.format("... GPU variant CRASHED (no .gkyl output)\n"))
      else
         -- Compare GPU output vs accepted baselines.
         if args.check then
            gpuStatus = gpuCheck_vs_accepted(test, runDir, testType, gpuTol)
            if gpuStatus == 1 then
               log("... GPU vs accepted: passed.\n")
            else
               log(string.format("... GPU vs accepted: %s FAILED!\n", test.name))
            end
         else
            gpuStatus = -1  -- bare run or create, no check
         end

         -- Compare CPU output vs GPU output.
         cpuGpuDiff = compareCpuGpu(test, cpuOutputDir, runDir, testType, gpuTol)
         if cpuGpuDiff == 1 then
            log("... CPU vs GPU: match.\n")
         else
            log(string.format("... CPU vs GPU: %s DIFFER!\n", test.name))
         end
      end

      -- Clean up temporary CPU output directory.
      os.execute(string.format("rm -rf '%s'", cpuOutputDir))
      return gpuStatus, gpuRuntime, cpuGpuDiff
   end

   if GPU_BUILD then
      if args.create then
         log("Creating accepted results on CPU\n")
      elseif args.check then
         if args.no_gpu then
            log("Only CPU tests requested (--no-gpu)\n")
         else
            log(string.format("GPU build: GPU tolerance = %g\n", gpuTol))
         end
      end
   end

   -- Determine concurrency level.  0 = auto-detect physical core count.
   local jobCount = args.jobs or 1
   if jobCount == 0 then
      jobCount = physicalCpuCount()
      log(string.format("--jobs 0: auto-detected %d physical CPUs\n", jobCount))
   end
   if jobCount > 1 then
      log(string.format("Running regression tests in parallel (--jobs %d) ...\n\n",
         jobCount))
   else
      log("Running regression tests ...\n\n")
   end
   local tmStart = Time.clock()

   if jobCount <= 1 then
      -- ================================================================
      -- SERIAL PATH (default): one test at a time, preserving all
      -- existing behaviour exactly.
      -- ================================================================

      -- ---- Lua tests --------------------------------------------------------
      if not args.c_only then
         for _, test in ipairs(luaTests) do
            layerCounts[test.layer].total = layerCounts[test.layer].total + 1
            local doGpu = shouldDoGpu(test, "lua")

            -- On a GPU build for a GPU-capable layer, force CPU mode via -G.
            -- For 'create', always force CPU to produce deterministic baselines.
            local cpuMode = (GPU_BUILD and GPU_LAYERS[test.layer]) and "cpu" or nil

            local runtm, runlog, runDir, timedOut = runLuaTest(test, timeoutSecs, cpuMode)

            if timedOut then
               table.insert(timedOutByLayer[test.layer].lua, stripext(basename(test.file)))
               insertRegressionData(
                  test.layer, runID, test.name, "lua", -3, runtm, "TIMED OUT")
               layerCounts[test.layer].failed = layerCounts[test.layer].failed + 1
            else
               local status, checkLog = postRun(test, runDir, "lua")
               checkLog = checkLog or ""

               -- GPU variant: only on check or bare run, not on create.
               local gpuStatus, gpuRuntime, cpuGpuDiff = -1, 0, -1
               if doGpu and not args.create then
                  layerCounts[test.layer].gpu_total = layerCounts[test.layer].gpu_total + 1

                  gpuStatus, gpuRuntime, cpuGpuDiff = runGpuVariant(
                     test, runDir, "lua",
                     runLuaTest, {test, timeoutSecs, "gpu"})

                  if gpuStatus == -3 then
                     table.insert(gpuTimedOutByLayer[test.layer].lua,
                        stripext(basename(test.file)))
                     layerCounts[test.layer].gpu_failed = layerCounts[test.layer].gpu_failed + 1
                  elseif gpuStatus == -5 or gpuStatus == 0 then
                     layerCounts[test.layer].gpu_failed = layerCounts[test.layer].gpu_failed + 1
                  elseif gpuStatus == 1 then
                     layerCounts[test.layer].gpu_passed = layerCounts[test.layer].gpu_passed + 1
                  end
               end

               insertRegressionData(
                  test.layer, runID, test.name, "lua",
                  status, runtm,
                  runlog .. (checkLog ~= "" and "\n" .. checkLog or ""),
                  gpuStatus, gpuRuntime, cpuGpuDiff)
            end
         end
      end

      -- ---- C tests ----------------------------------------------------------
      if not args.lua_only then
         for _, test in ipairs(cTests) do
            layerCounts[test.layer].total = layerCounts[test.layer].total + 1
            local doGpu = shouldDoGpu(test, "c")

            -- CPU run: compile + run (no -g flag = CPU mode by default).
            -- When doGpu is true, keep the binary for the subsequent GPU re-run.
            local runtm, runlog, runDir, timedOut, compileFailed =
               runCTest(test, timeoutSecs, nil, false, doGpu)

            if compileFailed then
               insertRegressionData(
                  test.layer, runID, test.name, "c", -4, runtm, runlog)
               layerCounts[test.layer].failed = layerCounts[test.layer].failed + 1
            elseif timedOut then
               table.insert(timedOutByLayer[test.layer].c,
                  stripext(basename(test.src)))
               insertRegressionData(
                  test.layer, runID, test.name, "c", -3, runtm, "TIMED OUT")
               layerCounts[test.layer].failed = layerCounts[test.layer].failed + 1
            else
               local status, checkLog = postRun(test, runDir, "c")
               checkLog = checkLog or ""

               -- GPU variant: only on check or bare run, not on create.
               local gpuStatus, gpuRuntime, cpuGpuDiff = -1, 0, -1
               if doGpu and not args.create then
                  layerCounts[test.layer].gpu_total = layerCounts[test.layer].gpu_total + 1

                  -- Re-run the same binary with -g (skipCompile=true).
                  gpuStatus, gpuRuntime, cpuGpuDiff = runGpuVariant(
                     test, runDir, "c",
                     runCTest, {test, timeoutSecs, "gpu", true})

                  if gpuStatus == -3 then
                     table.insert(gpuTimedOutByLayer[test.layer].c,
                        stripext(basename(test.src)))
                     layerCounts[test.layer].gpu_failed = layerCounts[test.layer].gpu_failed + 1
                  elseif gpuStatus == -5 or gpuStatus == 0 then
                     layerCounts[test.layer].gpu_failed = layerCounts[test.layer].gpu_failed + 1
                  elseif gpuStatus == 1 then
                     layerCounts[test.layer].gpu_passed = layerCounts[test.layer].gpu_passed + 1
                  end
               end

               -- Clean up the compiled C binary now that both CPU and GPU runs are done.
               if doGpu then
                  local testname = stripext(basename(test.src))
                  os.execute(string.format("rm -f '%s/%s' '%s/%s.d' 2>/dev/null",
                     runDir, testname, runDir, testname))
                  os.execute(string.format("rm -rf '%s/%s.dSYM' 2>/dev/null",
                     runDir, testname))
               end

               insertRegressionData(
                  test.layer, runID, test.name, "c",
                  status, runtm,
                  runlog .. (checkLog ~= "" and "\n" .. checkLog or ""),
                  gpuStatus, gpuRuntime, cpuGpuDiff)
            end
         end
      end -- if not args.lua_only

   else
      -- ================================================================
      -- PARALLEL PATH: up to jobCount tests run concurrently per batch.
      --
      -- Phase 1 (C only): compile all C tests serially — fast, and avoids
      --   Makefile conflicts.  Compile failures are recorded immediately
      --   and excluded from the run batch.
      -- Phase 2a (Lua): build prep list; handle mpiSkip tests inline.
      --   Run in batches of jobCount via executeBatch.
      -- Phase 2b (C): run compiled-OK tests in batches via executeBatch.
      -- Collect phase (after each executeBatch): log results in original
      --   order, call postRun, run GPU variants (always serial), record DB.
      --
      -- GPU variants stay serial to avoid GPU memory contention.
      -- ================================================================

      -- Helper: collect results for a single Lua test prep + batch result.
      local function collectLua(prep, r)
         local test = prep.test
         if r.timedOut then
            log(string.format("\n[Lua] %s TIMED OUT (%g sec)\n", test.name, r.runtm))
            table.insert(timedOutByLayer[test.layer].lua, stripext(basename(test.file)))
            insertRegressionData(
               test.layer, runID, test.name, "lua", -3, r.runtm, "TIMED OUT")
            layerCounts[test.layer].failed = layerCounts[test.layer].failed + 1
         else
            log(string.format("\n[Lua] %s completed (%g sec)\n", test.name, r.runtm))
            verboseLog(r.runlog)
            local status, checkLog = postRun(test, prep.runDir, "lua")
            checkLog = checkLog or ""
            local gpuStatus, gpuRuntime, cpuGpuDiff = -1, 0, -1
            if prep.doGpu and not args.create then
               layerCounts[test.layer].gpu_total = layerCounts[test.layer].gpu_total + 1
               gpuStatus, gpuRuntime, cpuGpuDiff = runGpuVariant(
                  test, prep.runDir, "lua",
                  runLuaTest, {test, timeoutSecs, "gpu"})
               if gpuStatus == -3 then
                  table.insert(gpuTimedOutByLayer[test.layer].lua,
                     stripext(basename(test.file)))
                  layerCounts[test.layer].gpu_failed = layerCounts[test.layer].gpu_failed + 1
               elseif gpuStatus == -5 or gpuStatus == 0 then
                  layerCounts[test.layer].gpu_failed = layerCounts[test.layer].gpu_failed + 1
               elseif gpuStatus == 1 then
                  layerCounts[test.layer].gpu_passed = layerCounts[test.layer].gpu_passed + 1
               end
            end
            insertRegressionData(
               test.layer, runID, test.name, "lua",
               status, r.runtm,
               r.runlog .. (checkLog ~= "" and "\n" .. checkLog or ""),
               gpuStatus, gpuRuntime, cpuGpuDiff)
         end
      end

      -- Helper: collect results for a single C test prep + batch result.
      local function collectC(prep, r)
         local test     = prep.test
         local testname = stripext(basename(test.src))
         local runDir   = prep.runDir
         if r.timedOut then
            log(string.format("\n[C] %s TIMED OUT (%g sec)\n", test.name, r.runtm))
            table.insert(timedOutByLayer[test.layer].c, testname)
            insertRegressionData(
               test.layer, runID, test.name, "c", -3, r.runtm, "TIMED OUT")
            layerCounts[test.layer].failed = layerCounts[test.layer].failed + 1
         else
            log(string.format("\n[C] %s completed (%g sec)\n", test.name, r.runtm))
            verboseLog(r.runlog)
            local status, checkLog = postRun(test, runDir, "c")
            checkLog = checkLog or ""
            local gpuStatus, gpuRuntime, cpuGpuDiff = -1, 0, -1
            if prep.doGpu and not args.create then
               layerCounts[test.layer].gpu_total = layerCounts[test.layer].gpu_total + 1
               -- Re-run existing binary with -g (skipCompile=true).
               gpuStatus, gpuRuntime, cpuGpuDiff = runGpuVariant(
                  test, runDir, "c",
                  runCTest, {test, timeoutSecs, "gpu", true})
               if gpuStatus == -3 then
                  table.insert(gpuTimedOutByLayer[test.layer].c, testname)
                  layerCounts[test.layer].gpu_failed = layerCounts[test.layer].gpu_failed + 1
               elseif gpuStatus == -5 or gpuStatus == 0 then
                  layerCounts[test.layer].gpu_failed = layerCounts[test.layer].gpu_failed + 1
               elseif gpuStatus == 1 then
                  layerCounts[test.layer].gpu_passed = layerCounts[test.layer].gpu_passed + 1
               end
            end
            local fullLog = (prep.compileLog or "") .. "\n" .. r.runlog
               .. (checkLog ~= "" and "\n" .. checkLog or "")
            insertRegressionData(
               test.layer, runID, test.name, "c",
               status, r.runtm, fullLog, gpuStatus, gpuRuntime, cpuGpuDiff)
         end
         -- Clean up binary after both CPU and (optional) GPU runs.
         os.execute(string.format("rm -f '%s/%s' '%s/%s.d' 2>/dev/null",
            runDir, testname, runDir, testname))
         os.execute(string.format("rm -rf '%s/%s.dSYM' 2>/dev/null",
            runDir, testname))
      end

      -- Phase 1: compile all C tests serially.
      local cPreps = {}
      if not args.lua_only then
         for _, test in ipairs(cTests) do
            layerCounts[test.layer].total = layerCounts[test.layer].total + 1
            local doGpu = shouldDoGpu(test, "c")
            local prep  = prepareCRun(test, timeoutSecs, nil, false)
            prep.doGpu  = doGpu
            table.insert(cPreps, prep)
            if prep.compileFailed then
               insertRegressionData(
                  test.layer, runID, test.name, "c", -4,
                  prep.compileSecs, "COMPILE FAILED:\n" .. prep.compileLog)
               layerCounts[test.layer].failed = layerCounts[test.layer].failed + 1
            end
         end
      end

      -- Phase 2a: prepare Lua preps and handle mpiSkip tests inline.
      local luaPreps = {}
      if not args.c_only then
         for _, test in ipairs(luaTests) do
            layerCounts[test.layer].total = layerCounts[test.layer].total + 1
            local doGpu   = shouldDoGpu(test, "lua")
            local cpuMode = (GPU_BUILD and GPU_LAYERS[test.layer]) and "cpu" or nil
            local prep    = prepareLuaRun(test, timeoutSecs, cpuMode)
            prep.doGpu    = doGpu
            if prep.mpiSkip then
               log(string.format("**** NOT RUNNING PARALLEL TEST %s\n", test.name))
               insertRegressionData(test.layer, runID, test.name, "lua", -1, 0, "")
            else
               table.insert(luaPreps, prep)
            end
         end
      end

      -- Phase 2a continued: run Lua tests in batches.
      for bStart = 1, #luaPreps, jobCount do
         local batch = {}
         for i = bStart, math.min(bStart + jobCount - 1, #luaPreps) do
            table.insert(batch, luaPreps[i])
         end
         log(string.format("\n[Batch-Lua] Launching %d test(s) ...\n", #batch))
         local bResults = executeBatch(batch)
         for i, prep in ipairs(batch) do
            collectLua(prep, bResults[i])
         end
      end

      -- Phase 2b: run compiled-OK C tests in batches.
      local cRunPreps = {}
      for _, prep in ipairs(cPreps) do
         if not prep.compileFailed then
            table.insert(cRunPreps, prep)
         end
      end
      for bStart = 1, #cRunPreps, jobCount do
         local batch = {}
         for i = bStart, math.min(bStart + jobCount - 1, #cRunPreps) do
            table.insert(batch, cRunPreps[i])
         end
         log(string.format("\n[Batch-C] Launching %d test(s) ...\n", #batch))
         local bResults = executeBatch(batch)
         for i, prep in ipairs(batch) do
            collectC(prep, bResults[i])
         end
      end

   end -- if jobCount <= 1 / else

   log(string.format(
      "\nAll regression tests completed in %g secs\n", Time.clock() - tmStart))

   -- Update ignoretests.lua for any layers that had timed-out tests.
   if timeoutSecs > 0 then
      for _, layer in ipairs(LAYERS) do
         local to    = timedOutByLayer[layer.name]
         local gpuTo = gpuTimedOutByLayer[layer.name]
         if #to.lua > 0 or #to.c > 0
            or #gpuTo.lua > 0 or #gpuTo.c > 0 then
            updateIgnoreTests(layer, to.lua, to.c, gpuTo.lua, gpuTo.c)
         end
      end
   end

   -- Finalize: write per-layer summary rows into each layer's database.
   -- We only write a metadata row for layers that actually had tests.
   for _, layer in ipairs(LAYERS) do
      local cnt = layerCounts[layer.name]
      if cnt.total > 0 then
         insertRegressionMeta(layer.name, runID, runDate,
            cnt.total, cnt.passed, cnt.failed,
            GPU_BUILD and 1 or 0, cnt.gpu_passed, cnt.gpu_failed)
         local summary = string.format(
            "  Layer %-12s  total=%d  passed=%d  failed=%d",
            layer.name, cnt.total, cnt.passed, cnt.failed)
         if GPU_BUILD and GPU_LAYERS[layer.name] then
            summary = summary .. string.format(
               "  gpu_total=%d  gpu_pass=%d  gpu_fail=%d",
               cnt.gpu_total, cnt.gpu_passed, cnt.gpu_failed)
         end
         log(summary .. "\n")
      end
   end
end

-- 'listunit' command.
local function listunit_action(args, name)
   loadConfigure(args)
   local lua, cxx = list_unit_tests(args)
   lume.each(lua, print)
   lume.each(cxx, print)
end

-- 'rununit' command.
local function rununit_action(args, name)
   loadConfigure(args)
   local lua, cxx = list_unit_tests(args)
   log("Running unit tests ...\n\n")
   local tmStart = Time.clock()
   lume.each(lua, runLuaUnitTest)
   lume.each(cxx, runCxxUnitTest)
   log(string.format("All unit tests completed in %g secs\n", Time.clock() - tmStart))
end

-- ---- CLI parser -------------------------------------------------------------

local parser = argparse()
   :name("runregression")
   :require_command(true)
   :help_description_margin(30)
   :help_vertical_space(1)
   :description [[
Run Gkeyll regression tests across the hierarchical layer architecture
(moments -> vlasov -> gyrokinetic -> pkpm).

Each layer has Lua regression tests (luareg/) and C regression tests (creg/).
Results are stored in per-layer SQLite databases under gkeyll-results/.

Typical workflow:
  1. Build and install: make install -j4
  2. Configure: gkeyll runregression configure --source-dir /absolute/path/to/gkeyll/
  3. Create baselines: gkeyll runregression run create --timeout 120
  4. Check results:    gkeyll runregression run check  --timeout 120
  5. Layer-specific:   gkeyll runregression run moments check

C regression tests are compiled on-the-fly using the installed
share/Makefile (PREFIX/gkeyll/share/Makefile). No separate 'make regression'
step is needed.
]]

parser:flag("-v --verbose", "Print verbose messages as tests are run")

-- 'configure' command ---------------------------------------------------------
-- Sets up directories, databases, and writes ~/runregression.config.lua.
local c_conf = parser:command("configure", "Configure regression tests")
   :action(config_action)

c_conf:option("-s --source-dir",
   "REQUIRED. Absolute path to gkeyll source root\n"
   .. "(parent of moments/, vlasov/, gyrokinetic/, pkpm/).\n"
   .. "Example: --source-dir /path/to/gkeyll")
   :target("config_source_dir")
c_conf:option("-p --prefix",
   "Where to write gkeyll-results/.\n"
   .. "Auto-detected from config.mak if omitted.")
   :target("config_prefix")
c_conf:option("-m --mpiexec",
   "Full path to MPI executable.\n"
   .. "(MPI parallelism not yet implemented).")
   :target("config_mpiexec")
c_conf:flag("--drop-tables",
   "Drop and re-create all SQL tables\n"
   .. "(erases existing regression data).", false)

-- 'list' command --------------------------------------------------------------
-- Lists all regression tests that would be run (useful for inspection).
local c_list = parser:command("list", "List all regression tests")
   :action(list_action)
c_list:option("-r --run-only",
   "List only this test. Accepts a bare name (rt_foo), an absolute path, or a directory.")
c_list:flag("-a --all",
   "List all tests, bypassing ignore_lua_tests.lua / ignore_c_tests.lua.\n"
   .. "Combine with --lua-only or --c-only to bypass only that suite's ignore list.")
c_list:flag("-m --moat", "Only list MOAT (Mother Of All Tests) regression tests\n"
   .. "A condensed suite of the most comprehensive regression tests.")
c_list:flag("-c --c-only",   "Only list C regression tests (skip Lua tests)")
c_list:flag("-l --lua-only", "Only list Lua regression tests (skip C tests)")

-- 'run' command ---------------------------------------------------------------
-- Runs regression tests. Layer is extracted from the command line BEFORE
-- argparse (see the pre-processing block near the top of this file), so:
--   gkeyll runregression run check          -> all layers, check
--   gkeyll runregression run moments check  -> moments layer only, check
local c_run = parser:command("run",
   "Run regression tests (check or create).\n"
   .. "Prefix with a layer name to restrict: 'run moments check', 'run vlasov create'.\n"
   .. "Without check/create, tests run but results are not saved or compared.")
   :require_command(false)
   :action(run_action)
c_run:option("-r --run-only",
   "Run only these tests (comma-separated).\n"
   .. "Accepts bare names (rt_foo), absolute paths, or directories.\n"
   .. "Use --c-only/--lua-only to disambiguate when a C and Lua test share the same name.")
c_run:flag("-a --all",
   "Run all tests, bypassing ignore_lua_tests.lua / ignore_c_tests.lua.\n"
   .. "Combine with --lua-only or --c-only to bypass only that suite's ignore list.")
c_run:flag("-m --moat", "Only run MOAT (Mother Of All Tests) regression tests\n"
   .. "A condensed suite of the most comprehensive regression tests.")
c_run:flag("-c --c-only",   "Only run C regression tests (skip Lua tests)")
c_run:flag("-l --lua-only", "Only run Lua regression tests (skip C tests)")
c_run:option("-t --timeout",
   "Per-test timeout in seconds (0 = unlimited).\n"
   .. "Timed-out tests are added to ignoretests.lua automatically.")
   :convert(tonumber)
   :default(0)
c_run:option("--gpu-tol",
   "Tolerance for GPU-vs-accepted and CPU-vs-GPU comparisons (default 1e-7).\n"
   .. "Only used on GPU builds (CC=nvcc).")
   :convert(tonumber)
   :default(1e-7)
c_run:flag("--no-gpu",
   "Skip GPU testing even on a GPU build.")
c_run:option("-j --jobs",
   "Concurrent tests per batch (0 = physical core count, 1 = serial).\n"
   .. "C compilation is always serial; GPU variants always run serially.")
   :convert(tonumber)
   :default(1)

c_run:command("check",
   "Compare test output against accepted baselines;\n"
   .. "report pass/fail and write to the SQLite DB.")
c_run:command("create",
   "Run tests and save output as accepted baselines.\n"
   .. "On GPU builds, create always runs in CPU mode so baselines are deterministic.")

-- 'listunit' command ----------------------------------------------------------
parser:command("listunit", "List all unit tests")
   :action(listunit_action)

-- 'rununit' command -----------------------------------------------------------
local c_rununit = parser:command("rununit", "Run unit tests.")
   :require_command(false)
   :action(rununit_action)

-- ---- Parse and dispatch -----------------------------------------------------
-- Use filteredCmds (layer names removed) rather than GKYL_COMMANDS_L.
local _ = parser:parse(filteredCmds)
