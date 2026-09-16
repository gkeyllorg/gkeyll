-- Gkyl ------------------------------------------------------------------------
--
-- Query the per-layer regression databases created by runregression.
-- Updated for the hierarchical layer architecture: each layer (moments,
-- vlasov, gyrokinetic, pkpm) now has its own SQLite database stored at
--   gkeyll-results/<layer>/regressiondb
-- The --layer option selects which layer's database to query.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Database schema (one file per layer, created by runregression configure):
--
-- table RegressionMeta (
--   guid               text,
--   tstamp             text,
--   GKYL_EXEC          text,
--   GKYL_GIT_CHANGESET text,
--   GKYL_BUILD_DATE    text,
--   ntotal             integer,
--   npass              integer,
--   nfail              integer,
--   gpu_build          integer,   -- 1=GPU build, 0=CPU-only
--   ngpu_pass          integer,
--   ngpu_fail          integer
-- );
--
-- table RegressionData (
--   guid         text,
--   name         text,
--   test_type    text,   -- 'lua' or 'c'
--   status       integer,
--   runtime      real,
--   runlog       text,
--   gpu_status   integer,  -- 1=pass, 0=fail, -1=skip, -3=timeout, -5=crash
--   (status uses the same -5=crash for a CPU leg that exited non-zero)
--   gpu_runtime  real,
--   cpu_gpu_diff integer   -- 1=match, 0=differ, -1=n/a
-- );

local argparse = require "Lib.argparse"
local Logger   = require "Lib.Logger"
local sql      = require "sqlite3"

-- SQLite connection handle (opened by configure()).
local sqlConn = nil

-- Human-readable status strings for the integer status codes stored in the DB.
local statusToString = {
   [-5] = "crash", [-4] = "compile_fail", [-3] = "timeout",
   [-2] = "create", [-1] = "skip", [0] = "fail", [1] = "pass",
}
-- Human-readable CPU-vs-GPU diff strings.
local cpuGpuDiffToString = { [-1] = "n/a", [0] = "differ", [1] = "match" }

-- GKYL_OUT_PREFIX is required by some Lua modules at load time.
GKYL_OUT_PREFIX = lfs.currentdir() .. "/" .. "queryrdb"

local log = Logger { logToFile = true }

-- Path of the configuration file written by 'runregression configure'.
local confFile = os.getenv("HOME") .. "/runregression.config.lua"

-- ---- Database connection ----------------------------------------------------
-- Opens the correct database based on the parsed args:
--   --db <path>    → open the specified file directly (for ad-hoc inspection)
--   --layer <name> → derive path from the config file:
--                    results_dir/<layer>/regressiondb
-- Having a --layer option (rather than requiring the user to know the database
-- path) keeps the interface consistent with how runregression names databases.
local function configure(args)
   if args.db then
      -- Direct path override: useful for inspecting an arbitrary DB file.
      sqlConn = sql.open(args.db)
   else
      if not args.layer then
         print("ERROR: --layer <name> is required (or use --db <path> for a specific file).")
         print("Valid layers: moments, vlasov, gyrokinetic, pkpm")
         os.exit(1)
      end

      local f = loadfile(confFile)
      if not f then
         print("Regression tests not configured. Run 'runregression configure' first.")
         os.exit(1)
      end
      local configVals = f()

      local dbPath = string.format("%s/%s/regressiondb",
         configVals.results_dir, args.layer)
      if not lfs.attributes(dbPath) then
         print(string.format(
            "Database for layer '%s' not found at:\n  %s\n"
            .. "Run 'runregression configure' and then 'runregression run create' first.",
            args.layer, dbPath))
         os.exit(1)
      end
      sqlConn = sql.open(dbPath)
   end
end

-- ---- Database read helpers --------------------------------------------------

-- Returns all rows from RegressionMeta as a Lua list of tables.
local function read_metatable()
   local t, nrow = sqlConn:exec("select * from RegressionMeta")
   local dbMeta = {}
   for i = 1, nrow do
      dbMeta[i] = {
         guid      = t['guid'][i],
         tstamp    = t['tstamp'][i],
         changeset = t['GKYL_GIT_CHANGESET'][i],
         builddate = t['GKYL_BUILD_DATE'][i],
         ntotal    = t['ntotal'][i],
         npass     = t['npass'][i],
         nfail     = t['nfail'][i],
         gpu_build = t['gpu_build'] and t['gpu_build'][i] or 0,
         ngpu_pass = t['ngpu_pass'] and t['ngpu_pass'][i] or 0,
         ngpu_fail = t['ngpu_fail'] and t['ngpu_fail'][i] or 0,
      }
   end
   return dbMeta
end

-- Returns all RegressionData rows for a specific guid.
-- Also returns the maximum name length (for column-aligned output).
local function read_tests_with_id(guid)
   local t, nrow = sqlConn:exec(
      string.format("select * from RegressionData where guid=='%s'", guid))
   local maxNm = 0
   local dbData = {}
   for i = 1, nrow do
      local gpuStat = t['gpu_status'] and tonumber(t['gpu_status'][i]) or -1
      local cgDiff  = t['cpu_gpu_diff'] and tonumber(t['cpu_gpu_diff'][i]) or -1
      dbData[i] = {
         name         = t['name'][i],
         test_type    = t['test_type'] and t['test_type'][i] or "?",
         status       = statusToString[tonumber(t['status'][i])],
         runtime      = t['runtime'][i],
         runlog       = t['runlog'][i],
         gpu_status   = statusToString[gpuStat] or "skip",
         gpu_runtime  = t['gpu_runtime'] and t['gpu_runtime'][i] or 0,
         cpu_gpu_diff = cpuGpuDiffToString[cgDiff] or "n/a",
      }
      maxNm = math.max(maxNm, string.len(t['name'][i]))
   end
   return dbData, maxNm
end

-- ---- Command actions --------------------------------------------------------

-- 'summary': prints a one-line summary for every run stored in the database.
-- Rows are displayed in reverse insertion order (most recent run = ID 1).
local function summary_action(args, name)
   configure(args)

   local dbMeta = read_metatable()
   local nrow   = #dbMeta

   -- Check if any run was a GPU build; if so, show GPU columns.
   local hasGpu = false
   for _, d in ipairs(dbMeta) do
      if tonumber(d.gpu_build) == 1 then hasGpu = true; break end
   end

   if hasGpu then
      local fmt = "%-4s: %-20s %-30s %-5s %-5s %-5s %-4s %-7s %-7s"
      print(string.format(fmt, "ID", "Time-Stamp", "Changeset",
         "Total", "Pass", "Fail", "GPU", "GPUPass", "GPUFail"))
      for i, d in ipairs(dbMeta) do
         print(string.format(fmt,
            nrow - i + 1, d.tstamp, d.changeset,
            tonumber(d.ntotal), tonumber(d.npass), tonumber(d.nfail),
            tonumber(d.gpu_build) == 1 and "yes" or "no",
            tonumber(d.ngpu_pass), tonumber(d.ngpu_fail)
         ))
      end
   else
      local fmt = "%-4s: %-20s %-30s %-5s %-5s %-5s"
      print(string.format(fmt, "ID", "Time-Stamp", "Changeset", "Total", "Pass", "Fail"))
      for i, d in ipairs(dbMeta) do
         print(string.format(fmt,
            nrow - i + 1, d.tstamp, d.changeset,
            tonumber(d.ntotal), tonumber(d.npass), tonumber(d.nfail)
         ))
      end
   end
end

-- 'query': prints the test results for a specific run (identified by ID from
-- the summary command). Supports filtering by pass/fail and comma-list output.
local function query_action(args, name)
   configure(args)

   local dbMeta = read_metatable()
   -- The summary command displays most-recent first (ID 1 = last row in DB),
   -- so invert the index here to select the correct DB row.
   local idx = #dbMeta - args.id + 1

   if args.net_time then
      -- Sum wall-clock runtimes for all tests in this run.
      local t, nrow = sqlConn:exec(string.format(
         "select runtime from RegressionData where guid=='%s'",
         dbMeta[idx].guid))
      local nettm = 0.0
      for i = 1, nrow do nettm = nettm + t['runtime'][i] end
      print(string.format("%.4g", nettm))
      return
   end

   local dbData, maxNm = read_tests_with_id(dbMeta[idx].guid)

   -- Filtering predicate: apply --fail-only, --pass-only, or --gpu-fail-only.
   local function shouldShow(d)
      if args.gpu_fail_only then
         -- Show tests where GPU failed but CPU passed (most interesting diagnostic).
         return d.status == "pass" and (d.gpu_status == "fail" or d.gpu_status == "crash")
      end
      if args.fail_only then return d.status == "fail" end
      if args.pass_only then return d.status == "pass" end
      return true
   end

   if args.test and args.test ~= "" then
      -- --test accepts either a row number (integer) or a name / substring.
      local matches = {}
      local asNum = tonumber(args.test)
      if asNum and asNum > 0 then
         -- Numeric: index directly into the result list.
         local tidx = math.min(asNum, #dbData)
         matches = { dbData[tidx] }
      else
         -- Name / substring: collect all tests whose stored name contains the value.
         for _, d in ipairs(dbData) do
            if string.find(d.name, args.test, 1, true) then
               table.insert(matches, d)
            end
         end
         if #matches == 0 then
            print(string.format("No test matching '%s' found in run %d.", args.test, args.id))
         end
      end
      for _, d in ipairs(matches) do
         print("=== ")
         print(string.format("=== Log for test %s [%s]", d.name, d.test_type))
         print("=== ")
         print(d.runlog)
      end
   elseif args.comma_list then
      -- Comma-separated list of matching test names.
      for _, d in pairs(dbData) do
         if shouldShow(d) then io.write(d.name .. ',') end
      end
      io.write('\n')
   else
      -- Tabular output with test type column and GPU columns.
      local nm  = maxNm + 2
      local fmt  = "%-4s: %-5s %-" .. nm .. "s %-7s %-9s %-10s %-9s %-7s"
      local fmt1 = "%-4s: %-5s %-" .. nm .. "s %-7s %.4g      %-10s %.4g      %-7s"
      print(string.format(fmt, "ID", "Type", "Name", "Status", "Run-Time",
         "GPU-Status", "GPU-Time", "CPU=GPU"))
      for i, d in pairs(dbData) do
         if shouldShow(d) then
            print(string.format(fmt1,
               i, d.test_type, d.name, d.status, d.runtime,
               d.gpu_status, d.gpu_runtime, d.cpu_gpu_diff))
         end
      end
   end
end

-- 'delete': removes a run's rows from both tables (identified by summary ID).
local function delete_action(args, name)
   configure(args)
   if tonumber(args.id) < 1 then
      print("No deletions done. Specify a run ID as shown by the summary command.")
      return
   end

   local dbMeta = read_metatable()
   local idx  = #dbMeta - args.id + 1
   local guid = dbMeta[idx].guid

   sqlConn:exec(string.format(
      "delete from RegressionData where guid=='%s'", guid))
   sqlConn:exec(string.format(
      "delete from RegressionMeta where guid=='%s'", guid))
   print(string.format("Deleted run %d (guid=%s).", args.id, guid))
end

-- 'history': shows the history of a single named test across all stored runs.
-- Useful for spotting trends in run time or pass/fail transitions.
local function history_action(args, name)
   configure(args)

   if args.regression == nil then return end

   -- Strip a leading "./" from the test name if present (for compatibility
   -- with the old system that stored names relative to the regression dir).
   local function trimname(nm)
      local s, e = string.find(nm, "./")
      if s and s == 1 then return string.sub(nm, e + 1) end
      return nm
   end

   local tNm = trimname(args.regression)

   if args.time_only then
      -- Print only runtimes for runs where the test passed.
      local dat, nrow = sqlConn:exec(string.format(
         "select * from RegressionData where name=='%s' and status==1", tNm))
      local fmt = "## %-20s %-4s"
      print(string.format(fmt, "Time-Stamp", "Run-Time"))
      for i = 1, nrow do
         local guid = dat['guid'][i]
         local tstamp, changeset = sqlConn:rowexec(string.format(
            "select tstamp, GKYL_GIT_CHANGESET from RegressionMeta where guid='%s'", guid))
         print(string.format("%.4g", dat['runtime'][i]))
      end
   else
      local dat, nrow = sqlConn:exec(string.format(
         "select * from RegressionData where name=='%s'", tNm))
      local fmt  = "%-20s %-30s %-5s %-7s %-9s %-10s %-7s"
      local fmt1 = "%-20s %-30s %-5s %-7s %.4g      %-10s %-7s"
      print(string.format(fmt, "Time-Stamp", "Changeset", "Type",
         "Status", "Run-Time", "GPU-Status", "CPU=GPU"))
      for i = 1, nrow do
         local guid = dat['guid'][i]
         local tstamp, changeset = sqlConn:rowexec(string.format(
            "select tstamp, GKYL_GIT_CHANGESET from RegressionMeta where guid='%s'", guid))
         local stat    = statusToString[tonumber(dat['status'][i])]
         local ttype   = dat['test_type'] and dat['test_type'][i] or "?"
         local gpuStat = dat['gpu_status']
            and (statusToString[tonumber(dat['gpu_status'][i])] or "skip")
            or "skip"
         local cgDiff  = dat['cpu_gpu_diff']
            and (cpuGpuDiffToString[tonumber(dat['cpu_gpu_diff'][i])] or "n/a")
            or "n/a"
         print(string.format(fmt1, tstamp, changeset, ttype,
            stat, dat['runtime'][i], gpuStat, cgDiff))
      end
   end
end

-- ---- CLI parser -------------------------------------------------------------

local parser = argparse()
   :name("queryrdb")
   :require_command(true)
   :description [[
Query the per-layer regression databases created by the runregression tool.

Each layer (moments, vlasov, gyrokinetic, pkpm) has its own SQLite database.
Use --layer to select which layer's database to query, or --db to open an
arbitrary database file directly.

Examples:
  gkeyll queryrdb --layer moments summary
  gkeyll queryrdb --layer vlasov query --id 1 --fail-only
  gkeyll queryrdb --layer gyrokinetic history -r moments/luareg/rt_gk_sheath_1x2v_p1.lua
]]

-- Global options (apply to all sub-commands).
parser:option("--layer",
   "Layer whose database to query (moments/vlasov/gyrokinetic/pkpm)")
parser:option("--db",
   "Open this specific database file instead of using --layer")

-- 'summary' command.
parser:command("summary", "Print a one-line summary of every stored run")
   :action(summary_action)

-- 'query' command.
local c_query = parser:command("query",
   "Print detailed results for a specific run (identified by summary ID)")
   :action(query_action)
c_query:option("-i --id", "ID of the run to query (from summary command)", 1)
c_query:flag("-f --fail-only", "Show only failed tests", false)
c_query:flag("-p --pass-only", "Show only passed tests", false)
c_query:flag("-g --gpu-fail-only",
   "Show only tests where GPU failed but CPU passed", false)
c_query:flag("-l --comma-list", "Output test names as a comma-separated list", false)
c_query:option("-t --test",
   "Print the full run log for a test. Accepts a row number (from the ID column)\n"
   .. "or a name / substring (e.g. rt_gk_sheath_2x2v_p1).", "")
c_query:flag("--net-time", "Print total wall-clock time for the run", false)

-- 'delete' command.
local c_delete = parser:command("delete",
   "Delete a run's data from the database")
   :action(delete_action)
c_delete:option("-i --id",
   "ID of the run to delete (from summary command)", 0)

-- 'history' command.
local c_history = parser:command("history",
   "Show the history of a single named test across all stored runs")
   :action(history_action)
c_history:option("-r --regression", "Name of the test to query history for")
c_history:flag("--time-only", "Print only run times (for passed runs)", false)

-- ---- Parse and dispatch -----------------------------------------------------
local _ = parser:parse(GKYL_COMMANDS_L)
