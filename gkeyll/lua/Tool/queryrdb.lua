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
--   run_mode            text
-- );
--
-- table RegressionData (
--   guid         text,
--   name         text,
--   test_type    text,   -- 'lua' or 'c'
--   status       integer,
--   runtime      real,
--   runlog       text,
-- );

local argparse = require "Lib.argparse"
local Logger   = require "Lib.Logger"
local sql      = require "sqlite3"

-- SQLite connection handle (opened by configure()).
local sqlConn = nil

-- Human-readable status strings for the integer status codes stored in the DB.
local statusToString = {
   [-6] = "crash", [-5] = "no_output", [-4] = "compile_fail", [-3] = "timeout",
   [-2] = "create", [-1] = "skip", [0] = "fail", [1] = "pass",
}

-- GKYL_OUT_PREFIX is required by some Lua modules at load time.
GKYL_OUT_PREFIX = lfs.currentdir() .. "/" .. "queryrdb"

local log = Logger { logToFile = true }

-- Path of the configuration file written by 'runregression configure'.
local function computeConfFile()
   local gkeyllDir = GKYL_EXEC_PATH and GKYL_EXEC_PATH:match("^(.+)/bin$")
   if gkeyllDir then
      local f = io.open(gkeyllDir .. "/share/config.mak", "r")
      if f then
         for line in f:lines() do
            local prefix = line:match("^PREFIX%s*=%s*(.+)%s*$")
            if prefix then f:close(); return prefix .. "/gkeyll-results/runregression.config.lua" end
         end
         f:close()
      end
   end
   return os.getenv("HOME") .. "/runregression.config.lua"
end
local confFile = computeConfFile()

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
   local columns, nrow = sqlConn:exec("pragma table_info(RegressionMeta)")
   local hasRunMode = false
   for i = 1, nrow do
      if columns.name[i] == "run_mode" then hasRunMode = true; break end
   end
   if not hasRunMode then
      print("ERROR: legacy regression database schema; rerun runregression configure --drop-tables.")
      os.exit(1)
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
         run_mode  = t['run_mode'][i],
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
      dbData[i] = {
         name         = t['name'][i],
         test_type    = t['test_type'] and t['test_type'][i] or "?",
         status       = statusToString[tonumber(t['status'][i])],
         runtime      = t['runtime'][i],
         runlog       = t['runlog'][i],
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

   local fmt = "%-4s: %-20s %-30s %-13s %-5s %-5s %-5s"
   print(string.format(fmt, "ID", "Time-Stamp", "Changeset", "Mode", "Total", "Pass", "Fail"))
   for i, d in ipairs(dbMeta) do
      print(string.format(fmt, nrow - i + 1, d.tstamp, d.changeset,
         d.run_mode, tonumber(d.ntotal), tonumber(d.npass), tonumber(d.nfail)))
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

   -- Filtering predicate: every non-pass execution failure is a failure.
   local function shouldShow(d)
      if args.fail_only then
         return d.status == "fail" or d.status == "timeout"
            or d.status == "compile_fail" or d.status == "no_output" or d.status == "crash"
      end
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
      -- Tabular output with test type column.
      local nm  = maxNm + 2
      local fmt  = "%-4s: %-5s %-" .. nm .. "s %-13s %-9s"
      local fmt1 = "%-4s: %-5s %-" .. nm .. "s %-13s %.4g"
      print(string.format(fmt, "ID", "Type", "Name", "Status", "Run-Time"))
      for i, d in pairs(dbData) do
         if shouldShow(d) then
            print(string.format(fmt1, i, d.test_type, d.name, d.status, d.runtime))
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
      local fmt  = "%-20s %-30s %-13s %-5s %-13s %-9s"
      local fmt1 = "%-20s %-30s %-13s %-5s %-13s %.4g"
      print(string.format(fmt, "Time-Stamp", "Changeset", "Mode", "Type", "Status", "Run-Time"))
      for i = 1, nrow do
         local guid = dat['guid'][i]
         local tstamp, changeset, mode = sqlConn:rowexec(string.format(
            "select tstamp, GKYL_GIT_CHANGESET, run_mode from RegressionMeta where guid='%s'", guid))
         local stat    = statusToString[tonumber(dat['status'][i])]
         local ttype   = dat['test_type'] and dat['test_type'][i] or "?"
         print(string.format(fmt1, tstamp, changeset, mode, ttype, stat, dat['runtime'][i]))
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
