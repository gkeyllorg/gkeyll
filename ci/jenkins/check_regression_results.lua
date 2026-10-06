-- Gkyl ------------------------------------------------------------------------
--
-- Jenkins helper: evaluate regression results written by
-- 'gkeyll runregression run ... check' and fail (os.exit(1)) if any test
-- did not pass, unless it's listed in the candidate acknowledgment file.
--
-- Usage: gkeyll ci/jenkins/check_regression_results.lua <resultsDir> [ackFile] [summaryFile]
--   <resultsDir>  the gkeyll-results/ directory written by 'runregression
--                 configure' (i.e. <prefix>/gkeyll-results).
--   [ackFile]     optional path to a text file listing tests (one per line,
--                 '#' comments allowed) whose diff is expected/acknowledged for
--                 this candidate and should not fail the build. A test may be
--                 named as stored in the database ("moments/creg/rt_euler_sodshock"),
--                 as "<layer>/<basename>" ("moments/rt_euler_sodshock"), or by
--                 basename alone. See ci/jenkins/expected_regression_diffs.txt.
--   [summaryFile] optional machine-readable pass/acknowledged/failure counts.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local sql = require "sqlite3"

local LAYERS = { "moments", "vlasov", "gyrokinetic", "pkpm" }

local statusToString = {
   [-6] = "crash", [-5] = "no_output", [-4] = "compile_fail", [-3] = "timeout", [-2] = "create", [-1] = "skip",
   [0] = "fail", [1] = "pass",
}
-- Statuses that block a candidate build unless explicitly acknowledged.
local BAD_STATUSES = { [0] = true, [-3] = true, [-4] = true, [-5] = true, [-6] = true }

local resultsDir = GKYL_COMMANDS_L[1]
local ackFile = GKYL_COMMANDS_L[2]
local summaryFile = GKYL_COMMANDS_L[3]

if not resultsDir then
   print("Usage: gkeyll check_regression_results.lua <resultsDir> [ackFile] [summaryFile]")
   os.exit(1)
end

-- Parse the acknowledgment file into a set of "<layer>/<name>" entries.
-- Blank lines and '#' comments (including trailing '# reason' text) are
-- ignored.
local acked = {}
if ackFile then
   local f = io.open(ackFile, "r")
   if f then
      for line in f:lines() do
         local entry = line:gsub("#.*$", ""):match("^%s*(.-)%s*$")
         if entry and entry ~= "" then acked[entry] = true end
      end
      f:close()
   end
end

local npass, ackedHits, unacked = 0, {}, {}
local layerCounts = {}  -- layerCounts[layer] = { passed, acked, failed }

-- Names are stored layer-qualified ("moments/creg/rt_x"); accept the shorter
-- spellings in the acknowledgment file too.
local function isAcked(layer, name)
   local base = name:match("([^/]+)$") or name
   return acked[name] or acked[layer .. "/" .. base] or acked[base]
end

-- Per-file comparison lines from the stored run log, for the CI report.
local function comparisonLines(runlog, maxLines)
   local lines = {}
   local body = runlog and runlog:match("%-%-%- Comparison failures %-%-%-\n(.*)$")
   if body then
      for line in body:gmatch("[^\n]+") do
         if #lines >= maxLines then table.insert(lines, "..."); break end
         table.insert(lines, line)
      end
   end
   return lines
end

for _, layer in ipairs(LAYERS) do
   local dbPath = string.format("%s/%s/regressiondb", resultsDir, layer)
   local f = io.open(dbPath, "r")
   if not f then
      print(string.format("WARNING: no regressiondb for layer '%s' at %s (skipped)", layer, dbPath))
   else
      f:close()
      local conn = sql.open(dbPath)
      local cols, ncols = conn:exec("pragma table_info(RegressionMeta)")
      local hasRunMode = false
      for i = 1, ncols do
         if cols.name[i] == "run_mode" then hasRunMode = true; break end
      end
      if not hasRunMode then
         error("Legacy regression database schema at " .. dbPath
            .. "; rerun runregression configure --drop-tables.")
      end
      local guid = conn:rowexec("select guid from RegressionMeta order by rowid desc limit 1")
      layerCounts[layer] = { passed = 0, acked = 0, failed = 0 }
      if guid then
         local t, nrow = conn:exec(string.format(
            "select name, status, runlog from RegressionData where guid=='%s'", guid))
         for i = 1, nrow do
            local status = tonumber(t.status[i])
            local key = t.name[i]
            if status == 1 then
               npass = npass + 1
               layerCounts[layer].passed = layerCounts[layer].passed + 1
            elseif BAD_STATUSES[status] then
               if isAcked(layer, key) then
                  table.insert(ackedHits, key)
                  layerCounts[layer].acked = layerCounts[layer].acked + 1
               else
                  layerCounts[layer].failed = layerCounts[layer].failed + 1
                  table.insert(unacked, {
                     key = key, status = statusToString[status] or tostring(status),
                     detail = comparisonLines(t.runlog[i], 10),
                  })
               end
            end
         end
      end
      conn:close()
   end
end

print(string.format(
   "Regression results: %d passed, %d acknowledged diff(s), %d unacknowledged failure(s)",
   npass, #ackedHits, #unacked))

if summaryFile then
   local summary, message = io.open(summaryFile, "w")
   if not summary then
      print(string.format("ERROR: cannot write regression summary %s: %s", summaryFile, message))
      os.exit(1)
   end
   summary:write(string.format(
      "c_regression_passed=%d\nc_regression_acknowledged=%d\nc_regression_unacknowledged=%d\n",
      npass, #ackedHits, #unacked))
   -- One line per test so the GitHub report (ci/jenkins/github_report.py)
   -- can list failures without opening the SQLite database.
   for _, layer in ipairs(LAYERS) do
      local c = layerCounts[layer]
      if c then
         summary:write(string.format("c_regression_layer=%s:%d:%d:%d\n", layer, c.passed, c.acked, c.failed))
      end
   end
   for _, u in ipairs(unacked) do
      summary:write(string.format("c_regression_failure=%s:%s\n", u.key, u.status))
      for _, line in ipairs(u.detail) do
         summary:write(string.format("c_regression_failure_detail=%s|%s\n", u.key, line))
      end
   end
   for _, key in ipairs(ackedHits) do
      summary:write(string.format("c_regression_acknowledged_test=%s\n", key))
   end
   summary:close()
end

if #ackedHits > 0 then
   print(string.format("Acknowledged (per %s):", tostring(ackFile)))
   for _, key in ipairs(ackedHits) do print("  " .. key) end
end

if #unacked > 0 then
   print("")
   print("UNACKNOWLEDGED FAILURES:")
   for _, u in ipairs(unacked) do
      print(string.format("  %s [%s]", u.key, u.status))
      for _, line in ipairs(u.detail) do print("      " .. line) end
   end
   print("")
   print("If any of these are an expected/intentional change (e.g. a physics")
   print("or algorithm change), add a line for each to")
   print("ci/jenkins/expected_regression_diffs.txt with a short reason, e.g.:")
   for _, u in ipairs(unacked) do
      print(string.format("  %s   # <why this changed>", u.key))
   end
   os.exit(1)
end

os.exit(0)
