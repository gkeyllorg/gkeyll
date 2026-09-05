-- Gkyl ------------------------------------------------------------------------
--
-- Jenkins helper: evaluate the MOAT regression results written by
-- 'gkeyll runregression run --moat check' and fail (os.exit(1)) if any test
-- did not pass, unless it's listed in the per-PR acknowledgment file.
--
-- Usage: gkeyll ci/jenkins/check_regression_results.lua <resultsDir> [ackFile]
--   <resultsDir>  the gkeyll-results/ directory written by 'runregression
--                 configure' (i.e. <prefix>/gkeyll-results).
--   [ackFile]     optional path to a text file listing "<layer>/<name>"
--                 entries (one per line, '#' comments allowed) whose diff is
--                 expected/acknowledged for this PR and should not fail the
--                 build. See ci/jenkins/expected_regression_diffs.txt.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local sql = require "sqlite3"

local LAYERS = { "moments", "vlasov", "gyrokinetic", "pkpm" }

local statusToString = {
   [-4] = "compile_fail", [-3] = "timeout", [-2] = "create", [-1] = "skip",
   [0] = "fail", [1] = "pass",
}
-- Statuses that block a PR build unless explicitly acknowledged.
local BAD_STATUSES = { [0] = true, [-3] = true, [-4] = true }

local resultsDir = GKYL_COMMANDS_L[1]
local ackFile = GKYL_COMMANDS_L[2]

if not resultsDir then
   print("Usage: gkeyll check_regression_results.lua <resultsDir> [ackFile]")
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

for _, layer in ipairs(LAYERS) do
   local dbPath = string.format("%s/%s/regressiondb", resultsDir, layer)
   local f = io.open(dbPath, "r")
   if not f then
      print(string.format("WARNING: no regressiondb for layer '%s' at %s (skipped)", layer, dbPath))
   else
      f:close()
      local conn = sql.open(dbPath)
      local guid = conn:rowexec("select guid from RegressionMeta order by rowid desc limit 1")
      if guid then
         local t, nrow = conn:exec(string.format(
            "select name, status from RegressionData where guid=='%s'", guid))
         for i = 1, nrow do
            local status = tonumber(t.status[i])
            local key = layer .. "/" .. t.name[i]
            if status == 1 then
               npass = npass + 1
            elseif BAD_STATUSES[status] then
               if acked[key] then
                  table.insert(ackedHits, key)
               else
                  table.insert(unacked, { key = key, status = statusToString[status] or tostring(status) })
               end
            end
         end
      end
      conn:close()
   end
end

print(string.format(
   "MOAT results: %d passed, %d acknowledged diff(s), %d unacknowledged failure(s)",
   npass, #ackedHits, #unacked))

if #ackedHits > 0 then
   print(string.format("Acknowledged (per %s):", tostring(ackFile)))
   for _, key in ipairs(ackedHits) do print("  " .. key) end
end

if #unacked > 0 then
   print("")
   print("UNACKNOWLEDGED FAILURES:")
   for _, u in ipairs(unacked) do
      print(string.format("  %s [%s]", u.key, u.status))
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
