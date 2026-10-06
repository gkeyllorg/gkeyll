-- Forward Jenkins CI commands to the shell script.

local function shell_quote(value)
   return "'" .. value:gsub("'", "'\\''") .. "'"
end

local commands = {
   scan = true, start = true, run = true, active = true,
   recent = true, info = true, artifact = true, status = true, follow = true, abort = true,
}

local is_help =
   (#GKYL_COMMANDS_L == 1
      and (GKYL_COMMANDS_L[1] == "-h" or GKYL_COMMANDS_L[1] == "--help"))
   or (#GKYL_COMMANDS_L == 2 and commands[GKYL_COMMANDS_L[1]]
      and (GKYL_COMMANDS_L[2] == "-h" or GKYL_COMMANDS_L[2] == "--help"))

local command = "./ci/jenkins/gkeyll-ci.sh"
for _, value in ipairs(GKYL_COMMANDS_L) do
   command = command .. " " .. shell_quote(value)
end

if is_help then
   command = command .. " | sed 's/gkeyll-ci\\.sh/gkeyll jeci/g'"
end

local result, _, status = os.execute(command)
if type(result) == "number" then
   -- LuaJIT follows Lua 5.1 and returns the raw POSIX wait status.
   if result ~= 0 then os.exit(math.floor(result / 256)) end
elseif result ~= true then
   os.exit(status or 1)
end
