-- Garbage-collector test for the App's ownership of its inputs. The App keeps
-- the species/field userdata alive through a private anchor table built from
-- the input table at construction, so editing or dropping entries of the
-- caller's table afterwards must not free the userdata a live app was built
-- from. The weak reference below detects collection directly. Keep this test:
-- the app's callback contexts are its own copies, and the Lua functions are
-- rooted by registry references, so the userdata anchor is what this pins.

local Vlasov = G0.Vlasov

t_end = 0.2 -- Final simulation time (a few steps).
num_frames = 1 -- Number of output frames.

local inp = {
  tEnd = t_end,
  nFrame = num_frames,
  lower = { 0.0 },
  upper = { 1.0 },
  cells = { 8 },
  polyOrder = 1,
  basis = "serendipity",
  periodicDirs = { 1 },
  skipField = true,

  fluid = Vlasov.FluidSpecies.new {
    equation = Vlasov.Eq.LinearAdvection.new { },
    init = function (t, xn)
      return 1.0 + 0.1 * math.sin(2.0 * math.pi * xn[1])
    end,
    appAdvect = function (t, xn)
      return 1.0, 0.0, 0.0
    end,
  },
}

local weak = setmetatable({ inp.fluid }, { __mode = "v" })

vlasovApp = Vlasov.App.new(inp)

-- Edit the caller's table after construction: the app must not depend on it.
inp.fluid = nil
inp = nil
collectgarbage("collect")
collectgarbage("collect")
assert(weak[1] ~= nil, "live app lost the userdata owning its callbacks")

vlasovApp:run()
