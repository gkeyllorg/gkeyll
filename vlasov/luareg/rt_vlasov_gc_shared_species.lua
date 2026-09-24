-- Test that a species definition can be shared by several apps without the
-- apps affecting each other: the callback contexts (including the coordinate
-- dimension a callback is called with) are copied into each app, so building a
-- 2D app from the same fluid species must not change how the live 1D app
-- calls its initial-condition function. Keep this test: it pins the per-app
-- ownership of callback contexts.

local Vlasov = G0.Vlasov

t_end = 0.2 -- Final simulation time (a few steps).
num_frames = 1 -- Number of output frames.

local observed_ndim = 0

local fluid = Vlasov.FluidSpecies.new {
  equation = Vlasov.Eq.LinearAdvection.new { },
  init = function (t, xn)
    observed_ndim = #xn
    return 1.0
  end,
  appAdvect = function (t, xn)
    return 1.0, 0.0, 0.0
  end,
}

vlasovApp = Vlasov.App.new {
  tEnd = t_end,
  nFrame = num_frames,
  lower = { 0.0 },
  upper = { 1.0 },
  cells = { 8 },
  polyOrder = 1,
  basis = "serendipity",
  periodicDirs = { 1 },
  skipField = true,
  fluid = fluid,
}
vlasovApp:apply_ic()
assert(observed_ndim == 1, "1D app called its initial condition with " .. observed_ndim .. " coordinates")

local app2d = Vlasov.App.new {
  tEnd = t_end,
  nFrame = num_frames,
  lower = { 0.0, 0.0 },
  upper = { 1.0, 1.0 },
  cells = { 4, 4 },
  polyOrder = 1,
  basis = "serendipity",
  periodicDirs = { 1, 2 },
  skipField = true,
  fluid = fluid,
}
app2d:apply_ic()
assert(observed_ndim == 2, "2D app called its initial condition with " .. observed_ndim .. " coordinates")

-- The first app must be unaffected by the second.
vlasovApp:apply_ic()
assert(observed_ndim == 1, "constructing a 2D app changed the live 1D app's callback to " .. observed_ndim .. "D")

vlasovApp:run()
