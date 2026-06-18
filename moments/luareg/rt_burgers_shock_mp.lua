-- Sinusoidal initial data test, with monotonicity-preserving reconstruction, for the inviscid Burgers' equation.
-- Input parameters correspond to a sine wave, equivalent to the initial conditions in Section 4, Example 2, but with periodicity 1 rather than 2, from the article:
-- X-d. Liu and S. Osher (1996), "Nonoscillatory High Order Accurate Self-Similar Maximum Principle Satisfying Shock Capturing Schemes I",
-- SIAM Journal on Numerical Analysis, Volume 33 (2): 760-779.
-- https://epubs.siam.org/doi/10.1137/0733038

local Moments = G0.Moments
local Burgers = G0.Moments.Eq.Burgers

-- Mathematical constants (dimensionless).
pi = math.pi

-- Simulation parameters.
Nx = 128 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
cfl_frac = 0.9 -- CFL coefficient.

t_end = 0.4 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

momentApp = Moments.App.new {
  
  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.0 },
  upper = { Lx },
  cells = { Nx },
  cflFrac = cfl_frac,
  schemeType = G0.SchemeType.MP,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).
    
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).
  
  -- Advected quantity.
  fluid = Moments.Species.new {
    equation = Burgers.new { },
    
    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local f = math.sin(2.0 * pi * x) -- Advected quantity.
      
      return f
    end,
  
    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (x-direction).
  }
}

-- Run application.
momentApp:run()