-- Discontinuous initial data test, with monotonicity-preserving reconstruction, for the linear advection equation.
-- Input parameters match the initial conditions in Section 4.1, Example 2, from the article:
-- A. Suresh and H. T. Huynh (1997), "Accurate Monotonicity-Preserving Schemes with Runge-Kutta Time Stepping",
-- Journal of Computational Physics, Volume 136 (1): 83-99.
-- https://www.sciencedirect.com/science/article/pii/S0021999197957454

local Moments = G0.Moments
local LinearAdvection = G0.Moments.Eq.LinearAdvection

-- Simulation parameters.
Nx = 200 -- Cell count (x-direction).
Lx = 2.0 -- Domain size (x-direction).
v_advect = 1.0 -- Advection velocity.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 20.0 -- Final simulation time.
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
  lower = { -0.5 * Lx },
  upper = { 0.5 * Lx },
  cells = { Nx },
  cflFrac = cfl_frac,
  schemeType = G0.SchemeType.MP,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).
    
  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).
  
  -- Advected quantity.
  fluid = Moments.Species.new {
    equation = LinearAdvection.new { advectionSpeed = v_advect },
    
    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local f = 0.0

      if -0.8 <= x and x <= -0.6 then
        f = math.exp(-math.log(2.0) * (x + 0.7) * (x + 0.7) / 0.0009) -- Advected quantity (between -0.8 and -0.6).
      elseif -0.4 <= x and x <= -0.2 then
        f = 1.0 -- Advected quantity (between -0.4 and -0.2).
      elseif 0 <= x and x <= 0.2 then
        f = 1.0 - math.abs(10.0 * (x - 0.1)) -- Advected quantity (between 0 and 0.2).
      elseif 0.4 <= x and x <= 0.6 then
        f = math.sqrt(1.0 - 100.0 * (x - 0.5) * (x - 0.5)) -- Advected quantity (between 0.4 and 0.6).
      end
      
      return f
    end,
  
    evolve = true -- Evolve species?
  }
}

-- Run application.
momentApp:run()