-- Smooth traveling wave problem for the 5-moment (Euler) equations.
-- Input parameters match the initial conditions found in entry JE22 of Ammar's Simulation Journal (https://ammar-hakim.org/sj/je/je22/je22-euler-2d.html#smooth-periodic-problem),
-- adapted from Section 4.1, from the article:
-- R. Liska and B. Wendroff (2003), "Comparison of Several Difference Schemes on 1D and 2D Test Problems for the Euler Equations",
-- SIAM Journal on Scientific Computing, Volume 25 (3): 995-1017.
-- https://epubs.siam.org/doi/10.1137/S1064827502402120

local Moments = G0.Moments
local Euler = G0.Moments.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
gas_gamma = 1.4 -- Adiabatic index.

u = 1.0 -- Fluid x-velocity.
v = -0.5 -- Fluid y-velocity
p = 1.0 -- Fluid pressure.

-- Simulation parameters.
Nx = 50 -- Cell count (x-direction).
Ny = 50 -- Cell count (y-direction).
Lx = 2.0 -- Domain size (x-direction).
Ly = 2.0 -- Domain size (y-direction).
cfl_frac = 0.9 -- CFL coefficient.

t_end = 4.0 -- Final simulation time.
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
  lower = { 0.0, 0.0 },
  upper = { Lx, Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction (x- and y-directions).
    
  -- Boundary conditions for configuration space.
  periodicDirs = { 1, 2 }, -- Periodic directions (x- and y-directions only).
  
  -- Fluid.
  fluid = Moments.Species.new {
    equation = Euler.new { gasGamma = gas_gamma },
    
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local rho = 1.0 + 0.2 * math.sin(pi * (x + y)) -- Fluid mass density.
      
      local mom_x = rho * u -- Fluid momentum density (x-direction).
      local mom_y = rho * v -- Fluid momentum density (y-direction).
      local mom_z = 0.0 -- Fluid momentum density (z-direction).
      local Etot = (p / (gas_gamma - 1.0)) + (0.5 * rho * (u * u + v * v)) -- Fluid total energy density.
      
      return rho, mom_x, mom_y, mom_z, Etot
    end,
  
    evolve = true -- Evolve species?
  }
}

-- Run application.
momentApp:run()
