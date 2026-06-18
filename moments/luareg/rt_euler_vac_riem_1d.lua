-- Sod-type shock tube test, with right-hand-side replaced with vacuum, for the 5-moment (Euler) equations.
-- Input parameters match the initial conditions in Section 2.6.2, with the contact discontinuity placed at x = 0.75 and the right-hand-side fluid replaced with vacuum, from the thesis:
-- A. Hakim (2006), "High Resolution Wave Propagation Schemes for Two-Fluid Plasma Simulations",
-- PhD Thesis, University of Washington.
-- https://www.aa.washington.edu/sites/aa/files/research/cpdlab/docs/PhDthesis_hakim.pdf

local Moments = G0.Moments
local Euler = G0.Moments.Eq.Euler

-- Physical constants (using normalized code units).
gas_gamma = 1.4 -- Adiabatic index.

rhol = 3.0 -- Left fluid mass density.
ul = 0.0 -- Left fluid velocity.
pl = 3.0 -- Left fluid pressure.

rhor = 1.0e-6 -- Right fluid mass density.
ur = 0.0 -- Right fluid velocity.
pr = 1.0e-6 -- Right fluid pressure.

-- Simulation parameters.
Nx = 512 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
cfl_frac = 0.9 -- CFL coefficient.

t_end = 0.05 -- Final simulation time.
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
  lower = { 0.25 },
  upper = { 0.25 + Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).
    
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).
  
  -- Fluid.
  fluid = Moments.Species.new {
    equation = Euler.new { gasGamma = gas_gamma },
    
    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local rho = 0.0
      local u = 0.0
      local p = 0.0

      if x < 0.75 then
        rho = rhol -- Fluid mass density (left).
        u = ul -- Fluid velocity (left).
        p = pl -- Fluid pressure (left).
      else
        rho = rhor -- Fluid mass density (right).
        u = ur -- Fluid velocity (right).
        p = pr -- Fluid pressure (right).
      end
  
      local mom_x = rho * u -- Fluid momentum density (x-direction).
      local mom_y = 0.0 -- Fluid momentum density (y-direction).
      local mom_z = 0.0 -- Fluid momentum density (z-direction).
      local Etot = (p / (gas_gamma - 1.0)) + (0.5 * rho * u * u) -- Fluid total energy density.
      
      return rho, mom_x, mom_y, mom_z, Etot
    end,
  
    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (x-direction).
  }
}

-- Run application.
momentApp:run()