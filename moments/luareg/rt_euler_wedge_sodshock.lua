-- 2D Sod-type shock tube test, with a wedge boundary condition, for the 5-moment (Euler) equations.
-- Input parameters are adapted from the initial conditions in Section 2.6.2, with the contract discontinuity placed at x = 0.75 rather than x = 0.5, from the thesis:
-- A. Hakim (2006), "High Resolution Wave Propagation Schemes for Two-Fluid Plasma Simulations",
-- PhD Thesis, University of Washington.
-- https://www.aa.washington.edu/sites/aa/files/research/cpdlab/docs/PhDthesis_hakim.pdf

local Moments = G0.Moments
local Euler = G0.Moments.Eq.Euler

-- Physical constants (using normalized code units).
gas_gamma = 1.4 -- Adiabatic index.

rhol = 3.0 -- Left/inner fluid mass density.
ul = 0.0 -- Left/inner fluid velocity.
pl = 3.0 -- Left/inner fluid pressure.

rhor = 1.0 -- Right/outer fluid mass density.
ur = 0.0 -- Right/outer fluid velocity.
pr = 1.0 -- Right/outer fluid pressure.

theta = 0.01 -- Wedge angle.

-- Simulation parameters.
Nr = 64 -- Cell count (radial direction).
Ntheta = 2 -- Cell count (angular direction).
Lr = 1.0 -- Domain size (radial direction).
Ltheta = theta -- Domain size (angular direction).
cfl_frac = 0.9 -- CFL coefficient.

t_end = 0.1 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

rloc = 0.5 * (0.25 + 1.25) -- Fluid boundary (radial coordinate).

momentApp = Moments.App.new {
  
  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.25, -0.5 * Ltheta },
  upper = { 0.25 + Lr, 0.5 * Ltheta },
  cells = { Nr, Ntheta },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction (radial and angular directions).
    
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Computational coordinates (r, theta) to physical coordinates (x, y).
  mapc2p = function (t, zc)
    local r, theta = zc[1], zc[2]

    local xp = { }
    xp[1] = r * math.cos(theta)
    xp[2] = r * math.sin(theta)

    return xp[1], xp[2]
  end,
  
  -- Fluid.
  fluid = Moments.Species.new {
    equation = Euler.new { gasGamma = gas_gamma },
    
    -- Initial conditions function.
    init = function (t, zc)
      local r = zc[1]

      local rho = 0.0
      local u = 0.0
      local p = 0.0
    
      if r < rloc then
        rho = rhol -- Fluid mass density (left/inner).
        u = ul -- Fluid velocity (left/inner).
        p = pl -- Fluid pressure (left/inner).
      else
        rho = rhor -- Fluid mass density (right/outer).
        u = ur -- Fluid velocity (right/outer).
        p = pr -- Fluid pressure (right/outer).
      end
  
      local mom_x = rho * u -- Fluid momentum density (x-direction).
      local mom_y = 0.0 -- Fluid momentum density (y-direction).
      local mom_z = 0.0 -- Fluid momentum density (z-direction).
      local Etot = (p / (gas_gamma - 1.0)) + (0.5 * rho * u * u) -- Fluid total energy density.
      
      return rho, mom_x, mom_y, mom_z, Etot
    end,
  
    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }, -- Copy boundary conditions (radial-direction).
    bcy = { G0.SpeciesBc.bcWedge, G0.SpeciesBc.bcWedge } -- Wedge boundary conditions (radial-direction).
  }
}

-- Run application.
momentApp:run()
