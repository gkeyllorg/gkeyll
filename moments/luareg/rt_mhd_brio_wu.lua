local Moments = G0.Moments
local MHD = G0.Moments.Eq.MHD

-- Physical constants (using normalized code units).
gas_gamma = 2.0 -- Adiabatic index.

rhol = 1.0 -- Left fluid mass density.
ul = 0.0 -- Left fluid velocity.
pl = 1.0 -- Left fluid pressure.
Bx_l = 0.75 -- Left magnetic field (x-direction).
By_l = 1.0 -- Left magnetic field (y-direction).

rhor = 0.125 -- Right fluid mass density.
ur = 0.0 -- Right fluid velocity.
pr = 0.1 -- Right fluid pressure.
Bx_r = 0.75 -- Right magnetic field (x-direction).
By_r = -1.0 -- Right magnetic field (y-direction).

-- Simulation parameters.
Nx = 400 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
cfl_frac = 0.95 -- CFL coefficient.

t_end = 0.1 -- Final simulation time.
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

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).
    
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).
  
  -- Fluid.
  fluid = Moments.Species.new {
    equation = MHD.new {
      gasGamma = gas_gamma,
      divergenceConstraint = G0.DivB.None
    },
    
    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local rho = 0.0
      local u = 0.0
      local p = 0.0

      local Bx = 0.0
      local By = 0.0
      local Bz = 0.0 -- Magnetic field (x-direction).

      if x < 0.5 then
        rho = rhol -- Fluid mass density (left).
        u = ul -- Fluid velocity (left).
        p = pl -- Fluid pressure (left).

        Bx = Bx_l -- Magnetic field (x-direction, left).
        By = By_l -- Magnetic field (y-direction, left).
      else
        rho = rhor -- Fluid mass density (right).
        u = ur -- Fluid velocity (right).
        p = pr -- Fluid pressure (right).

        Bx = Bx_r -- Magnetic field (x-direction, right).
        By = By_r -- Magnetic field (y-direction, right).
      end
  
      local mom_x = rho * u -- Fluid momentum density (x-direction).
      local mom_y = 0.0 -- Fluid momentum density (y-direction).
      local mom_z = 0.0 -- Fluid momentum density (z-direction).
      local Etot = (p / (gas_gamma - 1.0)) + (0.5 * rho * u * u) + (0.5 * ((Bx * Bx) + (By * By) + (Bz * Bz))) -- Fluid total energy density.
      
      return rho, mom_x, mom_y, mom_z, Etot, Bx, By, Bz, 0.0
    end,
  
    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (x-direction).
  }
}

-- Run application.
momentApp:run()