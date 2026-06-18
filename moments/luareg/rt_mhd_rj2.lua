local Moments = G0.Moments
local MHD = G0.Moments.Eq.MHD

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.

rhol = 1.08 -- Left fluid mass density.
ul = 1.2 -- Left fluid velocity (x-direction).
vl = 0.01 -- Left fluid velocity (y-direction).
wl = 0.5 -- Left fluid velocity (z-direction).
pl = 0.95 -- Left fluid pressure.

Bx_l = 2.0 / math.sqrt(4.0 * pi) -- Left magnetic field (x-direction).
By_l = 3.6 / math.sqrt(4.0 * pi) -- Left magnetic field (y-direction).
Bz_l = 2.0 / math.sqrt(4.0 * pi) -- Left magnetic field (z-direction).

rhor = 1.0 -- Right fluid mass density.
ur = 0.0 -- Right fluid velocity (x-direction).
vr = 0.0 -- Right fluid velocity (y-direction).
wr = 0.0 -- Right fluid velocity (z-direction).
pr = 1.0 -- Right fluid pressure.

Bx_r = 2.0 / math.sqrt(4.0 * pi) -- Right magnetic field (x-direction).
By_r = 4.0 / math.sqrt(4.0 * pi) -- Right magnetic field (y-direction).
Bz_r = 2.0 / math.sqrt(4.0 * pi) -- Right magnetic field (z-direction).

-- Simulation parameters.
 Nx = 512 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
cfl_frac = 0.8 -- CFL coefficient.

t_end = 0.2 -- Final simulation time.
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
      local v = 0.0
      local w = 0.0
      local p = 0.0
    
      local Bx = 0.0
      local By = 0.0
      local Bz = 0.0 -- Magnetic field (z-direction).

      if x < 0.5 then
        rho = rhol -- Fluid mass density (left).
        u = ul -- Fluid velocity (x-direction, left).
        v = vl -- Fluid velocity (y-direction, left).
        w = wl -- Fluid velocity (z-direction, left).
        p = pl -- Fluid pressure (left).
    
        Bx = Bx_l -- Magnetic field (x-direction, left).
        By = By_l -- Magnetic field (y-direction, left).
        Bz = Bz_l -- Magnetic field (z-direction, left).
      else
        rho = rhor -- Fluid mass density (right).
        u = ur -- Fluid velocity (x-direction, right).
        v = vr -- Fluid velocity (y-direction, right).
        w = wr -- Fluid velocity (z-direction, right).
        p = pr -- Fluid pressure (right).
    
        Bx = Bx_r -- Magnetic field (x-direction, right).
        By = By_r -- Magnetic field (y-direction, right).
        Bz = Bz_r -- Magnetic field (z-direction, right).
      end
    
      local mom_x = rho * u -- Fluid momentum density (x-direction).
      local mom_y = rho * v -- Fluid momentum density (y-direction).
      local mom_z = rho * w -- Fluid momentum density (z-direction).
      local Etot = (p / (gas_gamma - 1.0)) + (0.5 * rho * ((u * u) + (v * v) + (w * w))) + (0.5 * ((Bx * Bx) + (By * By) + (Bz * Bz))) -- Fluid total energy density.
      
      return rho, mom_x, mom_y, mom_z, Etot, Bx, By, Bz, 0.0
    end,
  
    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (x-direction).
  }
}

-- Run application.
momentApp:run()