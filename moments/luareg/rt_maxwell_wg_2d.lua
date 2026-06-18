local Moments = G0.Moments

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.

-- Derived physical quantities (using normalized code units).
omega = math.sqrt(2.0) * pi -- Frequency of mode.
t_period = 2.0 * pi / omega -- Period of mode.

-- Simulation parameters.
Nx = 70 -- Cell count (x-direction).
Ny = 50 -- Cell count (y-direction).
Lx = 7.0 -- Domain size (x-direction).
Ly = 5.0 -- Domain size (y-direction).
cfl_frac = 1.0 -- CFL coefficient.

t_end = 10.0 * t_period -- Final simulation time.
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
  periodicDirs = { }, -- Periodic directions (none).
      
  -- Field.
  field = Moments.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = math.sin(pi * x) * math.sin(pi * y) -- Total electric field (z-direction).
    
      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = math.cos(pi * x) * math.cos(pi * y) -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
    bcx = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall }, -- PEC wall boundary conditions (x-direction).
    bcy = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall }, -- PEC wall boundary conditions (y-direction).
    limiter = G0.WaveLimiter.NoLimiter
  }
}
  
-- Run application.
momentApp:run()
  