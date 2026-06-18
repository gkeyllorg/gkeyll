local Moments = G0.Moments

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.

-- Simulation parameters.
Nx = 128 -- Cell count (x-direction).
Ny = 128 -- Cell count (y-direction).
Lx = 2.0 -- Domain size (x-direction).
Ly = 2.0 -- Domain size (y-direction).
cfl_frac = 1.0 -- CFL coefficient.

t_end = 3.0 -- Final simulation time.
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
  lower = { -0.5 * Lx, -0.5 * Ly },
  upper = { 0.5 * Lx, 0.5 * Ly },
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

      local r_sq = (x * x) + (y * y)

      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = math.exp(-25.0 * r_sq) -- Total electric field (z-direction).
      
      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
    limiter = G0.WaveLimiter.NoLimiter,
    bcx = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall }, -- PEC wall boundary conditions (x-direction).
    bcy = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall } -- PEC wall boundary conditions (y-direction).
  }
}
  
-- Run application.
momentApp:run()
  