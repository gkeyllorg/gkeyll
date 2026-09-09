local Moments = G0.Moments
local MHD = G0.Moments.Eq.MHD

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.

rho = 25.0 / (36.0 * pi) -- Fluid mass density.
p = 5.0 / (12.0 * pi) -- Fluid pressure.
B0 = 1.0 / math.sqrt(4.0 * pi) -- Reference magnetic field strength.

-- Simulation parameters.
Nx = 128 -- Cell count (x-direction).
Ny = 128 -- Cell count (y-direction).
Lx = 1.0 -- Domain size (x-direction).
Ly = 1.0 -- Domain size (y-direction).
cfl_frac = 0.95 -- CFL coefficient.

t_end = 0.5 -- Final simulation time.
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
  periodicDirs = { 1, 2 }, -- Periodic directions (x- and y-directions only).
  
  -- Fluid.
  fluid = Moments.Species.new {
    equation = MHD.new {
      gasGamma = gas_gamma,
      divergenceConstraint = G0.DivB.EightWaves
    },
    
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local u = math.sin(2.0 * pi * y) -- Fluid velocity (x-direction).
      local v = -math.sin(2.0 * pi * x) -- Fluid velocity (y-direction).
      local w = 0.0 -- Fluid velocity (z-direction).
      
      local Bx = B0 * math.sin(2.0 * pi * y) -- Magnetic field (x-direction).
      local By = B0 * math.sin(4.0 * pi * x) -- Magnetic field (y-direction).
      local Bz = 0.0 -- Magnetic field (z-direction).
    
      local mom_x = rho * u -- Fluid momentum density (x-direction).
      local mom_y = rho * v -- Fluid momentum density (y-direction).
      local mom_z = rho * w -- Fluid momentum density (z-direction).
      local Etot = (p / (gas_gamma - 1.0)) + (0.5 * rho * ((u * u) + (v * v) + (w * w))) + (0.5 * ((Bx * Bx) + (By * By) + (Bz * Bz))) -- Fluid total energy density.
      
      return rho, mom_x, mom_y, mom_z, Etot, Bx, By, Bz, 0.0
    end,
  
    evolve = true -- Evolve species?
  }
}

-- Run application.
momentApp:run()