local Moments = G0.Moments
local Euler = G0.Moments.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
gas_gamma = 1.4 -- Adiabatic index.
grav = 0.1 -- Gravitational acceleration.

rho_top = 2.0 -- Top fluid mass density.
rho_bot = 1.0 -- Bottom fluid mass density.

p_ref = 0.01 -- Reference fluid pressure.
pert_max = 0.01 -- Maximum amplitude of initial perturbation.

-- Simulation parameters.
Nx = 50 -- Cell count (x-direction).
Ny = 200 -- Cell count (y-direction).
Lx = 1.0 / 6.0 -- Domain size (x-direction).
Ly = 1.0 -- Domain size (y-direction).
cfl_frac = 0.95 -- CFL coefficient.

t_end = 8.5 -- Final simulation time.
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
  
  -- Fluid.
  fluid = Moments.Species.new {
    equation = Euler.new { gasGamma = gas_gamma },
    
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local yloc = 0.5 + (pert_max * math.cos(2.0 * pi * (x / Lx)))

      local rho = 0.0
      local p = 0.0
    
      if y > yloc then
        rho = rho_top -- Fluid mass density (top).
        p = rho_top * grav * (1.0 - y) -- Fluid pressure (top).
      else
        rho = rho_bot -- Fluid mass density (bottom).
        p = (rho_top * grav * (1.0 - yloc)) + (rho_bot * grav * (yloc - y)) -- Fluid pressure (bottom).
      end
      
      local mom_x = 0.0 -- Fluid momentum density (x-direction).
      local mom_y = 0.0 -- Fluid momentum density (y-direction).
      local mom_z = 0.0 -- Fluid momentum density (z-direction).
      local Etot = (p + p_ref) / (gas_gamma - 1.0) -- Fluid total energy density.
      
      return rho, mom_x, mom_y, mom_z, Etot
    end,

    -- Applied acceleration function.
    appliedAcceleration = function (t, xn)
      local accel_x = 0.0 -- Applied acceleration (x-direction).
      local accel_y = -grav -- Applied acceleration (y-direction).
      local accel_z = 0.0 -- Applied acceleration (z-direction).

      return accel_x, accel_y, accel_z
    end,
  
    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcWall, G0.SpeciesBc.bcWall }, -- Wall boundary conditions (x-direction).
    bcy = { G0.SpeciesBc.bcWall, G0.SpeciesBc.bcWall } -- Wall boundary conditions (y-direction).
  }
}

-- Run application.
momentApp:run()