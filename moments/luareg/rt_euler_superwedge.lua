local Moments = G0.Moments
local Euler = G0.Moments.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
gas_gamma = 1.4 -- Adiabatic index.

rhol = 1.0 -- Left fluid mass density.
ul = 8.0 * math.sqrt(1.4) -- Left fluid velocity.
pl = 1.0 -- Left fluid pressure.

rhor = 1.0e-5 -- Right fluid mass density.
ur = 0.0 -- Right fluid velocity.
pr = 1.0e-5 -- Right fluid pressure.

theta = 15.0 * pi / 180.0 -- Wedge angle.

-- Simulation parameters.
Nx = 200 -- Cell count (x-direction).
Ny = 100 -- Cell count (y-direction).
Lx = 1.1 -- Domain size (x-direction).
Ly = 0.55 -- Domain size (y-direction).
cfl_frac = 0.9 -- CFL coefficient.

t_end = 1.0 -- Final simulation time.
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
  lower = { -0.1, 0.0 },
  upper = { -0.1 + Lx, Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction (x- and y-directions).
    
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Computational coordinates (x, y) to physical coordinates (x, y).
  mapc2p = function (t, zc)
    local x, y = zc[1], zc[2]

    local xp = { }

    xp[1] = x
    xp[2] = y
  
    if x > 0.0 then
      local yb = math.tan(theta) * x
      local a = (Ly - yb) / Ly
  
      xp[2] = a * y + yb
    end

    return xp[1], xp[2]
  end,
  
  -- Fluid.
  fluid = Moments.Species.new {
    equation = Euler.new { gasGamma = gas_gamma },
    
    -- Initial conditions function.
    init = function (t, zc)
      local x = zc[1]

      local rho = 0.0
      local u = 0.0
      local p = 0.0
      
      if x < -0.05 then
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
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }, -- Copy boundary conditions (x-direction).
    bcy = { G0.SpeciesBc.bcWall, G0.SpeciesBc.bcWall } -- Wall boundary conditions (y-direction).
  }
}

-- Run application.
momentApp:run()
