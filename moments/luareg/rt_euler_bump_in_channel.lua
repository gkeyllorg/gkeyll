-- Supersonic flow over a blunt body test for the Euler equations.
-- Input parameters match the initial conditions found in Entry JE24 of Ammar's Simulation Journal (https://ammar-hakim.org/sj/je/je24/je24-euler-embedded-bc.html#supersonic-flow-over-blunt-body).
-- Problem setup is similar to Section V. C. of the article:
-- P. V. Tota and Z. J. Wang (2007), "Meshfree Euler Solver using local Radial Basis Functions for inviscid Compressible Flows",
-- 18th AIAA Computational Fluid Dynamics Conference (25th June 2007 - 28th June 2007, Miami, Florida).
-- https://arc.aiaa.org/doi/10.2514/6.2007-4581

local Moments = G0.Moments
local Euler = G0.Moments.Eq.Euler

-- Physical constants (using normalized code units).
gas_gamma = 1.4 -- Adiabatic index.

rhol = 1.0 -- Left fluid mass density.
pl = 1.0 -- Left fluid pressure.

rhor = 1e-5 -- Right fluid mass density.
pr = 1e-5 -- Right fluid pressure.

-- Derived physical quantities (using normalized code units).
cs = math.sqrt(gas_gamma) -- Fluid sound speed

ul = 2.0 * cs -- Left fluid velocity
ur = 0.0 -- Right fluid velocity

-- Simulation parameters.
Nx = 150 -- Cell count (x-direction).
Ny = 75 -- Cell count (y-direction).
cfl_frac = 0.9 -- CFL coefficient.

t_end = 10.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

inlet_xlen = 1.0 -- Length of inlet (x-direction).
bump_xlen = 2.0 -- Length of bump (x-direction).
height = 2.0 -- Height of channel.
R = 5.0 -- Radius of bump.

half_xlen = inlet_xlen + 0.5 * bump_xlen -- Mid-poof bump (x-direction).

Lx = 2.0 * half_xlen -- Domain size (x-direction).
Ly = height -- Domain size (y-direction).

momentApp = Moments.App.new {
  
  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { -0.5 * Lx, 0.0 },
  upper = { 0.5 * Lx, Ly },
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
    local zeta_min = math.sqrt((R * R) - ((0.5 * bump_xlen) * (0.5 * bump_xlen)))

    xp[1] = x
    xp[2] = y
  
    if math.abs(x) < 0.5 * bump_xlen then
      local eta = x
      local zeta = math.sqrt((R * R) - (eta * eta))
      local yb = zeta - zeta_min
  
      xp[2] = (height - yb) * y / height + yb
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
    
      if x < -half_xlen * (1.0 - 0.25) then
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
    bcy = { G0.SpeciesBc.bcWall, G0.SpeciesBc.bcCopy } -- Wall/copy boundary conditions (y-direction).
  }
}

-- Run application.
momentApp:run()
