-- Advection in specified electromagnetic fields for the 5-moment equations.
-- Input parameters match the initial conditions found in entry JE32 of Ammar's Simulation Journal (https://ammar-hakim.org/sj/je/je32/je32-vlasov-test-ptcl.html)
-- but with a rotation so that the oscillating electric field is in the z_hat direction and the background magnetic field in the x_hat direction. 
-- Solution is given by the resonant case, omega = Omega_c where Omega_c = q B/m is the cyclotron frequency. 

local Moments = G0.Moments
local Euler = G0.Moments.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using non-normalized physical units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

-- Simulation parameters.
n0 = 1.0 -- Reference density. 
vt = 1.0 -- Reference thermal velocity. 

-- External EM field parameters.
omega = 1.0 -- Oscillating electric field frequency normalized to cyclotron frequency.
B0 = 1.0 -- Reference magnetic field strength.

-- Derived physical quantities (using non-normalized physical units).
light_speed = 1.0 / math.sqrt(mu0 * epsilon0) -- Speed of light.

-- Simulation parameters.
Nx = 2 -- Cell count (x-direction).
Lx = 4.0 * math.pi -- Domain size (x-direction).
cfl_frac = 0.001 -- CFL coefficient. Set to be small to compare with analytic result. 
t_end = 100.0 -- Final simulation time.
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
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Electrons.
  elc = Moments.Species.new {
    charge = charge_elc, mass = mass_elc,
    equation = Euler.new { gasGamma = gas_gamma },

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]
    
      local rhoe = mass_elc * n0 -- Electron mass density.
      local mome_x = 0.0 -- Electron momentum density (x-direction).
      local mome_y = 0.0 -- Electron momentum density (y-direction).
      local mome_z = 0.0 -- Electron momentum density (z-direction).
      local Ee_tot = n0 * vt * vt * mass_elc / (gas_gamma - 1.0) -- Electron total energy density.
 
      return rhoe, mome_x, mome_y, mome_z, Ee_tot
    end,
    evolve = true, -- Evolve species?
  },

  -- Field.
  field = Moments.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,
    
    externalFieldInit = function (t, xn)
      local Ex = 0.0 -- External electric field (x-direction).
      local Ey = 0.0 -- External electric field (y-direction).
      local Ez = math.cos(omega * t) -- External electric field (z-direction).

      local Bx = B0 -- External magnetic field (x-direction).
      local By = 0.0 -- External magnetic field (y-direction).
      local Bz = 0.0 -- External magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz
    end,

    evolve = false, -- Evolve field?
    evolveExternalField = true, 
  }
}

-- Run application.
momentApp:run()
