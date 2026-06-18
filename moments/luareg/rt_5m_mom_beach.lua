-- Propagation into a plasma wave beach test for the 5-moment equations.
-- Input parameters match the initial conditions found in entry JE8 of Ammar's Simulation Journal (https://ammar-hakim.org/sj/je/je8/je8-plasmabeach.html), adapted from Section III. A. of the article:
-- D. N. Smithe (2007), "Finite-difference time-domain simulation of fusion plasmas at radiofrequency time scales",
-- Physics of Plasmas, Volume 14 (5): 056104.
-- https://pubs.aip.org/aip/pop/article/14/5/056104/929539/Finite-difference-time-domain-simulation-of-fusion

local Moments = G0.Moments
local Euler = G0.Moments.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using non-normalized physical units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.
epsilon0 = 8.854187817620389850536563031710750260608e-12 -- Permittivity of free space.
mu0 = 12.56637061435917295385057353311801153679e-7 -- Permeability of free space.
mass_elc = 9.10938215e-31 -- Electron mass.
charge_elc = -1.602176487e-19 -- Electron charge.

J0 = 1.0e-12 -- Reference current density (Amps / m^3).

-- Derived physical quantities (using non-normalized physical units).
light_speed = 1.0 / math.sqrt(mu0 * epsilon0) -- Speed of light.

-- Simulation parameters.
Nx = 400 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
Lx100 = Lx / 100.0 -- Domain size over 100 (x-direction).
x_last_edge = Lx - Lx / Nx -- Location of center of last upper cell (low density side).
cfl_frac = 0.95 -- CFL coefficient.
t_end = 5.0e-9 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

deltaT = Lx100 / light_speed -- Arbitrary constant, with units of time.
factor = deltaT * deltaT * charge_elc * charge_elc / (mass_elc * epsilon0) -- Numerical factor for calculation of electron number density.
omega_drive = pi / 10.0 / deltaT -- Drive current angular frequency.

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

  -- Electrons.
  elc = Moments.Species.new {
    charge = charge_elc, mass = mass_elc,
    equation = Euler.new { gasGamma = gas_gamma },

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local omegaPdt = 25.0 * (1.0 - x) * (1.0 - x) * (1.0 - x) * (1.0 - x) * (1.0 - x) -- Plasma frequency profile.
      local ne = omegaPdt * omegaPdt / factor -- Electron number density.
    
      local rhoe = mass_elc * ne -- Electron mass density.
      local mome_x = 0.0 -- Electron momentum density (x-direction).
      local mome_y = 0.0 -- Electron momentum density (y-direction).
      local mome_z = 0.0 -- Electron momentum density (z-direction).
      local Ee_tot = ne * (-charge_elc) / (gas_gamma - 1.0) -- Electron total energy density.
 
      return rhoe, mome_x, mome_y, mome_z, Ee_tot
    end,

    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (x-direction).
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

    -- Applied current function.
    appliedCurrent = function (t, xn)
      local x = xn[1]

      local app_x = 0.0 -- Applied current (x-direction).
      local app_y = 0.0
      local app_z = 0.0 -- Applied current (z-direction).

      if x > x_last_edge then
        app_y = -J0 * math.sin(omega_drive * t) -- Applied current (y-direction, right).
      else
        app_y = 0.0 -- Applied current (y-direction, left).
      end

      return app_x, app_y, app_z
    end,
    evolveAppliedCurrent = true, -- Evolve applied current.

    evolve = true, -- Evolve field?
    bcx = { G0.FieldBc.bcCopy, G0.FieldBc.bcCopy } -- Copy boundary conditions (x-direction).
  }
}

-- Run application.
momentApp:run()
