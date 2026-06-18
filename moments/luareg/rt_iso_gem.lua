-- Geospace Environmental Modeling (GEM) magnetic reconnection test for the isothermal Euler equations.
-- Input parameters match the equilibrium and initial conditions in Section 2, from the article:
-- J. Birn et al. (2001), "Geospace Environmental Modeling (GEM) Magnetic Reconnection Challenge",
-- Journal of Geophysical Research: Space Physics, Volume 106 (A3): 3715-3719.
-- https://agupubs.onlinelibrary.wiley.com/doi/10.1029/1999JA900449

local Moments = G0.Moments
local IsoEuler = G0.Moments.Eq.IsoEuler

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_ion = 1.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.
mass_elc = 1.0 / 25.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
Ti_over_Te = 5.0 -- Ion temperature / electron temperature.
lambda = 0.5 -- Wavelength.
n0 = 1.0 -- Reference number density.
nb_over_n0 = 0.2 -- Background number density / reference number density.
B0 = 0.1 -- Reference magnetic field strength.
beta = 1.0 -- Plasma beta.

-- Derived physical quantities (using normalized code units).
psi0 = 0.1 * B0 -- Reference magnetic scalar potential.

Ti_frac = Ti_over_Te / (1.0 + Ti_over_Te) -- Fraction of total temperature from ions.
Te_frac = 1.0 / (1.0 + Ti_over_Te) -- Fraction of total temperature from electrons.
T_tot = beta * (B0 * B0) / 2.0 / n0 -- Total temperature.

vti = math.sqrt(T_tot * Ti_frac / mass_ion) -- Ion thermal velocity.
vte = math.sqrt(T_tot * Te_frac / mass_elc) -- Electron thermal velocity.

-- Simulation parameters.
Nx = 64 -- Cell count (x-direction).
Ny = 32 -- Cell count (y-direction).
Lx = 25.6 -- Domain size (x-direction).
Ly = 12.8 -- Domain size (y-direction).
cfl_frac = 1.0 -- CFL coefficient.

t_end = 250.0 -- Final simulation time.
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
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Electrons.
  elc = Moments.Species.new {
    charge = charge_elc, mass = mass_elc,
    equation = IsoEuler.new { vThermal = vte },

    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local sech_sq = (1.0 / math.cosh(y / lambda)) * (1.0 / math.cosh(y / lambda)) -- Hyperbolic secant squared.

      local n = n0 * (sech_sq + nb_over_n0) -- Total number density.
      local Jz = -(B0 / lambda) * sech_sq -- Total current density (z-direction).

      local rhoe = n * mass_elc -- Electron mass density.
      local mome_x = 0.0 -- Electron momentum density (x-direction).
      local mome_y = 0.0 -- Electron momentum density (y-direction).
      local mome_z = (mass_elc / charge_elc) * Jz * Te_frac -- Electron momentum density (z-direction).
 
      return rhoe, mome_x, mome_y, mome_z
    end,

    evolve = true, -- Evolve species?
    bcy = { G0.SpeciesBc.bcWall, G0.SpeciesBc.bcWall } -- Wall boundary conditions (y-direction).
  },

  -- Ions.
  ion = Moments.Species.new {
    charge = charge_ion, mass = mass_ion,
    equation = IsoEuler.new { vThermal = vti },
  
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local sech_sq = (1.0 / math.cosh(y / lambda)) * (1.0 / math.cosh(y / lambda)) -- Hyperbolic secant squared.

      local n = n0 * (sech_sq + nb_over_n0) -- Total number density.
      local Jz = -(B0 / lambda) * sech_sq -- Total current density (z-direction).

      local rhoi = n * mass_ion -- Ion mass density.
      local momi_x = 0.0 -- Ion momentum density (x-direction).
      local momi_y = 0.0 -- Ion momentum density (y-direction).
      local momi_z = (mass_ion / charge_ion) * Jz * Ti_frac -- Ion momentum density (z-direction).

      return rhoi, momi_x, momi_y, momi_z
    end,

    evolve = true, -- Evolve species?
    bcy = { G0.SpeciesBc.bcWall, G0.SpeciesBc.bcWall } -- Wall boundary conditions (y-direction).
  },

  -- Field.
  field = Moments.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,
    mgnErrorSpeedFactor = 1.0,

    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local Bxb = B0 * math.tanh(y / lambda) -- Total magnetic field strength.

      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = Bxb - psi0 * (pi / Ly) * math.cos(2.0 * pi * x / Lx) * math.sin(pi * y / Ly) -- Total magnetic field (x-direction).
      local By = psi0 * (2.0 * pi / Lx) * math.sin(2.0 * pi * x / Lx) * math.cos(pi * y / Ly) -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
    bcy = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall } -- Wall boundary conditions (y-direction).
  }
}

-- Run application.
momentApp:run()
