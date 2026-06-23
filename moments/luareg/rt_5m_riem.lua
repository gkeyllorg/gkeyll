-- Generalized Brio-Wu Riemann problem for the 5-moment equations.
-- Input parameters match the initial conditions found in entry JE4 of Ammar's Simulation Journal (https://ammar-hakim.org/sj/je/je4/je4-twofluid-shock.html), adapted from Section 7.1 of the article:
-- A. Hakim, J. Loverich and U. Shumlak (2006), "A high resolution wave propagation scheme for ideal Two-Fluid plasma equations",
-- Journal of Computational Physics, Volume 219 (1): 418-442.
-- https://www.sciencedirect.com/science/article/pii/S0021999106001707

local Moments = G0.Moments
local Euler = G0.Moments.Eq.Euler

-- Physical constants (using normalized code units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_ion = 1.0 -- Proton mass.
charge_ion = 1.0 -- Proton charge.
mass_elc = 1.0 / 1836.2 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

rhol_ion = 1.0 -- Left ion mass density.
rhor_ion = 0.125 -- Right ion mass density
pl = 5.0e-5 -- Left electron/ion pressure.
pr = 5.0e-6 -- Right electron/ion pressure.

Bx = 0.75e-2 -- Total magnetic field (x-direction).
Bzl = 1.0e-2 -- Left total magneic field (z-direction).
Bzr = -1.0e-2 -- Right total magnetic field (z-direction).

has_collision = false -- Whether to include collisions.
nu_base_ei = 0.5 -- Base electron-ion collision frequency.

-- Derived physical quantities (using normalized code units).
rhol_elc = rhol_ion * mass_elc / mass_ion -- Left electron mass density.
rhor_elc = rhor_ion * mass_elc / mass_ion -- Right electron mass density.

-- Simulation parameters.
Nx = 1024 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
cfl_frac = 0.95 -- CFL coefficient.

t_end = 10.0 -- Final simulation time.
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

  hasCollision = has_collision,
  nuBase = {
    { 0.0, nu_base_ei },
    { nu_base_ei, 0.0 }
  },

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

      local rho = 0.0
      local p = 0.0
    
      if x < 0.5 then
        rho = rhol_elc -- Electron mass density (left).
        p = pl -- Electron pressure (left).
      else
        rho = rhor_elc -- Electron mass density (right).
        p = pr -- Electron pressure (right).
      end
    
      local rhoe = rho -- Electron mass density.
      local mome_x = 0.0 -- Electron momentum density (x-direction).
      local mome_y = 0.0 -- Electron momentum density (y-direction).
      local mome_z = 0.0 -- Electron momentum density (z-direction).
      local Ee_tot = p / (gas_gamma - 1.0) -- Electron total energy density.
    
      return rhoe, mome_x, mome_y, mome_z, Ee_tot
    end,

    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (x-direction).
  },

  -- Ions.
  ion = Moments.Species.new {
    charge = charge_ion, mass = mass_ion,
    equation = Euler.new { gasGamma = gas_gamma },

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local rho = 0.0
      local p = 0.0
    
      if x < 0.5 then
        rho = rhol_ion -- Ion mass density (left).
        p = pl -- Ion pressure (left).
      else
        rho = rhor_ion -- Ion mass density (right).
        p = pr -- Ion pressure (right).
      end
    
      local rhoi = rho -- Ion mass density.
      local momi_x = 0.0 -- Ion momentum density (x-direction).
      local momi_y = 0.0 -- Ion momentum density (y-direction).
      local momi_z = 0.0 -- Ion momentum density (z-direction).
      local Ei_tot = p / (gas_gamma - 1.0) -- Ion total energy density.
    
      return rhoi, momi_x, momi_y, momi_z, Ei_tot
    end,

    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (x-direction).
  },

  -- Field.
  field = Moments.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).
    
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0
    
      if x < 0.5 then
        Bz = Bzl -- Total magnetic field (z-direction, left).
      else
        Bz = Bzr -- Total magnetic field (z-direction, right).
      end

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
    bcx = { G0.FieldBc.bcCopy, G0.FieldBc.bcCopy } -- Copy boundary conditions (x-direction).
  }
}

-- Run application.
momentApp:run()
