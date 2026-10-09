local Gyrokinetic = G0.Gyrokinetic

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using non-normalized physical units).
epsilon0 = 8.854187817620389850536563031710750260608e-12 -- Permittivity of free space.
mass_elc = 9.10938215e-31 -- Electron mass.
mass_ion = 1.672621637e-27 -- Proton mass.
charge_elc = -1.602176487e-19 -- Electron charge.
charge_ion = 1.602176487e-19 -- Proton charge.

Te = 30.0 * 1.602176487e-19 -- Electron temperature.
Ti = 30.0 * 1.602176487e-19 -- Ion temperature.
B0 = 1.0 -- Reference magnetic field strength (Tesla).
n0 = 1.0e19 --  Reference number density (1 / m^3).

nu_frac = 0.1 -- Collision frequency fraction.

k_perp_rho_s = 0.1 -- Product of perpendicular wavenumber and ion-sound gyroradius.

-- Derived physical constants (using non-normalized physical units).
log_lambda_elc = 6.6 - 0.5 * math.log(n0 / 1.0e20) + 1.5 * math.log(Te / charge_ion) -- Electron Coulomb logarithm.
log_lambda_ion = 6.6 - 0.5 * math.log(n0 / 1.0e20) + 1.5 * math.log(Ti / charge_ion) -- Ion Coulomb logarithm.
nu_elc = nu_frac * log_lambda_elc * math.pow(charge_ion, 4.0) * n0 /
  (6.0 * math.sqrt(2.0) * math.pow(pi, 3.0 / 2.0) * math.pow(epsilon0, 2.0) * math.sqrt(mass_elc) * math.pow(Te, 3.0 / 2.0)) -- Electron collision frequency.
nu_ion = nu_frac * log_lambda_ion * math.pow(charge_ion, 4.0) * n0 /
  (12.0 * math.pow(pi, 3.0 / 2.0) * math.pow(epsilon0, 2.0) * math.sqrt(mass_ion) * math.pow(Ti, 3.0 / 2.0)) -- Ion collision frequency.

c_s = math.sqrt(Te / mass_ion) -- Sound speed.
vte = math.sqrt(Te / mass_elc) -- Electron thermal velocity.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.
omega_ci = math.abs(charge_ion * B0 / mass_ion) -- Ion cyclotron frequency.
rho_s = c_s / omega_ci -- Ion-sound gyroradius.

k_perp = k_perp_rho_s / rho_s -- Perpendicular wavenumber (for Poisson solver).

-- Simulation parameters.
Nz = 2 -- Cell count (configuration space: z-direction).
Nvpar = 16 -- Cell count (velocity space: parallel velocity direction).
Nmu = 8 -- Cell count (velocity space: magnetic moment direction).
Lz = 100.0 * rho_s -- Domain size (configuration space: z-direction).
vpar_max_elc = 4.0 * vte -- Domain boundary (electron velocity space: parallel velocity direction).
mu_max_elc = 0.75 * mass_elc * math.pow(4.0 * vte, 2.0) / (2.0 * B0) -- Domain boundary (electron velocity space: magnetic moment direction).
vpar_max_ion = 4.0 * vti -- Domain boundary (ion velocity space: parallel velocity direction).
mu_max_ion = 0.75 * mass_ion * math.pow(4.0 * vti, 2.0) / (2.0 * B0) -- Domain boundary (ion velocity space: magnetic moment direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 1.0e-7 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

gyrokineticApp = Gyrokinetic.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { -0.5 * Lz },
  upper = { 0.5 * Lz },
  cells = { Nz },
  cflFrac = cfl_frac,

  --basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  geometry = {
    geometryID = G0.Geometry.MapC2P,
    world = { 0.0, 0.0, 0.0 },

    -- Computational coordinates (x, y, z) from physical coordinates (X, Y, Z).
    mapc2p = function (t, zc)
      local xp = { }
      
      xp[1] = zc[1]
      xp[2] = zc[2]
      xp[3] = zc[3]

      return xp[1], xp[2], xp[3]
    end,

    -- Magnetic field strength.
    bfieldFunc = function (t, zc)
      return 0.0, 0.0, B0
    end
  },

  -- Electrons.
  elc = Gyrokinetic.Species.new {
    balanceMoments = { G0.Moment.M0 },
    charge = charge_elc, mass = mass_elc,
    
    -- Velocity space grid.
    lower = { -vpar_max_elc, 0.0 },
    upper = { vpar_max_elc, mu_max_elc },
    cells = { Nvpar, Nmu },
    polarizationDensity = n0,

    -- Initial conditions.
    projection = {
      projectionID = G0.Projection.MaxwellianPrimitive,

      densityInit = function (t, xn)
        return n0 -- Electron total number density.
      end,
      temperatureInit = function (t, xn)
        return Te -- Electron isotropic temperature.
      end,
      parallelVelocityInit = function (t, xn)
        return 0.0 -- Electron parallel velocity.
      end
    },

    collisionless = {
      type = G0.GKCollisionless.GKCollisionlessES,
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_elc
      end,

      numCrossCollisions = 1,
      collideWith = { "ion" }
    },

    radiation = {
      radiationID = G0.Radiation.GKRadiation,

      numCrossCollisions = 1,
      collideWith = { "ion" },

      atomicZ = { 1 },
      chargeState = { 0 },
      numDensities = { 1 },

      TeMinModel = G0.TeMinModel.Const,
      TeMin = 12.0 * 1.602176487e-19
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2par, G0.Moment.M2perp, G0.Moment.M3par, G0.Moment.M3perp }
  },

  -- Ions.
  ion = Gyrokinetic.Species.new {
    balanceMoments = { G0.Moment.M0 },
    charge = charge_ion, mass = mass_ion,
    
    -- Velocity space grid.
    lower = { -vpar_max_ion, 0.0 },
    upper = { vpar_max_ion, mu_max_ion },
    cells = { Nvpar, Nmu },
    polarizationDensity = n0,

    -- Initial conditions.
    projection = {
      projectionID = G0.Projection.MaxwellianPrimitive,

      densityInit = function (t, xn)
        return n0 -- Ion total number density.
      end,
      temperatureInit = function (t, xn)
        return Ti -- Ion isotropic temperature.
      end,
      parallelVelocityInit = function (t, xn)
        return 0.0 -- Ion parallel velocity.
      end
    },

    collisionless = {
      type = G0.GKCollisionless.GKCollisionlessES,
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_ion
      end,

      numCrossCollisions = 1,
      collideWith = { "elc" }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2par, G0.Moment.M2perp, G0.Moment.M3par, G0.Moment.M3perp }
  },

  -- Field.
  field = Gyrokinetic.Field.new {
    femParBc = G0.ParProjBc.Periodic,
    kPerpSq = k_perp * k_perp,

    zeroInitField = true, -- Don't compute the field at t = 0.
    isStatic = true -- Don't evolve the field in time.
  }
}

gyrokineticApp:run()
