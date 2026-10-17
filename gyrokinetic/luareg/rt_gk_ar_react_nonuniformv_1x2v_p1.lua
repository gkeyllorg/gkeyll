local Gyrokinetic = G0.Gyrokinetic

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using non-normalized physical units).
epsilon0 = 8.854187817620389850536563031710750260608e-12 -- Permittivity of free space.
mass_elc = 9.10938215e-31 -- Electron mass.
mass_ion = 2.014 * 1.672621637e-27 -- Proton mass.
mass_Ar1 = 39.95 * 1.672621637e-27-- Ar1+ ion mass.
mass_Ar2 = 39.95 * 1.672621637e-27 -- Ar2+ ion mass.
charge_elc = -1.602176487e-19 -- Electron charge.
charge_ion = 1.602176487e-19 -- Proton charge.
charge_Ar1 = 1.602176487e-19 -- Ar1+ ion charge.
charge_Ar2 = 2.0 * 1.602176487e-19 -- Ar2+ ion charge.

Te = 200.0 * 1.602176487e-19 -- Electron temperature.
Ti = 200.0 * 1.602176487e-19 -- Ion temperature.
TAr1 = 200.0 * 1.602176487e-19 -- Ar1+ temperature.
TAr2 = 200.0 * 1.602176487e-19 -- Ar2+ temperature.
n0_elc = 1.0e21 --  Electron reference number density (1 / m^3).

B_axis = 0.5 -- Magnetic field axis (simple toroidal coordinates).
R0 = 0.85 -- Major radius (simple toroidal coordinates).
a0 = 0.15 -- Minor axis (simple toroidal coordinates).

nu_frac = 0.25 -- Collision frequency fraction.

k_perp_rho_s = 0.3 -- Product of perpendicular wavenumber and ion-sound gyroradius.

-- Derived physical quantities (using non-normalized physical units).
n0_Ar1 = n0_elc * 1.0e-1 -- Ar1+ reference number density (1 / m^3).
n0_Ar2 = n0_elc * 1.0e-2 -- Ar2+ reference number density (1 / m^3).
n0_ion = n0_elc - n0_Ar1 - (2.0 * n0_Ar2) -- Ion reference number density (1 / m^3).

R = R0 + a0 -- Radial coordinate (simple toroidal coordinates).
B0 = B_axis * (R0 / R) -- Reference magnetic field strength (Tesla).

log_lambda_elc = 6.6 - 0.5 * math.log(n0_elc / 1.0e20) + 1.5 * math.log(Te / charge_ion) -- Electron Coulomb logarithm.
log_lambda_ion = 6.6 - 0.5 * math.log(n0_elc / 1.0e20) + 1.5 * math.log(Ti / charge_ion) -- Ion Coulomb logarithm.
nu_elc = nu_frac * log_lambda_elc * math.pow(charge_ion, 4.0) * n0_elc /
  (6.0 * math.sqrt(2.0) * math.pow(pi, 3.0 / 2.0) * math.pow(epsilon0, 2.0) * math.sqrt(mass_elc) * math.pow(Te, 3.0 / 2.0)) -- Electron collision frequency.
nu_ion = nu_frac * log_lambda_ion * math.pow(charge_ion, 4.0) * n0_elc /
  (12.0 * math.pow(pi, 3.0 / 2.0) * math.pow(epsilon0, 2.0) * math.sqrt(mass_ion) * math.pow(Ti, 3.0 / 2.0)) -- Ion collision frequency.

c_s = math.sqrt(Te / mass_ion) -- Sound speed.
vte = math.sqrt(Te / mass_elc) -- Electron thermal velocity.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.
vtAr1 = math.sqrt(TAr1 / mass_Ar1) -- Ar1+ thermal velocity.
vtAr2 = math.sqrt(TAr2 / mass_Ar2) -- Ar2+ thermal velocity.
omega_ci = math.abs(charge_ion * B0 / mass_ion) -- Ion cyclotron frequency.
rho_s = c_s / omega_ci -- Ion-sound gyroradius.

k_perp = k_perp_rho_s / rho_s -- Perpendicular wavenumber (for Poisson solver).

-- Simulation parameters.
Nz = 2 -- Cell count (configuration space: z-direction).
Nvpar = 6 -- Cell count (velocity space: parallel velocity direction).
Nmu = 4 -- Cell count (velocity space: magnetic moment direction).
Lz = 4.0 -- Domain size (configuration space: z-direction).
vpar_max_elc = 6.0 * vte -- Domain boundary (electron velocity space: parallel velocity direction).
mu_max_elc = (3.0 / 2.0) * 0.5 * mass_elc * math.pow(4.0 * vte, 2.0) / (2.0 * B0) -- Domain boundary (electron velocity space: magnetic moment direction).
vpar_max_ion = 6.0 * vti -- Domain boundary (ion velocity space: parallel velocity direction).
mu_max_ion = (3.0 / 2.0) * 0.5 * mass_ion * math.pow(4.0 * vti, 2.0) / (2.0 * B0) -- Domain boundary (ion velocity space: magnetic moment direction).
vpar_max_Ar1 = 6.0 * vtAr1 -- Domain boundary (Ar1+ velocity space: parallel velocity direction).
mu_max_Ar1 = (3.0 / 2.0) * 0.5 * mass_Ar1 * math.pow(4.0 * vtAr1, 2.0) / (2.0 * B0) -- Domain boundary (Ar1+ velocity space: magnetic moment direction).
vpar_max_Ar2 = 6.0 * vtAr2 -- Domain boundary (Ar2+ velocity space: parallel velocity direction).
mu_max_Ar2 = (3.0 / 2.0) * 0.5 * mass_Ar2 * math.pow(4.0 * vtAr2, 2.0) / (2.0 * B0) -- Domain boundary (Ar2+ velocity space: magnetic moment direction).
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
    lower = { -1.0, 0.0 },
    upper = { 1.0, 1.0 },
    cells = { Nvpar, Nmu },
    polarizationDensity = n0_elc,

    mapc2p = {
      -- Rescaled electron velocity space coordinates (vpar, mu) from old velocity space coordinates (cpvar, cmu).
      mapping = function (t, vc)
        local cvpar, cmu = vc[1], vc[2]

        local vpar = 0.0
        local mu = 0.0

        if cvpar < 0.0 then
          vpar = -vpar_max_elc * (cvpar * cvpar)
        else
          vpar = vpar_max_elc * (cvpar * cvpar)
        end
        mu = mu_max_elc * (cmu * cmu)

        return vpar, mu
      end
    },

    -- Initial conditions.
    projection = {
      projectionID = G0.Projection.MaxwellianPrimitive,

      densityInit = function (t, xn)
        return n0_elc -- Electron total number density.
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

      normalizeNu = true,
      referenceDensity = n0_elc,
      referenceTemperature = Te,

      selfNu = function (t, xn)
        return nu_elc
      end,

      numCrossCollisions = 3,
      collideWith = { "ion", "Ar1", "Ar2" }
    },

    reaction = {
      numReactions = 2,
      
      reactionTypes = {
        {
          reactionID = G0.Reaction.Ionization,
          selfType = G0.Self.Electron,
          ionType = G0.Ion.Argon,

          electronName = "elc",
          ionName = "Ar2",
          donorName = "Ar1",
          
          chargeState = 1,
          ionMass = mass_Ar2,
          electronMass = mass_elc
        },
        {
          reactionID = G0.Reaction.Recombination,
          selfType = G0.Self.Electron,
          ionType = G0.Ion.Argon,

          electronName = "elc",
          ionName = "Ar2",
          receiverName = "Ar1",
        
          chargeState = 1,
          ionMass = mass_Ar2,
          electronMass = mass_elc
        }
      },
    },

    bcs = {
      { dir = 0, edge = 0, type = G0.GyrokineticBc.speciesZeroFlux },
      { dir = 0, edge = 1, type = G0.GyrokineticBc.speciesZeroFlux },
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2par, G0.Moment.M2perp }
  },

  -- Ions.
  ion = Gyrokinetic.Species.new {
    balanceMoments = { G0.Moment.M0 },
    charge = charge_ion, mass = mass_ion,
    
    -- Velocity space grid.
    lower = { -vpar_max_ion, 0.0 },
    upper = { vpar_max_ion, mu_max_ion },
    cells = { Nvpar, Nmu },
    polarizationDensity = n0_ion,

    -- Initial conditions.
    projection = {
      projectionID = G0.Projection.MaxwellianPrimitive,

      densityInit = function (t, xn)
        return n0_ion -- Ion total number density.
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

      normalizeNu = true,
      referenceDensity = n0_elc,
      referenceTemperature = Ti,

      selfNu = function (t, xn)
        return nu_ion
      end,

      numCrossCollisions = 3,
      collideWith = { "elc", "Ar1", "Ar2" }
    },

    bcs = {
      { dir = 0, edge = 0, type = G0.GyrokineticBc.speciesZeroFlux },
      { dir = 0, edge = 1, type = G0.GyrokineticBc.speciesZeroFlux },
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2par, G0.Moment.M2perp }
  },

  -- Ar1+ ions.
  Ar1 = Gyrokinetic.Species.new {
    balanceMoments = { G0.Moment.M0 },
    charge = charge_Ar1, mass = mass_Ar1,
    
    -- Velocity space grid.
    lower = { -vpar_max_Ar1, 0.0 },
    upper = { vpar_max_Ar1, mu_max_Ar1 },
    cells = { Nvpar, Nmu },
    polarizationDensity = n0_Ar1,

    -- Initial conditions.
    projection = {
      projectionID = G0.Projection.MaxwellianPrimitive,

      densityInit = function (t, xn)
        return n0_Ar1 -- Ar1+ ion total number density.
      end,
      temperatureInit = function (t, xn)
        return TAr1 -- Ar1+ ion isotropic temperature.
      end,
      parallelVelocityInit = function (t, xn)
        return 0.0 -- Ar1+ ion parallel velocity.
      end
    },

    collisionless = {
      type = G0.GKCollisionless.GKCollisionlessES,
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      normalizeNu = true,
      referenceDensity = n0_Ar1,
      referenceTemperature = TAr1,

      selfNu = function (t, xn)
        return nu_ion
      end,

      numCrossCollisions = 3,
      collideWith = { "elc", "ion", "Ar2" }
    },

    reaction = {
      numReactions = 2,
      
      reactionTypes = {
        {
          reactionID = G0.Reaction.Ionization,
          selfType = G0.Self.Donor,
          ionType = G0.Ion.Argon,

          electronName = "elc",
          ionName = "Ar2",
          donorName = "Ar1",
          
          chargeState = 1,
          ionMass = mass_Ar2,
          electronMass = mass_elc
        },
        {
          reactionID = G0.Reaction.Recombination,
          selfType = G0.Self.Receiver,
          ionType = G0.Ion.Argon,

          electronName = "elc",
          ionName = "Ar2",
          receiverName = "Ar1",
        
          chargeState = 1,
          ionMass = mass_Ar2,
          electronMass = mass_elc
        }
      },
    },

    bcs = {
      { dir = 0, edge = 0, type = G0.GyrokineticBc.speciesZeroFlux },
      { dir = 0, edge = 1, type = G0.GyrokineticBc.speciesZeroFlux },
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2par, G0.Moment.M2perp }
  },

  -- Ar2+ ions.
  Ar2 = Gyrokinetic.Species.new {
    balanceMoments = { G0.Moment.M0 },
    charge = charge_Ar2, mass = mass_Ar2,
    
    -- Velocity space grid.
    lower = { -vpar_max_Ar2, 0.0 },
    upper = { vpar_max_Ar2, mu_max_Ar2 },
    cells = { Nvpar, Nmu },
    polarizationDensity = n0_Ar2,

    -- Initial conditions.
    projection = {
      projectionID = G0.Projection.MaxwellianPrimitive,

      densityInit = function (t, xn)
        return n0_Ar2 -- Ar2+ ion total number density.
      end,
      temperatureInit = function (t, xn)
        return TAr2 -- Ar2+ ion isotropic temperature.
      end,
      parallelVelocityInit = function (t, xn)
        return 0.0 -- Ar2+ ion parallel velocity.
      end
    },

    collisionless = {
      type = G0.GKCollisionless.GKCollisionlessES,
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      normalizeNu = true,
      referenceDensity = n0_Ar2,
      referenceTemperature = TAr2,

      selfNu = function (t, xn)
        return nu_ion
      end,

      numCrossCollisions = 3,
      collideWith = { "elc", "ion", "Ar1" }
    },

    reaction = {
      numReactions = 2,
      
      reactionTypes = {
        {
          reactionID = G0.Reaction.Ionization,
          selfType = G0.Self.Ion,
          ionType = G0.Ion.Argon,

          electronName = "elc",
          ionName = "Ar2",
          donorName = "Ar1",
          
          chargeState = 1,
          ionMass = mass_Ar2,
          electronMass = mass_elc
        },
        {
          reactionID = G0.Reaction.Recombination,
          selfType = G0.Self.Ion,
          ionType = G0.Ion.Argon,

          electronName = "elc",
          ionName = "Ar2",
          receiverName = "Ar1",
        
          chargeState = 1,
          ionMass = mass_Ar2,
          electronMass = mass_elc
        }
      },
    },

    bcs = {
      { dir = 0, edge = 0, type = G0.GyrokineticBc.speciesZeroFlux },
      { dir = 0, edge = 1, type = G0.GyrokineticBc.speciesZeroFlux },
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2par, G0.Moment.M2perp }
  },

  -- Field.
  field = Gyrokinetic.Field.new {
    femParBc = G0.ParProjBc.None,
    kPerpSq = k_perp * k_perp,

    zeroInitField = true, -- Don't compute the field at t = 0.
    isStatic = true -- Don't evolve the field in time.
  }
}

gyrokineticApp:run()
