local Gyrokinetic = G0.Gyrokinetic

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
mass = 1.0 -- Top hat/bump mass.
charge = 1.0 -- Top hat/bump charge.

B0 = 1.0 -- Reference magnetic field strength.
n0 = 1.0 -- Reference number density.
u0 = 0.0 -- Reference velocity.
vt = 1.0 / 3.0 -- Top hat Maxwellian thermal velocity.
nu = 0.01 -- Collision frequency.

ab = math.sqrt(0.1) -- Bump Maxwellian amplitude.
sb = 0.12 -- Bump Maxwellian softening factor, to avoid divergence.
vtb = 1.0 -- Bump Maxwellian thermal velocity.

-- Derived physical quantities (using normalized code units).
ub = 4.0 * math.sqrt((math.pow(3.0 * vt / 2.0, 2.0)) / 3.0) -- Bump location (in velocity space).

-- Simulation parameters.
Nz = 2 -- Cell count (configuration space: z-direction).
Nvpar = 32 -- Cell count (velocity space: parallel velocity direction).
Nmu = 16 -- Cell count (velocity space: magnetic moment direction).
Lz = 1.0 -- Domain size (configuration space: x-direction).
vpar_max = 8.0 * vt -- Domain boundary (velocity space: parallel velocity direction).
mu_max = 0.5 * mass * math.pow(3.5 * vt, 2.0) / B0 -- Domain boundary (velocity space: magnetic moment direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.5 / nu -- Final simulation time.
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
  lower = { 0.0 },
  upper = { Lz },
  cells = { Nz },
  cflFrac = cfl_frac,

  basis = basis_type,
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

    -- Top hat species.
  square = Gyrokinetic.Species.new {
    charge = charge, mass = mass,
    
    -- Velocity space grid.
    lower = { -vpar_max, 0.0 },
    upper = { vpar_max, mu_max },
    cells = { Nvpar, Nmu },
    polarizationDensity = n0,

    -- Initial conditions.
    projection = {
      projectionID = G0.Projection.Func,

      init = function (t, xn)
        local vpar = xn[2]

        local v0 = math.sqrt(3.0) * vt

        local n = 0.0;
      
        if math.abs(vpar) < v0 then
          n = n0 / 2.0 / v0 -- Distribution function (low velocity).
        else
          n = 0.0 -- Distribution function (high velocity).
        end
        
        return n
      end,
    },

    collisionless = {
      type = G0.GKCollisionless.GKCollisionlessES,
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu
      end
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2par, G0.Moment.M2perp, G0.Moment.M3par, G0.Moment.M3perp }
  },

  -- Bump species.
  bump = Gyrokinetic.Species.new {
    charge = charge, mass = mass,
    
    -- Velocity space grid.
    lower = { -vpar_max, 0.0 },
    upper = { vpar_max, mu_max },
    cells = { Nvpar, Nmu },
    polarizationDensity = n0,

    -- Initial conditions.
    projection = {
      projectionID = G0.Projection.Func,

      init = function (t, xn)
        local vpar, mu = xn[2], xn[3]

        local v_sq = ((vpar - u0) / (math.sqrt(2.0) * vt)) * ((vpar - u0) / (math.sqrt(2.0) * vt)) + mu * B0
        local vb_sq = ((vpar - u0) / (math.sqrt(2.0) * vtb)) * ((vpar - u0) / (math.sqrt(2.0) * vtb)) + mu * B0
        
        local n = (n0 / math.sqrt(2.0 * pi * vt)) * math.exp(-v_sq) + (n0 / math.sqrt(2.0 * pi * vtb)) *
          math.exp(-vb_sq) * (ab * ab) / ((vpar - ub) * (vpar - ub) + sb * sb) -- Distribution function.
        
        return n
      end,
    },

    collisionless = {
      type = G0.GKCollisionless.GKCollisionlessES,
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu
      end
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2par, G0.Moment.M2perp, G0.Moment.M3par, G0.Moment.M3perp }
  },
  
  -- Field.
  field = Gyrokinetic.Field.new {
    fieldID = G0.GKField.Boltzmann,

    electronMass = mass,
    electronCharge = charge,
    electronTemperature = vt * vt * mass,
    femParBc = G0.ParProjBc.None,
    
    zeroInitField = true, -- Don't compute the field at t = 0.
    isStatic = true -- Don't evolve the field in time.
  }
}

gyrokineticApp:run()
