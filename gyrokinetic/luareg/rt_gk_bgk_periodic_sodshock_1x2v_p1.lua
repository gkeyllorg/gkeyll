local Gyrokinetic = G0.Gyrokinetic

-- Physical constants (using normalized code units).
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

nl = 1.0 -- Left number density.
Tl = 1.0 -- Left temperature.

nr = 0.125 -- Right number density.
Tr = math.sqrt(0.1 / 0.125) -- Right temperature.

B0 = 1.0 -- Reference magnetic field strength.
n0 = 1.0 -- Reference number density.
vt = 1.0 -- Reference thermal velocity. 
nu = 100.0 -- Collision frequency.

-- Simulation parameters.
Nz = 64 -- Cell count (configuration space: z-direction).
Nvpar = 16 -- Cell count (velocity space: parallel velocity direction).
Nmu = 16 -- Cell count (velocity space: magnetic moment direction).
Lz = 2.0 -- Domain size (configuration space: z-direction).
vpar_max = 6.0 * vt -- Domain boundary (velocity space: parallel velocity direction).
mu_max = 18.0 * (vt * vt) / 2.0 / B0 -- Domain boundary (velocity space: magnetic moment direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.1 -- Final simulation time.
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

  -- Neutral species.
  neut = Gyrokinetic.Species.new {
    balanceMoments = { G0.Moment.M0 },
    charge = charge, mass = mass,
    
    -- Velocity space grid.
    lower = { -vpar_max, 0.0 },
    upper = { vpar_max, mu_max },
    cells = { Nvpar, Nmu },
    polarizationDensity = n0,

    -- Initial conditions.
    projection = {
      projectionID = G0.Projection.MaxwellianPrimitive,

      densityInit = function (t, xn)
        local z = xn[1]

        local n = 0.0

        if math.abs(z) < 0.5 then
          n = nl -- Total number density (left).
        else
          n = nr -- Total number density (right).
        end
        
        return n
      end,
      temperatureInit = function (t, xn)
        local z = xn[1]

        local T = 0.0

        if math.abs(z) < 0.5 then
          T = Tl -- Isotropic temperature (left).
        else
          T = Tr -- Isotropic temperature (right).
        end

        return T
      end,
      parallelVelocityInit = function (t, xn)
        return 0.0 -- Parallel velocity.
      end,

      correctAllMoments = true
    },

    correct = {
      correctAllMoments = true,
      iterationEpsilon = 1.0e-12,
      maxIterations = 10,
      useLastConverged = true
    }, 

    collisionless = {
      type = G0.GKCollisionless.GKCollisionlessES,
    },

    collisions = {
      collisionID = G0.Collisions.BGK,

      selfNu = function (t, xn)
        return nu
      end
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2par, G0.Moment.M2perp, G0.Moment.MaxwellianMoments }
  },
  
  -- Field.
  field = Gyrokinetic.Field.new {
    fieldID = G0.GKField.Boltzmann,

    electronMass = mass,
    electronCharge = charge,
    electronTemperature = vt,
    femParBc = G0.ParProjBc.None,

    zeroInitField = true, -- Don't compute the field at t = 0.
    isStatic = true -- Don't evolve the field in time.
  }
}

gyrokineticApp:run()
