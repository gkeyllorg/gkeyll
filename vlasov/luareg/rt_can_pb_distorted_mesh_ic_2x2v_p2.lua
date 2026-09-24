local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

nl = 1.0 -- Left/inner number density.
Tl = 1.0 -- Left/inner temperature.
V_z1_drift_l = 0.0 -- Left/inner drift velocity (radial direction).
V_z2_drift_l = 0.0 -- Left/inner drift velocity (angular direction).

nr = 0.1 -- Right/outer number density.
Tr = 0.75 -- Right/outer temperature.
V_z1_drift_r = 0.0 -- Right/outer drift velocity (radial direction).
V_z2_drift_r = 0.0 -- Right/outer drift velocity (angular direction).

vt = 1.0 -- Thermal velocity.
nu = 15000.0 -- Collision frequency.

-- Simulation parameters.
Nz1 = 8 -- Cell count (configuration space: radial direction).
Nz2 = 8 -- Cell count (configuration space: angular direction).
Nvz1 = 8 -- Cell count (velocity space: radial direction).
Nvz2 = 8 -- Cell count (velocity space: angular direction).
Lz1 = 2.0 -- Domain size (configuration space: radial direction).
Lz2 = 2.0 -- Domain size (configuration space: angular direction).
vr_max = 12.0 * vt -- Domain boundary (velocity space: radial direction).
vtheta_max = 12.0 * vt -- Domain boundary (velocity space: angular direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 2e-3 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

midplane = 0.0 -- Radial midplane location designating jump in quantities.

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { - 1.0, - 1.0 },
  upper = { - 1.0 + Lz1, - 1.0 + Lz2 },
  cells = { Nz1, Nz2 },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 2 }, -- Periodic directions (angular direction only).

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.CanonicalPB,
    charge = charge, mass = mass,

    hamiltonian = function (t, xn)
      local q_z1, p_z1, p_z2 = xn[1], xn[3], xn[4]

      local inv_metric_z1_z1 = 1.0 -- Inverse metric tensor (z1-z1 component).
      local inv_metric_z1_z2 = 3.0 * pi * math.sin( pi * q_z1 ) / 2.0  -- Inverse metric tensor (z1-z2 component).
      local inv_metric_z2_z2 = ( 3.0 * pi * math.sin( pi * q_z1 ) / 2.0 ) * ( 3.0 * pi * math.sin( pi * q_z1 ) / 2.0 ) + 1.0  -- Inverse metric tensor (z2-z2 component).

      local hamiltonian = (0.5 * inv_metric_z1_z1 * p_z1 * p_z1) + (0.5 * (2.0 * inv_metric_z1_z2 * p_z1 * p_z2)) +
        (0.5 * inv_metric_z2_z2 * p_z2 * p_z2) -- Canonical Hamiltonian.

      return hamiltonian
    end,
    inverseMetric = function (t, xn)
      local q_z1 = xn[1]

      local inv_metric_z1_z1 = 1.0 -- Inverse metric tensor (z1-z1 component).
      local inv_metric_z1_z2 = 3.0 * pi * math.sin( pi * q_z1 ) / 2.0  -- Inverse metric tensor (z1-z2 component).
      local inv_metric_z2_z2 = ( 3.0 * pi * math.sin( pi * q_z1 ) / 2.0 ) * ( 3.0 * pi * math.sin( pi * q_z1 ) / 2.0 ) + 1.0  -- Inverse metric tensor (z2-z2 component).

      return inv_metric_z1_z1, inv_metric_z1_z2, inv_metric_z2_z2
    end,
    metric = function (t, xn)
      local q_z1 = xn[1]

      local metric_z1_z1 = 9.0 * pi * pi / 8.0 - 9.0 * pi * pi * math.cos( 2 * pi * q_z1 ) / 8.0 + 1.0 -- Metric tensor (z1-z1 component).
      local metric_z1_z2 = - 3.0 * pi * math.sin( pi * q_z1 ) / 2.0 -- Metric tensor (z1-z2 component).
      local metric_z2_z2 = 1.0 -- Metric tensor (z2-z2 component).

      return metric_z1_z1, metric_z1_z2, metric_z2_z2
    end,
    metricDeterminant = function (t, xn)

      local metric_det = 1 -- Metric tensor determinant.

      return metric_det
    end,
    
    -- Velocity space grid.
    lower = { -vr_max, -vtheta_max },
    upper = { vr_max, vtheta_max },
    cells = { Nvz1, Nvz2 },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local r = xn[1]

          local n = 0.0
          if r < midplane then
            n = nl -- Total number density (left/inner).
          else
            n = nr -- Total number density (right/outer).
          end

          local metric_det = 1

          return metric_det * n
        end,
        temperatureInit = function (t, xn)
          local r = xn[1]

          local T = 0.0
          if r < midplane then
            T = Tl -- Isotropic temperature (left/inner).
          else
            T = Tr -- Isotropic temperature (right/outer).
          end

          return T
        end,
        driftVelocityInit = function (t, xn)
          local r = xn[1]

          local V_z1_drift = 0.0
          local V_z2_drift = 0.0

          if r < midplane then
            V_z1_drift = V_z1_drift_l -- Radial drift velocity (left/inner).
            V_z2_drift = V_z2_drift_l -- Angular drift velocity (left/inner).
          else
            V_z1_drift = V_z1_drift_r -- Radial drift velocity (right/outer).
            V_z2_drift = V_z2_drift_r --Angular drift velocity (right/outer).
          end

          return V_z1_drift, V_z2_drift
        end,

        correctAllMoments = true,
        iterationEpsilon = 0.0,
        maxIterations = 0,
        useLastConverged = false
      }
    },

    collisions = {
      collisionID = G0.Collisions.BGK,

      selfNu = function (t, xn)
        return nu -- Collision frequency.
      end,
      
      useImplicitCollisionScheme = true,
      correctAllMoments = true,
      iterationEpsilon = 0.0,
      maxIterations = 0,
      useLastConverged = false
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcReflect
      },
      upper = {
        type = G0.SpeciesBc.bcReflect
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.LTEMoments, G0.Moment.EnergyMoment }
  },

  skipField = true,

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = 1.0, mu0 = 1.0,

    -- Initial conditions function.
    init = function (t, xn)
      return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0,

    isStatic = true
  }
}

vlasovApp:run()