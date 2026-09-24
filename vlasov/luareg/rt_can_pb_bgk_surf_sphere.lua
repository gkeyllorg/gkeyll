local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

n0 = 1.0 -- Reference number density.
T0 = 1.0 -- Reference temperature.
V_theta_drift = 0.0 -- Drift velocity (polar angular direction).
V_phi_drift = 0.0 -- Drift velocity (azimuthal angular direction).

vt = 1.0 -- Thermal velocity.
nu = 100.0 -- Collision frequency.

-- Simulation parameters.
Ntheta = 8 -- Cell count (configuration space: polar angular direction).
Nphi = 32 -- Cell count (configuration space: azimuthal angular direction).
Nvtheta = 8 -- Cell count (velocity space: polar angular direction).
Nvphi = 8 -- Cell count (velocity space: azimuthal angular direction).
Ltheta = pi / 2.0 -- Domain size (configuration space: polar angular direction).
Lphi = 2.0 * pi -- Domain size (configuration space: azimuthal angular direction).
vtheta_max = 5.0 * vt -- Domain boundary (velocity space: polar angular direction).
vphi_max = 5.0 * vt -- Domain boundary (velocity space: azimuthal angular direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.1 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

R = 1.0 -- Radius of the sphere.

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { pi / 4.0, 0.0 },
  upper = { (pi / 4.0) + Ltheta, Lphi },
  cells = { Ntheta, Nphi },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 2 }, -- Periodic directions (azimuthal angular direction only).

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.CanonicalPB,
    charge = charge, mass = mass,

    hamiltonian = function (t, xn)
      local q_theta, p_theta_dot, p_phi_dot = xn[1], xn[3], xn[4]

      local inv_metric_theta_theta = 1.0 / (R * R)
      local inv_metric_theta_phi = 0.0
      local inv_metric_phi_phi = 1.0 / ((R * math.sin(q_theta)) * (R * math.sin(q_theta)))
    
      local hamiltonian = (0.5 * inv_metric_theta_theta * p_theta_dot * p_theta_dot) + (0.5 * (2.0 * inv_metric_theta_phi * p_theta_dot * p_phi_dot)) +
        (0.5 * inv_metric_phi_phi * p_phi_dot * p_phi_dot) -- Canonical Hamiltonian.

      return hamiltonian
    end,
    inverseMetric = function (t, xn)
      local q_theta = xn[1]

      local inv_metric_theta_theta = 1.0 / (R * R) -- Inverse metric tensor (polar-polar component).
      local inv_metric_theta_phi = 0.0 -- Inverse metric tensor (polar-azimuthal component).
      local inv_metric_phi_phi = 1.0 / ((R * math.sin(q_theta)) * (R * math.sin(q_theta))) -- Inverse metric tensor (azimuthal-azimuthal component).

      return inv_metric_theta_theta, inv_metric_theta_phi, inv_metric_phi_phi
    end,
    metric = function (t, xn)
      local q_theta = xn[1]

      local metric_theta_theta = R * R -- Metric tensor (polar-polar component).
      local metric_theta_phi = 0.0 -- Metric tensor (polar-azimuthal component).
      local metric_phi_phi = (R * math.sin(q_theta)) * (R * math.sin(q_theta)) -- Metric tensor (azimuthal-azimuthal component).

      return metric_theta_theta, metric_theta_phi, metric_phi_phi
    end,
    metricDeterminant = function (t, xn)
      local q_theta = xn[1]

      local metric_det = (R * R) * math.sin(q_theta) -- Metric tensor determinant.

      return metric_det
    end,
    
    -- Velocity space grid.
    lower = { -vtheta_max, -vphi_max },
    upper = { vtheta_max, vphi_max },
    cells = { Nvtheta, Nvphi },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local theta = xn[1]
          
          local metric_det = (R * R) * math.sin(theta)

          return metric_det * n0 -- Total number density.
        end,
        temperatureInit = function (t, xn)
          return T0 -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return V_theta_drift, V_phi_drift -- Total drift velocity.
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
      
      useImplicitCollisionScheme = true
    },

    correct = {
      outputfLTE = true,
      correctAllMoments = true,
      iterationEpsilon = 1e-12,
      maxIterations = 100,
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
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.LTEMoments }
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