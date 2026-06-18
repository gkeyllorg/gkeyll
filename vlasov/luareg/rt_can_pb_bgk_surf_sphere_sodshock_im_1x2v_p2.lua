local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

nl = 2.0 -- Left/inner number density.
Tl = 1.25 -- Left/inner temperature.
V_theta_drift_l = 0.0 -- Left/inner drift velocity (polar angular direction).
V_phi_drift_l = -0.5 -- Left/inner drift velocity (azimuthal angular direction).

nr = 0.25 -- Right/outer number density.
Tr = math.sqrt(0.1 / 0.125) -- Right/outer temperature.
V_theta_drift_r = 0.0 -- Right/outer drift velocity (polar angular direction).
V_phi_drift_r = 0.5 -- Right/outer drift velocity (azimuthal angular direction).

vt = 1.0 -- Thermal velocity.
nu = 15000.0 -- Collision frequency.

-- Simulation parameters.
Ntheta = 128 -- Cell count (configuration space: polar angular direction).
Nvtheta = 12 -- Cell count (velocity space: polar angular direction).
Nvphi = 12 -- Cell count (velocity space: azimuthal angular direction).
Ltheta = (3.0 * pi) / 8.0 -- Domain size (configuration space: polar angular direction).
vtheta_max = 8.0 * vt -- Domain boundary (velocity space: polar angular direction).
vphi_max = 8.0 * vt -- Domain boundary (velocity space: azimuthal angular direction).
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
midplane = pi / 2.0 -- Polar angular midplane location.
theta_loc = pi / 8.0 -- Polar angular boundary location designating jump in quantities.

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { pi / 8.0 },
  upper = { (pi / 8.0) + Ltheta },
  cells = { Ntheta },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.CanonicalPB,
    charge = charge, mass = mass,

    hamiltonian = function (t, xn)
      local q_theta, p_theta_dot, p_phi_dot = xn[1], xn[2], xn[3]

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

    outputfLTE = true,

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local theta = xn[1]

          local n = 0.0
          if math.abs(theta - midplane) < theta_loc then
            n = nl -- Total number density (left/inner).
          else
            n = nr -- Total number density (right/outer).
          end

          local metric_det = (R * R) * math.sin(theta)

          return metric_det * n
        end,
        temperatureInit = function (t, xn)
          local theta = xn[1]

          local T = 0.0
          if math.abs(theta - midplane) < theta_loc then
            T = Tl -- Isotropic temperature (left/inner).
          else
            T = Tr -- Isotropic temperature (right/outer).
          end

          return T
        end,
        driftVelocityInit = function (t, xn)
          local theta = xn[1]

          local V_theta_drift = 0.0
          local V_phi_drift = 0.0

          if math.abs(theta - midplane) < theta_loc then
            V_theta_drift = V_theta_drift_l -- Polar angular drift velocity (left/inner).
            V_phi_drift = V_phi_drift_l -- Azimuthal angular drift velocity (left/inner).
          else
            V_theta_drift = V_theta_drift_r -- Polar angular drift velocity (right/outer).
            V_phi_drift = V_phi_drift_r -- Azimuthal angular drift velocity (right/outer).
          end

          return V_theta_drift, V_phi_drift
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