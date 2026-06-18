local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

nl = 1.0 -- Left/inner number density.
Tl = 1.0 -- Left/inner temperature.
V_r_drift_l = 0.0 -- Left/inner drift velocity (radial direction).
V_theta_drift_l = 0.0 -- Left/inner drift velocity (polar angular direction).
V_phi_drift_l = 0.0 -- Left/inner drift velocity (azimuthal angular direction).

nr = 0.125 -- Right/outer number density.
Tr = math.sqrt(0.1 / 0.125) -- Right/outer temperature.
V_r_drift_r = 0.0 -- Right/outer drift velocity (radial direction).
V_theta_drift_r = 0.0 -- Right/outer drift velocity (polar angular direction).
V_phi_drift_r = 0.0 -- Right/outer drift velocity (azimuthal angular direction).

vt = 1.0 -- Thermal velocity.
nu = 15000.0 -- Collision frequency.

-- Simulation parameters.
Nr = 32 -- Cell count (configuration space: radial direction).
Ntheta = 1 -- Cell count (configuration space: polar angular direction).
Nphi = 1 -- Cell count (configuration space: azimuthal angular direction).
Nvr = 8 -- Cell count (velocity space: radial direction).
Nvtheta = 8 -- Cell count (velocity space: polar angular direction).
Nvphi = 8 -- Cell count (velocity space: azimuthal angular direction).
Lr = 1.0 -- Domain size (configuration space: radial direction).
Ltheta = 2.0 * pi -- Domain size (configuration space: polar angular direction).
Lphi = 2.0 * pi -- Domain size (configuration space: azimuthal angular direction).
vr_max = 8.0 * vt -- Domain boundary (velocity space: radial direction).
vtheta_max = 8.0 * vt -- Domain boundary (velocity space: polar angular direction).
vphi_max = 8.0 * vt -- Domain boundary (velocity space: azimuthal angular direction).
poly_order = 1 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.1 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

R = 2.0 -- Major radius of the torus.
midplane = 1.0 -- Radial midplane location designating jump in quantities.

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.5, 0.0, 0.0 },
  upper = { 0.5 + Lr, Ltheta, Lphi },
  cells = { Nr, Ntheta, Nphi },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 2, 3 }, -- Periodic directions (polar and azimuthal directions only).

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.CanonicalPB,
    charge = charge, mass = mass,

    hamiltonian = function (t, xn)
      local q_r, q_theta, p_r_dot, p_theta_dot, p_phi_dot = xn[1], xn[2], xn[4], xn[5], xn[6]

      local inv_metric_r_r = 1.0
      local inv_metric_r_theta = 0.0
      local inv_metric_r_phi = 0.0
      local inv_metric_theta_theta = 1.0 / (q_r * q_r)
      local inv_metric_theta_phi = 0.0
      local inv_metric_phi_phi = 1.0 / ((R + (q_r * math.cos(q_theta))) * (R + (q_r * math.cos(q_theta))))

      local hamiltonian = (0.5 * inv_metric_r_r * p_r_dot * p_r_dot) + (0.5 * (2.0 * inv_metric_r_theta * p_r_dot * p_theta_dot)) +
        (0.5 * (2.0 * inv_metric_r_phi * p_r_dot * p_phi_dot)) + (0.5 * inv_metric_theta_theta * p_theta_dot * p_theta_dot) +
        (0.5 * (2.0 * inv_metric_theta_phi * p_theta_dot * p_phi_dot)) + (0.5 * inv_metric_phi_phi * p_phi_dot * p_phi_dot) -- Canonical Hamiltonian.

      return hamiltonian
    end,
    metric = function (t, xn)
      local q_r, q_theta = xn[1], xn[2]

      local metric_r_r = 1.0 -- Metric tensor (radial-radial component).
      local metric_r_theta = 0.0 -- Metric tensor (radial-polar component).
      local metric_r_phi = 0.0 -- Metric tensor (radial-azimuthal component).
      local metric_theta_theta = q_r * q_r -- Metric tensor (polar-polar component).
      local metric_theta_phi = 0.0 -- Metric tensor (polar-azimuthal component).
      local metric_phi_phi = (R + (q_r * math.cos(q_theta))) * (R + (q_r * math.cos(q_theta))) -- Metric tensor (azimuthal-azimuthal component).

      return metric_r_r, metric_r_theta, metric_r_phi, metric_theta_theta, metric_theta_phi, metric_phi_phi
    end,
    inverseMetric = function (t, xn)
      local q_r, q_theta = xn[1], xn[2]

      local inv_metric_r_r = 1.0 -- Inverse metric tensor (radial-radial component).
      local inv_metric_r_theta = 0.0 -- Inverse metric tensor (radial-polar component).
      local inv_metric_r_phi = 0.0 -- Inverse metric tensor (radial-azimuthal component).
      local inv_metric_theta_theta = 1.0 / (q_r * q_r) -- Inverse metric tensor (polar-polar component).
      local inv_metric_theta_phi = 0.0 -- Inverse metric tensor (polar-azimuthal component).
      local inv_metric_phi_phi = 1.0 / ((R + (q_r * math.cos(q_theta))) * (R + (q_r * math.cos(q_theta)))) -- Inverse metric tensor (azimuthal-azimuthal component).

      return inv_metric_r_r, inv_metric_r_theta, inv_metric_r_phi, inv_metric_theta_theta, inv_metric_theta_phi, inv_metric_phi_phi
    end,
    metricDeterminant = function (t, xn)
      local q_r, q_theta = xn[1], xn[2]

      local metric_det = q_r * (R + (q_r * math.cos(q_theta))) -- Metric tensor determinant.

      return metric_det
    end,
    
    -- Velocity space grid.
    lower = { -vr_max, -vtheta_max, -vphi_max },
    upper = { vr_max, vtheta_max, vphi_max },
    cells = { Nvr, Nvtheta, Nvphi },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local r, theta = xn[1], xn[2]

          local n = 0.0
          if r < midplane then
            n = nl -- Total number density (left/inner).
          else
            n = nr -- Total number density (right/outer).
          end

          local metric_det = r * (R+ (r * math.cos(theta)))

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

          local V_r_drift = 0.0
          local V_theta_drift = 0.0
          local V_phi_drift = 0.0

          if r < midplane then
            V_r_drift = V_r_drift_l -- Radial drift velocity (left/inner).
            V_theta_drift = V_theta_drift_l -- Polar angular drift velocity (left/inner).
            V_phi_drift = V_phi_drift_l -- Azimuthal angular drift velocity (left/inner).
          else
            V_r_drift = V_r_drift_r -- Radial drift velocity (right/outer).
            V_theta_drift = V_theta_drift_r -- Polar angular drift velocity (right/outer).
            V_phi_drift = V_phi_drift_r -- Azimuthal angular drift velocity (right/outer).
          end

          return V_r_drift, V_theta_drift, V_phi_drift
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