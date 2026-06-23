local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

n0 = 1.0 -- Reference number density.
T0 = 1.0 -- Reference temperature.
Vx_drift = 0.0 -- Drift velocity (x-direction).
Vy_drift = 0.0 -- Drift velocity (y-direction).

vt = 1.0 -- Thermal velocity.
nu = 15000.0 -- Collision frequency.

-- Simulation parameters.
Nx = 2 -- Cell count (configuration space: x-direction).
Ny = 2 -- Cell count (configuration space: y-direction).
Nvx = 8 -- Cell count (velocity space: vx-direction).
Nvy = 8 -- Cell count (velocity space: vy-direction).
Lx = 1.0 -- Domain size (configuration space: x-direction).
Ly = 1.0 -- Domain size (configuration space: y-direction).
vx_max = 6.0 * vt -- Domain boundary (velocity space: vx-direction).
vy_max = 6.0 * vt -- Domain boundary (velocity space: vy-direction).
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

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.0, 0.0 },
  upper = { Lx, Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1, 2 }, -- Periodic directions (x- and y-directions only).

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.CanonicalPB,
    charge = charge, mass = mass,

    hamiltonian = function (t, xn)
      local p_x_dot, p_y_dot = xn[3], xn[4]

      local inv_metric_x_x = 1.0
      local inv_metric_x_y = 0.0
      local inv_metric_y_y = 1.0
    
      local hamiltonian = (0.5 * inv_metric_x_x * p_x_dot * p_x_dot) + (0.5 * (2.0 * inv_metric_x_y * p_x_dot * p_y_dot)) +
        (0.5 * inv_metric_y_y * p_y_dot * p_y_dot) -- Canonical Hamiltonian.

      return hamiltonian
    end,
    inverseMetric = function (t, xn)
      local inv_metric_x_x = 1.0 -- Inverse metric tensor (x-x component).
      local inv_metric_x_y = 0.0 -- Inverse metric tensor (x-y component).
      local inv_metric_y_y = 1.0 -- Inverse metric tensor (y-y component).

      return inv_metric_x_x, inv_metric_x_y, inv_metric_y_y
    end,
    metric = function (t, xn)
      local metric_x_x = 1.0 -- Metric tensor (x-x component).
      local metric_x_y = 0.0 -- Metric tensor (x-y component).
      local metric_y_y = 1.0 -- Metric tensor (y-y component).

      return metric_x_x, metric_x_y, metric_y_y
    end,
    metricDeterminant = function (t, xn)
      local metric_det = 1.0 -- Metric tensor determinant.

      return metric_det
    end,
    
    -- Velocity space grid.
    lower = { -vx_max, -vy_max },
    upper = { vx_max, vy_max },
    cells = { Nvx, Nvy },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local metric_det = 1.0

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
      
      useImplicitCollisionScheme = false,
      correctAllMoments = true,
      iterationEpsilon = 0.0,
      maxIterations = 0,
      useLastConverged = false
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