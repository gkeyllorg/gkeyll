local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

nl = 1.0 -- Left number density.
Tl = 1.0 -- Left temperature.
Vx_drift_l = 0.0 -- Left drift velocity (x-direction).
Vy_drift_l = 0.0 -- Left drift velocity (y-direction).

nr = 0.125 -- Right number density.
Tr = math.sqrt(0.1 / 0.125) -- Right temperature.
Vx_drift_r = 0.0 -- Right drift velocity (x-direction).
Vy_drift_r = 0.0 -- Right drift velocity (y-direction).

vt = 1.0 -- Thermal velocity.
nu = 15000.0 -- Collision frequency.

-- Simulation parameters.
Nx = 128 -- Cell count (configuration space: x-direction).
Nvx = 16 -- Cell count (velocity space: vx-direction).
Nvy = 16 -- Cell count (velocity space: vy-direction).
Lx = 1.0 -- Domain size (configuration space: x-direction).
vx_max = 6.0 * vt -- Domain boundary (velocity space: vx-direction).
vy_max = 6.0 * vt -- Domain boundary (velocity space: vy-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

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
  lower = { 0.0 },
  upper = { Lx },
  cells = { Nx },
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
      local p_x_dot, p_y_dot = xn[2], xn[3]

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

    outputfLTE = true,

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]

          local n = 0.0
          if x < 0.5 then
            n = nl -- Total number density (left).
          else
            n = nr -- Total number density (right).
          end

          local metric_det = 1.0

          return metric_det * n
        end,
        temperatureInit = function (t, xn)
          local x = xn[1]

          local T = 0.0
          if x < 0.5 then
            T = Tl -- Isotropic temperature (left).
          else
            T = Tr -- Isotropic temperature (right).
          end

          return T
        end,
        driftVelocityInit = function (t, xn)
          local x = xn[1]

          local Vx_drift = 0.0
          local Vy_drift = 0.0

          if x < 0.5 then
            Vx_drift = Vx_drift_l -- Left drift velocity (x-direction).
            Vy_drift = Vy_drift_l -- Left drift velocity (y-direction).
          else
            Vx_drift = Vx_drift_r -- Right drift velocity (x-direction).
            Vy_drift = Vy_drift_r -- Right drift velocity (y-direction).
          end

          return Vx_drift, Vy_drift
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