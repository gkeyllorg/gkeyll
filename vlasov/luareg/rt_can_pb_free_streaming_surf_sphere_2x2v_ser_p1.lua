-- Collisionless free streaming of a Maxwellian gas on the unit sphere (2x2v). A density pattern n =
-- 0.3 + 2 sin^4 theta sin^4(3 phi / 2) at temperature 1 streams along great circles between
-- reflecting walls at theta = pi/4 and 3 pi/4 (periodic in phi); the reference density at t = 0.5
-- is the exact solution by characteristics, n(x, t) = int n_0(x_-t(x, v)) M(v) d^2v with the
-- great-circle transport of the velocity. Figures of merit (serendipity p1, 8 x 16 x 8 x 8 cells):
-- at t = 0.5 the density differs from the characteristics solution by at most 0.26 (L1 0.065) while
-- the pattern has decayed from the range 0.3-2.26 to 0.66-1.06; mass and energy conserved to
-- round-off.
local Vlasov = G0.Vlasov

mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.
n0 = 0.3 -- Background number density.
T0 = 1.0 -- Temperature.
vt = 1.0 -- Thermal velocity.

Ntheta = 8 -- Cell count (configuration space: polar direction).
Nphi = 16 -- Cell count (configuration space: azimuthal direction).
Nvtheta = 8 -- Cell count (velocity space: polar momentum).
Nvphi = 8 -- Cell count (velocity space: azimuthal momentum).
Ltheta = math.pi / 2.0 -- Domain size (configuration space: polar direction).
Lphi = 2.0 * math.pi -- Domain size (configuration space: azimuthal direction).
vtheta_max = 5.0 * vt -- Domain boundary (velocity space: polar momentum).
vphi_max = 5.0 * vt -- Domain boundary (velocity space: azimuthal momentum).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.5 -- Final simulation time.
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
  lower = { math.pi / 4.0, 0.0 },
  upper = { (math.pi / 4.0) + Ltheta, Lphi },
  cells = { Ntheta, Nphi },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = { 2 }, -- Periodic directions.

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.CanonicalPB,
    charge = charge, mass = mass,
    hamiltonian = function (t, xn)
      local q_theta = xn[1]
      local p_theta_dot, p_phi_dot = xn[3], xn[4]

      local inv_metric_aa = 1.0
      local inv_metric_ab = 0.0
      local inv_metric_bb = 1.0 / (math.sin(q_theta) * math.sin(q_theta))

      -- Canonical Hamiltonian H = (1/2) g^ij p_i p_j.
      return (0.5 * inv_metric_aa * p_theta_dot * p_theta_dot) +
        (inv_metric_ab * p_theta_dot * p_phi_dot) +
        (0.5 * inv_metric_bb * p_phi_dot * p_phi_dot)
    end,
    inverseMetric = function (t, xn)
      local q_theta = xn[1]
      -- Inverse metric tensor (aa, ab, bb components).
      return 1.0, 0.0, 1.0 / (math.sin(q_theta) * math.sin(q_theta))
    end,
    metric = function (t, xn)
      local q_theta = xn[1]
      -- Metric tensor (aa, ab, bb components).
      return 1.0, 0.0, math.sin(q_theta) * math.sin(q_theta)
    end,
    metricDeterminant = function (t, xn)
      local q_theta = xn[1]
      -- Metric tensor determinant (square root of det g).
      return math.sin(q_theta)
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
          local theta, phi = xn[1], xn[2]
          local s4 = math.pow(math.sin(theta), 4.0)

          -- Density pattern n0 + 2 sin^4(theta) sin^4(3 phi / 2), times the metric determinant.
          return math.sin(theta) * (n0 + (2.0 * s4 * math.pow(math.sin(1.5 * phi), 4.0)))
        end,
        temperatureInit = function (t, xn)
          -- Isotropic temperature.
          return T0
        end,
        driftVelocityInit = function (t, xn)
          -- Total drift velocity (no flow).
          return 0.0, 0.0
        end,
        correctAllMoments = true,
        iterationEpsilon = 0.0,
        maxIterations = 0,
        useLastConverged = false
      }
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
}

vlasovApp:run()
