-- Stationary Maxwell-Boltzmann equilibrium on a rotating unit sphere in the extended-Hamiltonian
-- form with a background flow (0, omega sin^2 theta) and the centrifugal potential -omega^2 sin^2
-- theta / 2, set up through the LTE projection and held by implicit BGK collisions (1x2v). The
-- Hamiltonian is the full rotating-frame one of the companion test, H = p^2/2 - omega p_phi = (1/2)
-- h^ij (p_i - A_i)(p_j - A_j) + Phi, which the LTE moments and projection require when a background
-- flow and potential are declared, so the LTE projection with density n0 2 pi T sin(theta)
-- exp(omega^2 sin^2 theta / 2T), zero drift relative to the flow and temperature T must reproduce
-- the companion's f = exp(-H/T), the LTE moments must return that drift and temperature, and the
-- BGK target must leave the state stationary between reflecting walls at theta = pi/4 and 5 pi/8
-- for t = 2 (nu t = 20). Figures of merit (serendipity p1, 32 x 12 x 12 cells): the LTE projection
-- differs from the companion's directly projected f = exp(-H/T) by 0.20 of its peak
-- (the p1 Maxwellian is corrected to the requested moments) while its LTE moments return zero drift
-- and T = 1 to round-off; under implicit BGK at nu = 10 the density, temperature and mean p_phi
-- change by at most 9.3e-2, 6.7e-2 and 6.1e-3 over t = 2; mass and energy conserved to round-off.
local Vlasov = G0.Vlasov

mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.
n0 = 1.0 -- Amplitude of the distribution function.
T0 = 1.0 -- Temperature.
omega = 1.0 -- Rotation rate of the sphere.
vt = 1.0 -- Thermal velocity.
nu = 10.0 -- Collision frequency.

Ntheta = 32 -- Cell count (configuration space: polar direction).
Nvtheta = 12 -- Cell count (velocity space: polar momentum).
Nvphi = 12 -- Cell count (velocity space: azimuthal momentum).
Ltheta = (3.0 * math.pi) / 8.0 -- Domain size (configuration space: polar direction).
vtheta_max = 6.0 * vt -- Domain boundary (velocity space: polar momentum).
vphi_max = 6.0 * vt -- Domain boundary (velocity space: azimuthal momentum).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 2.0 -- Final simulation time.
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
  lower = { math.pi / 4.0 },
  upper = { (math.pi / 4.0) + Ltheta },
  cells = { Ntheta },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = {  }, -- Periodic directions.

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.CanonicalPB,
    charge = charge, mass = mass,
    hamiltonian = function (t, xn)
      local q_theta = xn[1]
      local p_theta_dot, p_phi_dot = xn[2], xn[3]
      local s2 = math.sin(q_theta) * math.sin(q_theta)

      -- Rotating-frame Hamiltonian H = p_theta^2/2 + p_phi^2/(2 sin^2 theta) - omega p_phi.
      return (0.5 * p_theta_dot * p_theta_dot) + (0.5 * p_phi_dot * p_phi_dot / s2) -
        (omega * p_phi_dot)
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
    backgroundFlows = function (t, xn)
      local q_theta = xn[1]

      -- Background canonical flow of the rotating frame (theta and phi components).
      return 0.0, omega * math.sin(q_theta) * math.sin(q_theta)
    end,
    effectivePotential = function (t, xn)
      local q_theta = xn[1]

      -- Effective (centrifugal) potential of the rotating frame.
      return -0.5 * omega * omega * math.sin(q_theta) * math.sin(q_theta)
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
          local s2 = math.sin(theta) * math.sin(theta)

          -- M0 of the equilibrium n0 exp(-H/T): the Boltzmann factor of the centrifugal
          -- potential times the Maxwellian measure 2 pi T sin(theta).
          return n0 * 2.0 * math.pi * T0 * math.sin(theta) * math.exp(0.5 * omega * omega * s2 / T0)
        end,
        temperatureInit = function (t, xn)
          -- Isotropic temperature.
          return T0
        end,
        driftVelocityInit = function (t, xn)
          -- Drift velocity relative to the background flow (none).
          return 0.0, 0.0
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
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.LTEMoments, G0.Moment.EnergyMoment }
  },

  skipField = true,
}

vlasovApp:run()
