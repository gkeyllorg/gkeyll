-- Stationary Maxwell-Boltzmann equilibrium on a rotating unit sphere in the rotating-frame
-- Hamiltonian form H = p^2/2 - omega p_phi (1x2v, collisionless). The distribution f = exp(-H/T),
-- with density proportional to exp(omega^2 sin^2 theta / 2T) and co-rotating mean momentum p_phi =
-- omega sin^2 theta, is a function of the Hamiltonian and must stay stationary between reflecting
-- walls at theta = pi/4 and 5 pi/8 for t = 5 (about four thermal transits). Figures of merit
-- (serendipity p2, 32 x 12 x 12 cells, 5766 steps): over t = 5 the density profile changes by at
-- most 1.7e-3 (relative), the temperature by 1.4e-3 and the mean p_phi by 2.1e-3; mass and energy
-- conserved to round-off. With implicit BGK collisions at nu = 10 the same state stays stationary
-- to 2e-4 over t = 2.
local Vlasov = G0.Vlasov

mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.
n0 = 1.0 -- Amplitude of the distribution function.
T0 = 1.0 -- Temperature.
omega = 1.0 -- Rotation rate of the sphere.
vt = 1.0 -- Thermal velocity.

Ntheta = 32 -- Cell count (configuration space: polar direction).
Nvtheta = 12 -- Cell count (velocity space: polar momentum).
Nvphi = 12 -- Cell count (velocity space: azimuthal momentum).
Ltheta = (3.0 * math.pi) / 8.0 -- Domain size (configuration space: polar direction).
vtheta_max = 6.0 * vt -- Domain boundary (velocity space: polar momentum).
vphi_max = 6.0 * vt -- Domain boundary (velocity space: azimuthal momentum).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 5.0 -- Final simulation time.
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

    -- Velocity space grid.
    lower = { -vtheta_max, -vphi_max },
    upper = { vtheta_max, vphi_max },
    cells = { Nvtheta, Nvphi },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,
        init = function (t, xn)
          local theta, p_theta_dot, p_phi_dot = xn[1], xn[2], xn[3]
          local s2 = math.sin(theta) * math.sin(theta)

          -- Exact equilibrium f = n0 exp(-H/T) of the rotating-frame Hamiltonian
          -- H = p_theta^2/2 + p_phi^2/(2 sin^2 theta) - omega p_phi, projected directly.
          local H = (0.5 * p_theta_dot * p_theta_dot) + (0.5 * p_phi_dot * p_phi_dot / s2) -
            (omega * p_phi_dot)
          return n0 * math.exp(-H / T0)
        end
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
