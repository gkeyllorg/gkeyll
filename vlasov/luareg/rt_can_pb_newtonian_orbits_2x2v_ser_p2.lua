-- Ring of circular Kepler orbits in the Newtonian point-mass potential, H = p_r^2/2 +
-- p_theta^2/(2 r^2) - 1/r, through the canonical Poisson-bracket solver (2x2v, collisionless). The
-- ring f = exp(-(H - H_circ(L)) / T) exp(-(L - 1.05)^2 / 2 (0.04)^2) with T = 0.02
-- (circular-orbit radius L^2 = 1.1, radial width 0.2) is a smooth function of the integrals of
-- motion, hence an exact stationary state; its density profile must be unchanged after t = 2, a
-- quarter of the orbital period 2 pi r^(3/2) at the ring centre, with absorbing walls at r = 0.5
-- and 3 where the Boltzmann factor is below 2e-4. Figures of merit
-- (serendipity p2, 20 x 2 x 10 x 8 cells, 705 steps): at t = 2 the azimuthally averaged density
-- profile differs from its initial projection by 2.0e-3 in relative L1 (3.0e-3 at the peak) and
-- stays azimuthally uniform to 4e-5; mass and energy change by 8e-5 and 4e-5
-- (absorption of the Boltzmann tail at the outer wall).
local Vlasov = G0.Vlasov

mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.
L_0 = 1.05 -- Central angular momentum of the ring (circular-orbit radius L^2 = 1.1).
sigma_L = 0.04 -- Width of the Gaussian angular-momentum window.
T0 = 0.02 -- Temperature of the Maxwell-Boltzmann factor around each circular orbit.

Nr = 20 -- Cell count (configuration space: radial direction).
Ntheta = 2 -- Cell count (configuration space: angular direction).
Nvr = 10 -- Cell count (velocity space: radial momentum).
Nvtheta = 8 -- Cell count (velocity space: angular momentum).
r_min = 0.5 -- Lower radius of the domain.
r_max = 3.0 -- Upper radius of the domain.
vr_max = 0.6 -- Domain boundary (velocity space: radial momentum).
ptheta_lo = 0.85 -- Lower boundary of the angular-momentum domain.
ptheta_hi = 1.25 -- Upper boundary of the angular-momentum domain.
poly_order = 2 -- Polynomial order.
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
  lower = { r_min, 0.0 },
  upper = { r_max, 2.0 * math.pi },
  cells = { Nr, Ntheta },
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
      local q_r = xn[1]
      local p_r_dot, p_theta_dot = xn[3], xn[4]

      -- Newtonian point-mass Hamiltonian: H = p_r^2/2 + p_theta^2/(2 r^2) - 1/r.
      return (0.5 * p_r_dot * p_r_dot) + (0.5 * p_theta_dot * p_theta_dot / (q_r * q_r)) -
        (1.0 / q_r)
    end,
    inverseMetric = function (t, xn)
      local q_r = xn[1]
      -- Inverse metric tensor (aa, ab, bb components).
      return 1.0, 0.0, 1.0 / (q_r * q_r)
    end,
    metric = function (t, xn)
      local q_r = xn[1]
      -- Metric tensor (aa, ab, bb components).
      return 1.0, 0.0, q_r * q_r
    end,
    metricDeterminant = function (t, xn)
      local q_r = xn[1]
      -- Metric tensor determinant (square root of det g).
      return q_r
    end,

    -- Velocity space grid.
    lower = { -vr_max, ptheta_lo },
    upper = { vr_max, ptheta_hi },
    cells = { Nvr, Nvtheta },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,
        init = function (t, xn)
          local q_r = xn[1]
          local p_r_dot, p_theta_dot = xn[3], xn[4]
          local H = (0.5 * p_r_dot * p_r_dot) + (0.5 * p_theta_dot * p_theta_dot / (q_r * q_r)) -
            (1.0 / q_r)
          local H_circ = -0.5 / (p_theta_dot * p_theta_dot)

          -- Stationary ring f = exp(-(H - H_circ(L)) / T) exp(-(L - L_0)^2 / 2 sigma_L^2): a
          -- smooth function of the integrals H and L, peaked on the circular orbits r = L^2.
          -- No metric factor: f is a phase-space density in the canonical coordinates
          -- (r, theta, p_r, p_theta), so any function of H and L alone is exactly stationary.
          local dL = p_theta_dot - L_0
          return math.exp(-(H - H_circ) / T0) * math.exp(-(dL * dL) / (2.0 * sigma_L * sigma_L))
        end
      }
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcAbsorb
      },
      upper = {
        type = G0.SpeciesBc.bcAbsorb
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.LTEMoments, G0.Moment.EnergyMoment }
  },

  skipField = true,
}

vlasovApp:run()
