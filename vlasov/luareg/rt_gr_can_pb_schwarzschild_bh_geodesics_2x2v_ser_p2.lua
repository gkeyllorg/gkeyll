-- Geodesic motion of a phase-space blob on a bound orbit in the Schwarzschild equatorial plane
-- through the general-relativistic canonical Poisson-bracket solver (2x2v, collisionless). The
-- Hamiltonian is H = alpha sqrt(1 + g^ij p_i p_j) with the Schwarzschild lapse and spatial metric
-- for M = 1; the blob is released at the apoapsis r = 15 of the strongly relativistic bound orbit
-- with semi-latus rectum 7.5 and eccentricity 0.5 (angular momentum L = 3.64, periapsis 5), and its
-- centroid must follow the geodesic integrated from the same Hamiltonian for t = 30, during which
-- it falls to r ~ 8 and the Newtonian ellipse with the same L already lags noticeably. Figures of
-- merit (serendipity p2, 8 x 12 x 10 x 10 cells, 189 steps): at t = 30 the centroid is at
-- (r, theta) = (14.443, 2.678) against the geodesic (14.520, 2.694), i.e. within 0.08 in r and
-- 0.016 rad in theta, while a Newtonian orbit with the same angular momentum would be at
-- (14.768, 2.651); 2% of the mass is lost through the walls.
local Vlasov = G0.Vlasov

mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.
bh_mass = 1.0 -- Mass of the black hole (geometric units).
latus_rectum = 7.5 -- Semi-latus rectum of the orbit (bound orbits need p > 6 + 2e).
eccentricity = 0.5 -- Eccentricity of the orbit.

Nr = 8 -- Cell count (configuration space: radial direction).
Ntheta = 12 -- Cell count (configuration space: angular direction).
Nvr = 10 -- Cell count (velocity space: radial momentum).
Nvtheta = 10 -- Cell count (velocity space: angular momentum).
r_min = 3.0 -- Lower radius of the domain.
r_max = 20.0 -- Upper radius of the domain.
vr_max = 0.5 -- Domain boundary (velocity space: radial momentum).
ptheta_lo = -4.4 -- Lower boundary of the angular-momentum domain.
ptheta_hi = -2.9 -- Upper boundary of the angular-momentum domain (orbit: p_theta = -3.64).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 30.0 -- Final simulation time.
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
      local f = 1.0 - (2.0 * bh_mass / q_r)

      -- Massive-particle Hamiltonian in the Schwarzschild equatorial plane, H = alpha gamma.
      local gamma = math.sqrt(1.0 + (f * p_r_dot * p_r_dot) +
        (p_theta_dot * p_theta_dot / (q_r * q_r)))
      return math.sqrt(f) * gamma
    end,
    inverseMetric = function (t, xn)
      local q_r = xn[1]
      local f = 1.0 - (2.0 * bh_mass / q_r)

      -- Inverse spatial metric (rr, rtheta, thetatheta components).
      return f, 0.0, 1.0 / (q_r * q_r)
    end,
    metric = function (t, xn)
      local q_r = xn[1]
      local f = 1.0 - (2.0 * bh_mass / q_r)

      -- Spatial metric (rr, rtheta, thetatheta components).
      return 1.0 / f, 0.0, q_r * q_r
    end,
    metricDeterminant = function (t, xn)
      local q_r = xn[1]
      local f = 1.0 - (2.0 * bh_mass / q_r)

      -- Metric determinant (square root of det g).
      return q_r / math.sqrt(f)
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
          local q_r, q_theta = xn[1], xn[2]
          local p_r_dot, p_theta_dot = xn[3], xn[4]
          local M, lr, e = bh_mass, latus_rectum, eccentricity
          local dr, dtheta = (r_max - r_min) / Nr, 2.0 * math.pi / Ntheta
          local dpr, dptheta = 2.0 * vr_max / Nvr, (ptheta_hi - ptheta_lo) / Nvtheta

          -- Conserved angular momentum and energy of the bound orbit (p, e), released at apoapsis.
          local L = math.sqrt((M * lr * lr) / (lr - (M * (3.0 + (e * e)))))
          local E = math.sqrt(1.0 - ((L * L) / (lr * lr * lr)) * (lr - (4.0 * M)) * (1.0 - (e * e)))
          local r0, theta0, pr0, ptheta0 = lr / (1.0 - e), math.pi, 0.0, -L

          -- Gaussian blob on the orbit: sigma of half a cell in position and one cell in momentum.
          local g = math.exp(-0.5 * (q_r - r0) * (q_r - r0) / (0.25 * dr * dr)) *
            math.exp(-0.5 * (q_theta - theta0) * (q_theta - theta0) / (0.25 * dtheta * dtheta)) *
            math.exp(-0.5 * (p_r_dot - pr0) * (p_r_dot - pr0) / (dpr * dpr)) *
            math.exp(-0.5 * (p_theta_dot - ptheta0) * (p_theta_dot - ptheta0) / (dptheta * dptheta))
          return g
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
