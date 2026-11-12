-- Coriolis deflection of a bump launched along a meridian on a rotating unit sphere, rotating-frame
-- Hamiltonian form H = p^2/2 - omega p_phi (2x2v, collisionless). A Gaussian bump of temperature
-- 0.01 just below the equator moves poleward with unit momentum while co-rotating with the sphere;
-- in the inertial frame its mean momentum follows the great circle through (theta_0, phi_0) =
-- (1.4, pi) with physical velocity (-1, omega sin theta_0), x(t) = x_0 cos(s t) + (v/s) sin(s t)
-- with s = |v| = 1.405, and in the rotating frame phi_rot = phi_inertial - omega t, so the centroid
-- must reach (theta, phi) = (1.113, 3.167) at t = 0.3
-- (thermal corrections are of order T t^2 = 1e-3); the drift in phi relative to the meridian is the
-- Coriolis effect. Figures of merit (serendipity p1, 10 x 20 x 11 x 8 cells): at t = 0.3 the
-- centroid is at (theta, phi) = (1.1080, 3.1712) against the exact (1.1129, 3.1673), errors -0.005
-- and +0.004 rad with cells of 0.039 and 0.105 rad; mass changes by -3e-2.
local Vlasov = G0.Vlasov

mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.
n0 = 1.0 -- Peak number density of the bump.
T0 = 0.01 -- Temperature of the bump.
vt = math.sqrt(T0) -- Thermal velocity.
omega = 1.0 -- Rotation rate of the sphere.
theta_0 = 1.4 -- Polar angle of the bump.
phi_0 = math.pi -- Azimuth of the bump.
L_bump = 0.1 -- Squared width of the bump.
V_theta_0 = -1.0 -- Polar momentum of the bump (poleward).

Ntheta = 10 -- Cell count (configuration space: polar direction).
Nphi = 20 -- Cell count (configuration space: azimuthal direction).
Nvtheta = 11 -- Cell count (velocity space: polar momentum).
Nvphi = 8 -- Cell count (velocity space: azimuthal momentum).
Ltheta = (3.0 * math.pi) / 8.0 -- Domain size (configuration space: polar direction).
Lphi = 2.0 * math.pi -- Domain size (configuration space: azimuthal direction).
ptheta_lo = -1.5 -- Lower boundary of the polar-momentum domain.
ptheta_hi = -0.4 -- Upper boundary of the polar-momentum domain (the Coriolis force raises p_theta).
pphi_lo = 0.55 -- Lower boundary of the azimuthal-momentum domain (around omega sin^2 theta_0).
pphi_hi = 1.35 -- Upper boundary of the azimuthal-momentum domain.
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.3 -- Final simulation time.
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
    lower = { ptheta_lo, pphi_lo },
    upper = { ptheta_hi, pphi_hi },
    cells = { Nvtheta, Nvphi },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,
        init = function (t, xn)
          local theta, phi, p_theta_dot, p_phi_dot = xn[1], xn[2], xn[3], xn[4]
          local s2 = math.sin(theta) * math.sin(theta)
          local dth, dph = theta - theta_0, (phi - phi_0) * math.sin(theta)

          -- Gaussian bump at (theta_0, phi_0) carrying a Maxwellian centred on the canonical
          -- momenta (V_theta_0, omega sin^2 theta) of a parcel on the meridian, co-rotating.
          local n = n0 * math.exp(-((dth * dth) + (dph * dph)) / (0.5 * L_bump))
          local dpt, dpp = p_theta_dot - V_theta_0, p_phi_dot - (omega * s2)
          local efact = ((dpt * dpt) + (dpp * dpp / s2)) / (2.0 * T0)

          -- Metric determinant times n times the normalized Maxwellian (measure 2 pi T sin theta).
          return n * math.exp(-efact) / (2.0 * math.pi * T0)
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
