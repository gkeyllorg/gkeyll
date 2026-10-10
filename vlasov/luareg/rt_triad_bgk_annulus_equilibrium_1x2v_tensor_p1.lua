-- Rotating equilibrium of a neutral gas in an annulus (polar coordinates, orthonormal triad) with
-- implicit BGK collisions through the triad Vlasov solver (1x2v). With uniform density, temperature
-- T = r^(3/2) and azimuthal velocity V_theta = sqrt(3 r^(3/2) / 2), the centrifugal force balances
-- the pressure gradient exactly, n V_theta^2 / r = d(nT)/dr, so the Euler state is stationary
-- between reflecting walls at r = 0.5 and 1.5; the kinetic solution must keep the density,
-- temperature and rotation profiles fixed over t = 0.5 (a sound crossing). Figures of merit
-- (tensor p1, 32 x 12 x 12 cells): over t = 0.5 the density changes by at most 1.7e-2 (relative),
-- the temperature by 1.9e-2 and the rotation velocity by 3.4e-3; mass and energy conserved to
-- round-off.
local Vlasov = G0.Vlasov

mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.
n0 = 1.0 -- Number density.
vt = 1.0 -- Thermal velocity.
nu = 15000.0 -- Collision frequency.

Nr = 32 -- Cell count (configuration space: radial direction).
Nvr = 12 -- Cell count (velocity space: radial direction).
Nvtheta = 12 -- Cell count (velocity space: angular direction).
Lr = 1.0 -- Domain size (configuration space: radial direction).
vr_max = 6.0 * vt -- Domain boundary (velocity space: radial direction).
vtheta_lo = -4.0 * vt -- Lower boundary of the azimuthal velocity domain.
vtheta_hi = 8.0 * vt -- Upper boundary of the azimuthal velocity domain.
poly_order = 1 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
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
  lower = { 0.5 },
  upper = { 0.5 + Lr },
  cells = { Nr },
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
    modelID = G0.Model.Triad,
    charge = charge, mass = mass,

    -- Geometry from the covariant tangent basis e_i and the orthonormal triad sigma_a.
    covTangentBasis = function (t, xn)
      local q_r = xn[1]
      local q_theta = 0.0

      -- Covariant tangent basis e_i . sigma_j of polar coordinates (row i, column j).
      return math.cos(q_theta), math.sin(q_theta), -q_r * math.sin(q_theta), q_r * math.cos(q_theta)
    end,
    triadBasis = function (t, xn)
      local q_r = xn[1]
      local q_theta = 0.0

      -- Orthonormal triad sigma_a . sigma_j (unit radial and azimuthal vectors).
      return math.cos(q_theta), math.sin(q_theta), -math.sin(q_theta), math.cos(q_theta)
    end,
    triadBasisGradient = function (t, xn)
      local q_r = xn[1]
      local q_theta = 0.0

      -- d(sigma_a . sigma_j)/dx^k for x^k = r, theta (four components each).
      return 0.0, 0.0, 0.0, 0.0, -math.sin(q_theta), math.cos(q_theta), -math.cos(q_theta),
        -math.sin(q_theta)
    end,

    -- Velocity space grid.
    lower = { -vr_max, vtheta_lo },
    upper = { vr_max, vtheta_hi },
    cells = { Nvr, Nvtheta },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,
        densityInit = function (t, xn)
          local r = xn[1]

          -- Uniform density times the metric determinant.
          return r * n0
        end,
        temperatureInit = function (t, xn)
          local r = xn[1]

          -- Temperature profile T = r^(3/2).
          return r * math.sqrt(r)
        end,
        driftVelocityInit = function (t, xn)
          local r = xn[1]

          -- Azimuthal rotation balancing the pressure gradient: n V_theta^2 / r = d(n T)/dr.
          return 0.0, math.sqrt(1.5 * r * math.sqrt(r) / n0)
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
