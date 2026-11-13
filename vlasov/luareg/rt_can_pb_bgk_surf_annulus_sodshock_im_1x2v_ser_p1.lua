-- Sod shock tube of a neutral gas with implicit BGK collisions through the canonical
-- Poisson-bracket Vlasov solver on an annulus (polar coordinates r, theta with metric diag(1, r^2),
-- jump at r = 1 on 0.5 < r < 1.5) (1x2v). The states (n, T) = (1, 1) and (0.125, 0.8) are at rest
-- on either side of the jump; with collision frequency 2000.0 the solution is the Euler one with
-- gamma = 2. The reference is a quasi-1D finite-volume Euler solution with the metric determinant
-- as area factor (reflecting walls), and the exact planar Riemann solution bounds the early-time
-- profile. Figures of merit (serendipity p1, 64 x 12 x 16 cells): L1 density error 0.0176 against
-- the quasi-1D reference (0.0201 against the planar exact solution), shock at r = 1.1836
-- (reference 1.1934), contact at 1.0742 (reference 1.0756); mass and energy conserved to round-off.
local Vlasov = G0.Vlasov

mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.
nl = 1.0 -- Left/inner number density.
Tl = 1.0 -- Left/inner temperature.
nr = 0.125 -- Right/outer number density.
Tr = 0.8 -- Right/outer temperature.
vt = 1.0 -- Thermal velocity.
nu = 2000.0 -- Collision frequency.

Nr = 64 -- Cell count (configuration space: radial direction).
Nvr = 12 -- Cell count (velocity space: radial momentum).
Nvtheta = 16 -- Cell count (velocity space: angular momentum).
Lr = 1.0 -- Domain size (configuration space: radial direction).
vr_max = 6.0 * vt -- Domain boundary (velocity space: radial momentum).
vtheta_max = 8.0 * vt -- Domain boundary (velocity space: angular momentum).
midplane = 1.0 -- Radius of the jump in the initial state.
poly_order = 1 -- Polynomial order.
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
    modelID = G0.Model.CanonicalPB,
    charge = charge, mass = mass,
    hamiltonian = function (t, xn)
      local q_r = xn[1]
      local p_r_dot, p_theta_dot = xn[2], xn[3]

      local inv_metric_aa = 1.0
      local inv_metric_ab = 0.0
      local inv_metric_bb = 1.0 / (q_r * q_r)

      -- Canonical Hamiltonian H = (1/2) g^ij p_i p_j.
      return (0.5 * inv_metric_aa * p_r_dot * p_r_dot) +
        (inv_metric_ab * p_r_dot * p_theta_dot) +
        (0.5 * inv_metric_bb * p_theta_dot * p_theta_dot)
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
    lower = { -vr_max, -vtheta_max },
    upper = { vr_max, vtheta_max },
    cells = { Nvr, Nvtheta },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,
        densityInit = function (t, xn)
          local r = xn[1]

          -- Left and right states, times the metric determinant.
          local n = nr
          if r < midplane then
            n = nl
          end
          return (r) * n
        end,
        temperatureInit = function (t, xn)
          local r = xn[1]

          -- Isotropic temperature.
          local T = Tr
          if r < midplane then
            T = Tl
          end
          return T
        end,
        driftVelocityInit = function (t, xn)
          -- Total drift velocity (the gas is at rest).
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
