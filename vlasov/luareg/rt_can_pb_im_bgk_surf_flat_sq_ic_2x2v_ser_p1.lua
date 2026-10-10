-- Relaxation of a square velocity distribution to a Maxwellian by implicit BGK collisions through
-- the canonical Poisson-bracket solver with a flat metric (2x2v, periodic). The uniform square
-- |p_x|, |p_y| < 1 of height 0.5 has density 2, zero flow and temperature 1/3; with nu dt >> 1 the
-- implicit scheme must return the Maxwellian with exactly these moments in every cell, conserving
-- density, momentum and energy. Figures of merit (serendipity p1, 2 x 2 x 16 x 16 cells): the
-- projected square has n = 2.5313 and T = 0.4939 (quadrature of the discontinuity); after
-- relaxation every cell carries exactly these moments and the cell averages of f agree with those
-- of the Maxwellian(n, T) to 7.0e-2 (cell averages of the p1 projection) of its peak; mass and
-- energy conserved to round-off.
local Vlasov = G0.Vlasov

mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.
n0 = 0.5 -- Distribution function value inside the square.
vt = 1.0 -- Thermal velocity.
nu = 15000.0 -- Collision frequency.

Nx = 2 -- Cell count (configuration space: x-direction).
Ny = 2 -- Cell count (configuration space: y-direction).
Nvx = 16 -- Cell count (velocity space: vx-direction).
Nvy = 16 -- Cell count (velocity space: vy-direction).
Lx = 1.0 -- Domain size (configuration space: x-direction).
Ly = 1.0 -- Domain size (configuration space: y-direction).
vx_max = 6.0 * vt -- Domain boundary (velocity space: vx-direction).
vy_max = 6.0 * vt -- Domain boundary (velocity space: vy-direction).
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
  lower = { 0.0, 0.0 },
  upper = { Lx, Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = { 1, 2 }, -- Periodic directions.

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.CanonicalPB,
    charge = charge, mass = mass,
    hamiltonian = function (t, xn)
      local p_x_dot, p_y_dot = xn[3], xn[4]

      local inv_metric_aa = 1.0
      local inv_metric_ab = 0.0
      local inv_metric_bb = 1.0

      -- Canonical Hamiltonian H = (1/2) g^ij p_i p_j.
      return (0.5 * inv_metric_aa * p_x_dot * p_x_dot) +
        (inv_metric_ab * p_x_dot * p_y_dot) +
        (0.5 * inv_metric_bb * p_y_dot * p_y_dot)
    end,
    inverseMetric = function (t, xn)
      -- Inverse metric tensor (aa, ab, bb components).
      return 1.0, 0.0, 1.0
    end,
    metric = function (t, xn)
      -- Metric tensor (aa, ab, bb components).
      return 1.0, 0.0, 1.0
    end,
    metricDeterminant = function (t, xn)
      -- Metric tensor determinant (square root of det g).
      return 1.0
    end,

    -- Velocity space grid.
    lower = { -vx_max, -vy_max },
    upper = { vx_max, vy_max },
    cells = { Nvx, Nvy },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,
        init = function (t, xn)
          local p_x_dot, p_y_dot = xn[3], xn[4]

          -- Square distribution |p_x| < 1, |p_y| < 1 of height n0 (metric determinant 1).
          if (math.abs(p_x_dot) < 1.0) and (math.abs(p_y_dot) < 1.0) then
            return n0
          end
          return 0.0
        end
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

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.LTEMoments, G0.Moment.EnergyMoment }
  },

  skipField = true,
}

vlasovApp:run()
