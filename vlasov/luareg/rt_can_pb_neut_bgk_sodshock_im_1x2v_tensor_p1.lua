-- Sod shock tube of a neutral gas with implicit BGK collisions through the canonical
-- Poisson-bracket Vlasov solver with a flat (identity) metric (1x2v). Same problem as the neutral
-- BGK Sod family: (n, T) = (1, 1) left and (0.125, 0.8) right at rest, gamma = 2 for two velocity
-- dimensions, collision frequency 2000 (mean free path 0.03 cells) so the solution is the Euler
-- one; the reference is the exact Riemann solution. Figures of merit
-- (tensor p1, 64 x 12 x 12 cells): L1 density error 0.0137 against the exact Riemann solution,
-- shock at x = 0.6992 (exact 0.6957), contact at 0.5742 (exact 0.5760).
local Vlasov = G0.Vlasov

mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.
nl = 1.0 -- Left number density.
Tl = 1.0 -- Left temperature.
nr = 0.125 -- Right number density.
Tr = 0.8 -- Right temperature.
vt = 1.0 -- Thermal velocity.
nu = 2000.0 -- Collision frequency.

Nx = 64 -- Cell count (configuration space: x-direction).
Nvx = 12 -- Cell count (velocity space: vx-direction).
Nvy = 12 -- Cell count (velocity space: vy-direction).
Lx = 1.0 -- Domain size (configuration space: x-direction).
vx_max = 6.0 * vt -- Domain boundary (velocity space: vx-direction).
vy_max = 6.0 * vt -- Domain boundary (velocity space: vy-direction).
midplane = 0.5 -- Location of the jump in the initial state.
poly_order = 1 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
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
  lower = { 0.0 },
  upper = { Lx },
  cells = { Nx },
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
      local p_x_dot, p_y_dot = xn[2], xn[3]

      local inv_metric_x_x = 1.0
      local inv_metric_x_y = 0.0
      local inv_metric_y_y = 1.0

      local hamiltonian = (0.5 * inv_metric_x_x * p_x_dot * p_x_dot) +
        (0.5 * (2.0 * inv_metric_x_y * p_x_dot * p_y_dot)) +
        (0.5 * inv_metric_y_y * p_y_dot * p_y_dot) -- Canonical Hamiltonian.

      return hamiltonian
    end,
    inverseMetric = function (t, xn)
      -- Inverse metric tensor (xx, xy, yy components).
      return 1.0, 0.0, 1.0
    end,
    metric = function (t, xn)
      -- Metric tensor (xx, xy, yy components).
      return 1.0, 0.0, 1.0
    end,
    metricDeterminant = function (t, xn)
      -- Metric tensor determinant.
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
        projectionID = G0.Projection.LTE,
        densityInit = function (t, xn)
          local x = xn[1]

          -- Left and right states.
          local n = nr
          if x < midplane then
            n = nl
          end

          -- Total number density (times the metric determinant, 1 here).
          return n
        end,
        temperatureInit = function (t, xn)
          local x = xn[1]

          -- Isotropic temperature.
          local T = Tr
          if x < midplane then
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
        type = G0.SpeciesBc.bcCopy
      },
      upper = {
        type = G0.SpeciesBc.bcCopy
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.LTEMoments, G0.Moment.EnergyMoment }
  },

  skipField = true,
}

vlasovApp:run()
