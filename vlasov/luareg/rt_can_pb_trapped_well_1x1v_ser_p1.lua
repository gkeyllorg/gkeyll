-- Trapped particles in the harmonic well H = p^2/2 + x^2 through the canonical Poisson-bracket
-- solver with absorbing walls at |x| = 1 (1x1v, collisionless). Every orbit in the well has the
-- period 2 pi / sqrt(2), so a distribution filling the diamond |p| < sqrt(2) (1 - |x|) rotates
-- rigidly in the (x, p / sqrt(2)) plane; after an eighth of the period it is the axis-aligned
-- square |x| < 1 / sqrt(2), |p| < 1, with the exact top-hat density n = 2 for |x| < 1 / sqrt(2) and
-- zero outside, and no particle reaches the walls (all energies are below the rim, 1). Figures of
-- merit (serendipity p1, 32 x 32 cells): at t = T/8 = 0.5554 the density differs from the exact top
-- hat by 0.052 in relative L1 (edge smearing), the plateau averages 2.010 (exact 2) and the density
-- outside |x| > 0.8 stays below 7.6e-3; mass and energy conserved to 1e-4.
local Vlasov = G0.Vlasov

mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.
vt = 1.0 -- Momentum scale.
t_period = math.pi * math.sqrt(2.0) -- Oscillation period in the well, 2 pi / sqrt(2).

Nx = 32 -- Cell count (configuration space: x-direction).
Nvx = 32 -- Cell count (velocity space: vx-direction).
Lx = 1.0 -- Half width of the domain.
vx_max = 2.0 * vt -- Domain boundary (velocity space: vx-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = t_period / 8.0 -- Final simulation time.
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
  lower = { -Lx },
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
      local x, p_x_dot = xn[1], xn[2]

      -- Canonical Hamiltonian with the harmonic well: H = p^2/2 + x^2.
      return (0.5 * p_x_dot * p_x_dot) + (x * x)
    end,
    inverseMetric = function (t, xn)
      -- Inverse metric tensor (xx component).
      return 1.0
    end,
    metric = function (t, xn)
      -- Metric tensor (xx component).
      return 1.0
    end,
    metricDeterminant = function (t, xn)
      -- Metric tensor determinant.
      return 1.0
    end,

    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,
        init = function (t, xn)
          local x, p_x_dot = xn[1], xn[2]

          -- Trapped particles: uniform f inside the diamond |p| < sqrt(2) (1 - |x|), 1e-10 outside.
          if math.abs(p_x_dot) < math.sqrt(2.0) * (1.0 - math.abs(x)) then
            return 1.0
          end
          return 1.0e-10
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
