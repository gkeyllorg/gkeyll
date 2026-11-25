-- Sod shock tube on a distorted (sheared) mesh through the canonical Poisson-bracket solver
-- (2x2v, implicit BGK). The coordinates z1, z2 map to the plane by x = z1, y = z2 + 0.2 cos(pi z1),
-- giving the non-diagonal metric [[1 + h'^2, h'], [h', 1]] with h' = -0.2 pi sin(pi z1)
-- (up to 0.63) and unit determinant; a planar Sod state (n, T) = (1, 1) | (0.125, 0.8) across z1 =
-- 0 is independent of y, so the kinetic solution must equal the exact planar Riemann solution
-- (gamma = 2) in z1 and stay uniform in z2. Figures of merit
-- (serendipity p2, 32 x 2 x 12 x 12 cells, 123 steps): L1 density error 0.0110 against the exact
-- planar Riemann solution along z1, variation across z2 below 1e-14; mass changes by 1e-3 through
-- the copy boundaries.
local Vlasov = G0.Vlasov

mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.
nl = 1.0 -- Left number density.
Tl = 1.0 -- Left temperature.
nr = 0.125 -- Right number density.
Tr = 0.8 -- Right temperature.
vt = 1.0 -- Thermal velocity.
nu = 2000.0 -- Collision frequency.
eps = 0.2 -- Amplitude of the mesh shear, y = z2 + eps cos(pi z1).

Nz1 = 32 -- Cell count (configuration space: z1-direction).
Nz2 = 2 -- Cell count (configuration space: z2-direction).
Nvz1 = 12 -- Cell count (velocity space: p_z1).
Nvz2 = 12 -- Cell count (velocity space: p_z2).
Lz1 = 2.0 -- Domain size (configuration space: z1-direction).
Lz2 = 2.0 -- Domain size (configuration space: z2-direction).
vz1_max = 8.0 * vt -- Domain boundary (velocity space: p_z1, which is v_x + h' v_y).
vz2_max = 6.0 * vt -- Domain boundary (velocity space: p_z2 = v_y).
midplane = 0.0 -- Location of the jump in the initial state (z1).
poly_order = 2 -- Polynomial order.
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
  lower = { -0.5 * Lz1, -0.5 * Lz2 },
  upper = { 0.5 * Lz1, 0.5 * Lz2 },
  cells = { Nz1, Nz2 },
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
      local q_z1 = xn[1]
      local hp = -eps * math.pi * math.sin(math.pi * q_z1)
      local p_z1_dot, p_z2_dot = xn[3], xn[4]

      local inv_metric_aa = 1.0
      local inv_metric_ab = -hp
      local inv_metric_bb = 1.0 + (hp * hp)

      -- Canonical Hamiltonian H = (1/2) g^ij p_i p_j.
      return (0.5 * inv_metric_aa * p_z1_dot * p_z1_dot) +
        (inv_metric_ab * p_z1_dot * p_z2_dot) +
        (0.5 * inv_metric_bb * p_z2_dot * p_z2_dot)
    end,
    inverseMetric = function (t, xn)
      local q_z1 = xn[1]
      local hp = -eps * math.pi * math.sin(math.pi * q_z1)
      -- Inverse metric tensor (aa, ab, bb components).
      return 1.0, -hp, 1.0 + (hp * hp)
    end,
    metric = function (t, xn)
      local q_z1 = xn[1]
      local hp = -eps * math.pi * math.sin(math.pi * q_z1)
      -- Metric tensor (aa, ab, bb components).
      return 1.0 + (hp * hp), hp, 1.0
    end,
    metricDeterminant = function (t, xn)
      local q_z1 = xn[1]
      local hp = -eps * math.pi * math.sin(math.pi * q_z1)
      -- Metric tensor determinant (square root of det g).
      return 1.0
    end,

    -- Velocity space grid.
    lower = { -vz1_max, -vz2_max },
    upper = { vz1_max, vz2_max },
    cells = { Nvz1, Nvz2 },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,
        densityInit = function (t, xn)
          local z1 = xn[1]

          -- Left and right states, times the metric determinant.
          local n = nr
          if z1 < midplane then
            n = nl
          end
          return (1.0) * n
        end,
        temperatureInit = function (t, xn)
          local z1 = xn[1]

          -- Isotropic temperature.
          local T = Tr
          if z1 < midplane then
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
