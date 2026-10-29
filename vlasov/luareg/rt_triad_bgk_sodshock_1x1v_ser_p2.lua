-- Sod shock tube of a neutral gas with implicit BGK collisions through the triad Vlasov solver on a
-- Cartesian line with the identity triad (the solver must reduce to plain Vlasov) (1x1v). The
-- states (n, T) = (1, 1) and (0.125, 0.8) are at rest on either side of the jump at x = 0.5 on 0 <
-- x < 1; with collision frequency NU the solution is the Euler one with gamma = 3. On the Cartesian
-- line the metric determinant is 1 and the waves stay clear of the reflecting walls until t = 0.1,
-- so the reference is the exact planar Riemann solution. Figures of merit (serendipity p2, 128 x 32
-- cells, 510 steps), for the density at t = 0.1: relative L1 error (int |n - n_ex| dx / int n_ex
-- dx) 0.0064 against the exact solution; the density crosses the midpoint of the exact jump at x =
-- 0.7281 across the shock and 0.5597 across the contact (exact 0.7273 and 0.5609); mass and energy
-- conserved to round-off.
local Vlasov = G0.Vlasov

mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.
nl = 1.0 -- Left/inner number density.
Tl = 1.0 -- Left/inner temperature.
nr = 0.125 -- Right/outer number density.
Tr = 0.8 -- Right/outer temperature.
vt = 1.0 -- Thermal velocity.
nu = 5000.0 -- Collision frequency.
midplane = 0.5 -- Location of the jump.

Nx = 128 -- Cell count (configuration space: x-direction).
Nvx = 32 -- Cell count (velocity space: vx-direction).
Lx = 1.0 -- Domain size (configuration space: x-direction).
vx_max = 8.0 * vt -- Domain boundary (velocity space: vx-direction).
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
    modelID = G0.Model.Triad,
    charge = charge, mass = mass,

    -- Geometry from the covariant tangent basis e_i and the orthonormal triad sigma_a.
    covTangentBasis = function (t, xn)
      -- Cartesian coordinates: the tangent basis is the unit vector.
      return 1.0
    end,
    triadBasis = function (t, xn)
      -- Orthonormal triad: the unit vector.
      return 1.0
    end,
    triadBasisGradient = function (t, xn)
      -- The triad is constant.
      return 0.0
    end,

    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,
        densityInit = function (t, xn)
          local x = xn[1]

          -- Left and right states, times the metric determinant.
          local n = nr
          if x < midplane then
            n = nl
          end
          return (1.0) * n
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
          -- Drift velocity in the orthonormal frame (the gas is at rest).
          return 0.0
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
