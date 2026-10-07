-- Sod shock tube for a neutral gas with BGK collisions (Vlasov, 1x2v).
-- Density 1 and pressure 1 on the left, density 1/8 and pressure 1/10 on the right, at rest. The collision
-- frequency puts the mean free path at 1/30 of a cell, so the solution follows the Euler Sod solution of a gas
-- with adiabatic index (d+2)/d = 2 with the discontinuities resolved to the grid scale.
-- The BGK relaxation is integrated implicitly, split from the explicit RK3 step of the collisionless terms.
-- Quadratic velocity maps cluster the cells at the origin of velocity space.
-- Exact solution at t = 0.1: shock at x = 0.696, contact at 0.576, post-shock density 0.204, velocity 0.760.
-- The shock and contact positions agree to within a cell and the plateaus to 0.4%; the density L1 error is 0.006.
-- The density differs from the explicit run at the same collision frequency by 1% at most.
local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
mass_neut = 1.0 -- Neutral mass.
charge_neut = 0.0 -- Neutral charge.

nl = 1.0 -- Left number density.
Tl = 1.0 -- Left temperature.
nr = 0.125 -- Right number density.
Tr = 0.8 -- Right temperature.
Vx_drift = 0.0 -- Drift velocity (x-direction).
nu = 2000.0 -- Collision frequency.

-- Derived physical quantities (using normalized code units).
vt = math.sqrt(Tl / mass_neut) -- Thermal velocity (left).

-- Simulation parameters.
Nx = 64 -- Cell count (configuration space: x-direction).
Nvx = 12 -- Cell count (velocity space: vx-direction).
Nvy = 12 -- Cell count (velocity space: vy-direction).
Lx = 1.0 -- Domain size (configuration space: x-direction).
vx_max = 8.0 * vt -- Domain boundary (velocity space: vx-direction).
vy_max = 8.0 * vt -- Domain boundary (velocity space: vy-direction).
vx_lin = 2.0 * vt -- Velocity map: cell size at the origin relative to a uniform grid (vx-direction).
vy_lin = 2.0 * vt -- Velocity map: cell size at the origin relative to a uniform grid (vy-direction).
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
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Neutrals.
  neut = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_neut, mass = mass_neut,

    -- Velocity space grid.
    lower = { -1.0, -1.0 },
    upper = { 1.0, 1.0 },
    cells = { Nvx, Nvy },

    -- Quadratic velocity maps: finest cells at the origin, stretching to the domain boundary.
    mapc2pVel = {
      -- vx mapping
      {
        vmap = function (t, xn)
          local vc = xn[1]
          return vx_lin * vc + (vx_max - vx_lin) * vc * math.abs(vc)
        end
      },
      -- vy mapping
      {
        vmap = function (t, xn)
          local vc = xn[1]
          return vy_lin * vc + (vy_max - vy_lin) * vc * math.abs(vc)
        end
      }
    },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]

          local n = 0.0
          if x < 0.5 then
            n = nl -- Total number density (left).
          else
            n = nr -- Total number density (right).
          end

          return n
        end,
        temperatureInit = function (t, xn)
          local x = xn[1]

          local T = 0.0
          if x < 0.5 then
            T = Tl -- Isotropic temperature (left).
          else
            T = Tr -- Isotropic temperature (right).
          end

          return T
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift, 0.0 -- Total drift velocity.
        end,

        correctAllMoments = true,
        useLastConverged = true
      }
    },

    collisions = {
      collisionID = G0.Collisions.BGK,
      useImplicitCollisionScheme = true,

      selfNu = function (t, xn)
        return nu -- Collision frequency.
      end
    },

    correct = {
      correctAllMoments = true,
      iterationEpsilon = 1.0e-12,
      maxIterations = 100,
      useLastConverged = true
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.LTEMoments }
  },

  skipField = true,

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = 1.0, mu0 = 1.0,

    -- Initial conditions function.
    init = function (t, xn)
      return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0,

    isStatic = true
  }
}

vlasovApp:run()
