-- Relaxation under implicit BGK collisions and the implicit BGK-form temperature source nu_Q
-- (f_M(n, u, T_Q) - f) together, with the Vlasov solver (1x2v):
-- rt_vlasov_temp_source_relax_im_1x2v.c with its collision frequency raised from nu_Q to 1000 nu_Q
-- (stiff self-collisions).
-- Two homogeneous neutral species start as bi-Maxwellians drifting along x with T_xx = 1.4 and T_yy
-- = 0.6 (T = 1). The source of the first has T_Q = 10 (heating) and that of the second T_Q = 0.1
-- (cooling). Source rate 100 (z = nu_Q dt between 2.45 and 2.55) and collision frequency 1e5 (nu_c
-- dt about 2500).
-- The two operators take one backward-Euler step together in each SSP-RK3 stage, f' = (f + dt nu_c
-- f_M(T') + dt nu_Q f_M(T_Q))/(1 + dt (nu_c + nu_Q)), where T' = (T + dt nu_Q T_Q)/(1 + dt nu_Q) is
-- the temperature the source gives the species over the step. The collisions therefore leave the
-- temperature to the source, which multiplies T - T_Q by R(z) = [1/3 + (2/3) (3/4 + (1/4)/(1 +
-- z))/(1 + z/4)]/(1 + 2 z/3) every step, the anisotropy decays by R((nu_c + nu_Q) dt), and the
-- stage ends on the Maxwellian of the new temperature up to the ratio nu_Q/nu_c.
-- Figures of merit (serendipity p2, 2 x 48 x 48 cells, 4 steps): density and momentum are conserved
-- to round-off (1e-14). Each step multiplies T - T_Q by R(z) to 2e-13. The first step multiplies
-- the pressure anisotropy P_xx - P_yy by R((nu_c + nu_Q) dt) to 2e-6 (heated) and 4e-8 (cooled).
-- After the first step the cell averages of the distribution differ from those of the Maxwellian
-- with its density, drift and temperature by 1.1e-3 (heated) and 1.5e-3 (cooled) of its peak, where
-- the initial state differs by 0.16.

local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
mass = 1.0 -- Species mass.
charge = 0.0 -- Species charge.

n0 = 1.0 -- Reference number density.
Txx = 1.4 -- Initial temperature along x.
Tyy = 0.6 -- Initial temperature along y.
T0 = 0.5 * (Txx + Tyy) -- Initial temperature.
ux0 = 0.5 -- Initial drift velocity (x-direction).
nu_Q = 100.0 -- Rate of the temperature source.
T_Q_heat = 10.0 -- Source temperature of the heated species.
T_Q_cool = 0.1 -- Source temperature of the cooled species.
nu_c = 1.0e5 -- Collision frequency.

-- Simulation parameters.
Nx = 2 -- Cell count (configuration space: x-direction).
Nvx = 48 -- Cell count (velocity space: vx-direction).
Nvy = 48 -- Cell count (velocity space: vy-direction).
Lx = 5.0 -- Domain size (configuration space: x-direction).
v_max_heat = 20.0 -- Domain boundary of the heated species (velocity space).
v_max_cool = 8.0 -- Domain boundary of the cooled species (velocity space).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.1 -- Final simulation time.
num_frames = 4 -- Number of output frames.
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
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Heated species.
  heat = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge, mass = mass,

    -- Velocity space grid.
    lower = { -v_max_heat, -v_max_heat },
    upper = { v_max_heat, v_max_heat },
    cells = { Nvx, Nvy },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local vx, vy = xn[2], xn[3]

          -- Drifting bi-Maxwellian.
          local f = n0 * mass / (2.0 * pi * math.sqrt(Txx * Tyy)) *
            math.exp(-mass * (vx - ux0) * (vx - ux0) / (2.0 * Txx) - mass * vy * vy / (2.0 * Tyy))

          return f
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.BGK,

      selfNu = function (t, xn)
        return nu_c -- Collision frequency.
      end,

      useImplicitCollisionScheme = true
    },

    sourceBGK = {
      rateProfile = function (t, xn)
        return nu_Q -- Rate of the temperature source.
      end,
      tempShape = function (t, xn)
        return T_Q_heat -- Source temperature.
      end,
      useImplicitScheme = true
    },

    correct = {
      correctAllMoments = true,
      iterationEpsilon = 1.0e-12,
      maxIterations = 100,
      useLastConverged = true
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2ij }
  },

  -- Cooled species.
  cool = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge, mass = mass,

    -- Velocity space grid.
    lower = { -v_max_cool, -v_max_cool },
    upper = { v_max_cool, v_max_cool },
    cells = { Nvx, Nvy },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local vx, vy = xn[2], xn[3]

          -- Drifting bi-Maxwellian.
          local f = n0 * mass / (2.0 * pi * math.sqrt(Txx * Tyy)) *
            math.exp(-mass * (vx - ux0) * (vx - ux0) / (2.0 * Txx) - mass * vy * vy / (2.0 * Tyy))

          return f
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.BGK,

      selfNu = function (t, xn)
        return nu_c -- Collision frequency.
      end,

      useImplicitCollisionScheme = true
    },

    sourceBGK = {
      rateProfile = function (t, xn)
        return nu_Q -- Rate of the temperature source.
      end,
      tempShape = function (t, xn)
        return T_Q_cool -- Source temperature.
      end,
      useImplicitScheme = true
    },

    correct = {
      correctAllMoments = true,
      iterationEpsilon = 1.0e-12,
      maxIterations = 100,
      useLastConverged = true
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2ij }
  },

  skipField = true
}

-- Run application.
vlasovApp:run()
