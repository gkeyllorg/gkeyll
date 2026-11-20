-- Relaxation of the temperature and of the pressure anisotropy under the BGK-form temperature
-- source nu_Q (f_M(n, u, T_Q) - f) with the Vlasov solver (1x2v), explicit update: the counterpart
-- of rt_pkpm_temp_source_relax_1x.c.
-- Two homogeneous neutral species start as bi-Maxwellians drifting along x with T_xx = 1.4 and T_yy
-- = 0.6 (T = 1). The source of the first has T_Q = 10 (heating) and that of the second T_Q = 0.1
-- (cooling). Source rate 1, run for 5 relaxation times.
-- Figures of merit (serendipity p2, 2 x 48 x 48 cells, 201 steps): density and momentum are
-- conserved to round-off (2e-14). The temperature follows T_Q + (T_0 - T_Q) exp(-nu_Q t) to 2.2e-6
-- (heated) and 2.2e-7 (cooled): 9.9394 and 0.10606 at t = 5, with T - T_Q decaying at rate
-- 1.000001. The pressure anisotropy P_xx - P_yy decays at rates 1.000003 and 1.000001.

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
nu_Q = 1.0 -- Rate of the temperature source.
T_Q_heat = 10.0 -- Source temperature of the heated species.
T_Q_cool = 0.1 -- Source temperature of the cooled species.

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

t_end = 5.0 -- Final simulation time.
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

    sourceBGK = {
      rateProfile = function (t, xn)
        return nu_Q -- Rate of the temperature source.
      end,
      tempShape = function (t, xn)
        return T_Q_heat -- Source temperature.
      end,
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

    sourceBGK = {
      rateProfile = function (t, xn)
        return nu_Q -- Rate of the temperature source.
      end,
      tempShape = function (t, xn)
        return T_Q_cool -- Source temperature.
      end,
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
