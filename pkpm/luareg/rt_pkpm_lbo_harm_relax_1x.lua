-- Relaxation of first-harmonic perturbations under LBO self-collisions with the PKPM system of
-- equations (1x, normalized parallel velocity w = v_par/xi), the counterpart of
-- vlasov/creg/rt_vlasov_lbo_harm_relax_1x3v.c.
-- Two homogeneous neutral species start from a unit Maxwellian, F = G = M(w), with a small first
-- gyroangle harmonic about a static magnetic field along x. Species drift carries M0_y = a w M and
-- N_y = 2 M0_y, the reduced moments of f = f_M (1 + a v_par v_y). Species heat carries only
-- N_y = 2 a M, the reduced moments of f = f_M (1 + a v_y (v_perp^2/2 - 2)), a perpendicular heat
-- flux with no perpendicular drift. Collision frequency 1, run for 2 collision times.
-- int M0_perp dw = 0 holds to 2e-17 and the other constraints to 4e-14. The LBO relaxes the Hermite
-- mode k (in w) of M0_perp at (k + 1) nu and that of N_perp - 2 M0_perp at (k + 3) nu: int w M0_y dw
-- and int w N_y dw (k = 1, with N_y = 2 M0_y) decay at rate 2.00000, and int N_y dw (k = 0, with
-- M0_y = 0) at 3.00000 (Vlasov 1x3v: the off-diagonal pressure at 2.00000 and the perpendicular heat
-- flux at 3.00006).
local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass = 1.0 -- Species mass.
charge = 0.0 -- Species charge.

n0 = 1.0 -- Reference number density.
nu = 1.0 -- Collision frequency.
a_drift = 0.1 -- Amplitude of the drift-like first harmonic.
a_heat = 0.05 -- Amplitude of the heat-flux-like first harmonic.
B0 = 1.0 -- Magnetic field along x (sets the field direction b).


-- Simulation parameters.
Nx = 2 -- Cell count (configuration space: x-direction).
Nvx = 32 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 5.0 -- Domain size (configuration space: x-direction).
vx_max = 8.0 -- Domain boundary (velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 2.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate integrated L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

pkpmApp = PKPM.App.new {

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
  decompCuts = { 1 }, -- Cuts in each coordinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Species seeded with the drift-like first harmonic.
  drift = PKPM.Species.new {
    charge = charge, mass = mass,

    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[2]
      local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * pi) -- Unit Maxwellian in w.

      return F, F
    end,

    initHarm = function (t, xn)
      local w = xn[2]
      local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * pi) -- Unit Maxwellian in w.

      local M0_y = a_drift * w * F -- Reduced moment of f = f_M (1 + a w z_y).

      return 0.0, M0_y, 0.0, 0.0, 2.0 * M0_y, 0.0
    end,

    initFluid = function (t, xn)
      local rho = mass * n0 -- Mass density.

      return rho, 0.0, 0.0, 0.0, 1.5 * rho
    end,

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu
      end
    },

    evolve = true -- Evolve species?
  },

  -- Species seeded with the heat-flux-like first harmonic.
  heat = PKPM.Species.new {
    charge = charge, mass = mass,

    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[2]
      local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * pi) -- Unit Maxwellian in w.

      return F, F
    end,

    initHarm = function (t, xn)
      local w = xn[2]
      local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * pi) -- Unit Maxwellian in w.

      -- The reduced moments of f = f_M (1 + a z_y (z^2/2 - 2)), which has no M0_perp.
      return 0.0, 0.0, 0.0, 0.0, 2.0 * a_heat * F, 0.0
    end,

    initFluid = function (t, xn)
      local rho = mass * n0 -- Mass density.

      return rho, 0.0, 0.0, 0.0, 1.5 * rho
    end,

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu
      end
    },

    evolve = true -- Evolve species?
  },

  -- Field.
  field = PKPM.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      return 0.0, 0.0, 0.0, B0, 0.0, 0.0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()
