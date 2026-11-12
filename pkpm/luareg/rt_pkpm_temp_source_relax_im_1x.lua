-- Relaxation of the temperature and of the pressure anisotropy under the BGK-form temperature
-- source nu_Q (f_M(n, u, T_Q) - f) with the PKPM system of equations (1x, normalized parallel
-- velocity w = v_par/xi), implicit update: the backward-Euler twin of
-- rt_pkpm_temp_source_relax_1x.c at 100 times its source rate.
-- Two homogeneous neutral species start as drifting bi-Maxwellians at T = 1 with T_par/T = 1.4 and
-- T_perp/T = 0.8. The source of the first has T_Q = 10 (heating, T_Q/T = 10 initially) and that of
-- the second T_Q = 0.1 (cooling). A static magnetic field along x sets the field direction. Source
-- rate 100, with nu_Q dt between 1.5 and 1.7 (the explicit update would need nu_Q dt below 0.1).
-- The source takes one backward-Euler step in each SSP-RK3 stage, over the time weight of the
-- stage, so with z = nu_Q dt a step multiplies a relaxing quantity by
-- R(z) = [1/3 + (2/3) (3/4 + (1/4)/(1 + z))/(1 + z/4)]/(1 + 2 z/3) (0.335 at z = 1.74).
-- The stages combine the distributions with weights that follow the temperature change of the
-- implicit steps, so the combination is linear in the pressure components although the temperature
-- changes seven-fold in the first step of the heated species.
-- Density and momentum are conserved to round-off and the constraints stay at [1, 0, 1]. Each step
-- multiplies T - T_Q by R(z) to 1e-9 and the pressure anisotropy p_par - p_perp by R(z) to 1e-6, and
-- the temperature stays between T_0 and T_Q: 6.9882 and 0.40118 after the first step.

local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass = 1.0 -- Species mass.
charge = 0.0 -- Species charge.

n0 = 1.0 -- Reference number density.
T0 = 1.0 -- Initial temperature.
Tpar_over_T = 1.4 -- Initial parallel temperature / temperature.
ux0 = 0.5 -- Flow velocity (x-direction, along the field).
uy0 = 0.3 -- Flow velocity (y-direction).
nu_Q = 100.0 -- Rate of the temperature source.
T_heat = 10.0 -- Temperature of the source of the heated species.
T_cool = 0.1 -- Temperature of the source of the cooled species.
B0 = 1.0 -- Magnetic field along x (sets the field direction b).

-- Derived physical quantities (using normalized code units).
Tperp_over_T = 0.5 * (3.0 - Tpar_over_T) -- Initial perpendicular temperature / temperature.

-- Simulation parameters.
Nx = 2 -- Cell count (configuration space: x-direction).
Nvx = 48 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 5.0 -- Domain size (configuration space: x-direction).
vx_max = 8.0 -- Domain boundary (velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.5 -- Final simulation time.
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

  -- Heated species.
  heat = PKPM.Species.new {
    charge = charge, mass = mass,

    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[2]

      -- Bi-Maxwellian in w = v_par/xi: F has variance T_par/T and G = (T_perp/T) F.
      local F = math.exp(-w * w / (2.0 * Tpar_over_T)) / math.sqrt(2.0 * pi * Tpar_over_T)

      return F, Tperp_over_T * F
    end,

    initFluid = function (t, xn)
      local rho = mass * n0 -- Mass density.
      local E = 0.5 * rho * (ux0 * ux0 + uy0 * uy0) + 1.5 * n0 * T0 -- Total energy density.

      return rho, rho * ux0, rho * uy0, 0.0, E
    end,

    sourceBGK = {
      rateProfile = function (t, xn)
        return nu_Q
      end,

      tempShape = function (t, xn)
        return T_heat
      end,

      useImplicitScheme = true,
    },

    evolve = true -- Evolve species?
  },

  -- Cooled species.
  cool = PKPM.Species.new {
    charge = charge, mass = mass,

    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[2]

      -- Bi-Maxwellian in w = v_par/xi: F has variance T_par/T and G = (T_perp/T) F.
      local F = math.exp(-w * w / (2.0 * Tpar_over_T)) / math.sqrt(2.0 * pi * Tpar_over_T)

      return F, Tperp_over_T * F
    end,

    initFluid = function (t, xn)
      local rho = mass * n0 -- Mass density.
      local E = 0.5 * rho * (ux0 * ux0 + uy0 * uy0) + 1.5 * n0 * T0 -- Total energy density.

      return rho, rho * ux0, rho * uy0, 0.0, E
    end,

    sourceBGK = {
      rateProfile = function (t, xn)
        return nu_Q
      end,

      tempShape = function (t, xn)
        return T_cool
      end,

      useImplicitScheme = true,
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
