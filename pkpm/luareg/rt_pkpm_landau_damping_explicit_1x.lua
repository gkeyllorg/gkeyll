-- Linear Landau damping of a Langmuir wave with the PKPM system of equations (1x, normalized
-- parallel velocity w = v_par/xi), the counterpart of vlasov/creg/rt_vlasov_landau_damping_1x1v.c.
-- A Maxwellian electron plasma (vt = 1, lambda_D = 1) on a neutralizing ion background is perturbed
-- in density at k lambda_D = 0.5 with amplitude 1e-4; the electric field follows from Gauss's law.
-- A uniform magnetic field along x gives PKPM its field direction and does not act on the parallel
-- electrostatic wave.
-- The fluid-field coupling is advanced explicitly within the SSP-RK3 step.
-- Collisionless linear theory at k lambda_D = 0.5 (Vlasov; the parallel kinetics of PKPM are exact
-- for this wave): omega = 1.4157 omega_pe, damping rate gamma = 0.1533 omega_pe.
-- Figures of merit (32 x 32 cells): the field-energy peaks over t = 0.5-12 decay at gamma = 0.1550
-- and recur at omega = 1.414; the field energy at t = 20 is 0.0009 of its initial value, as with
-- the implicit coupling (Vlasov, serendipity p2 on the same grid: 0.1549, 1.413, 0.00091).

local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

vt = 1.0 -- Thermal velocity.
B0 = 1.0 -- Magnetic field along x (sets the field direction b).
alpha = 1.0e-4 -- Applied perturbation amplitude.
k0 = 0.5 -- Perturbed wave number.

-- Simulation parameters.
Nx = 32 -- Cell count (configuration space: x-direction).
Nvx = 32 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 4.0 * pi -- Domain size (configuration space: x-direction).
vx_max = 6.0 -- Domain boundary (velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 20.0 -- Final simulation time.
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
  lower = { -0.5 * Lx },
  upper = { 0.5 * Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  useExplicitSource = true,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coordinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Electrons.
  elc = PKPM.Species.new {
    charge = charge_elc, mass = mass_elc,

    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[2]

      local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * pi) -- Normalized Maxwellian in w = v_par/xi.

      return F, F
    end,

    initFluid = function (t, xn)
      local x = xn[1]

      local n = 1.0 + alpha * math.cos(k0 * x) -- Electron number density.
      local rho = mass_elc * n -- Electron mass density.
      local E = 1.5 * n * mass_elc * vt * vt -- Electron total energy density (at rest).

      return rho, 0.0, 0.0, 0.0, E
    end,

    evolve = true -- Evolve species?
  },

  -- Field.
  field = PKPM.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local Ex = -alpha * math.sin(k0 * x) / k0 -- Electric field (x-direction).

      return Ex, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    externalFieldInit = function (t, xn)
      return 0.0, 0.0, 0.0, B0, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()
