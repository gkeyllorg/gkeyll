-- Advection in specified electromagnetic fields with the PKPM system of equations (1x, normalized
-- parallel velocity w = v_par/xi), the counterpart of
-- vlasov/creg/rt_vlasov_em_advect_resonant_1x3v.c. The initial conditions follow entry JE32 of
-- Ammar's Simulation Journal (https://ammar-hakim.org/sj/je/je32/je32-vlasov-test-ptcl.html),
-- rotated so that the oscillating electric field E_z = cos(omega t) is along z and the background
-- magnetic field B_0 = 1 along x; omega = 1.0 Omega_c with Omega_c = q B_0 / m the cyclotron
-- frequency, the resonant case. The fields are external and the species does not feed back on them.
-- The fluid-field coupling is advanced explicitly within the SSP-RK3 step.
-- The resonant drive grows the bulk velocity secularly, u ~ (q / m) E_0 t / 2: at t = 10, u_y =
-- -2.72010 and u_z = 4.46736.
-- Figures of merit (2 x 8 cells, Omega_c dt = 0.12): at t = 10, u_y = -2.71977 and u_z = 4.46565
-- (errors +3.3e-4 and -1.7e-3). SSP-RK3 slightly damps the gyration and the bulk kinetic energy it
-- removes appears as heat: the temperature ends at 1.0035. The Vlasov test is exact to 1e-5.

local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

vt = 1.0 -- Thermal velocity.
n0 = 1.0 -- Reference number density.
B0 = 1.0 -- Reference magnetic field strength.
omega = 1.0 -- Oscillating field frequency normalized to cyclotron frequency.

-- Simulation parameters.
Nx = 2 -- Cell count (configuration space: x-direction).
Nvx = 8 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 4.0 * pi -- Domain size (configuration space: x-direction).
vx_max = 8.0 -- Domain boundary (velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 10.0 -- Final simulation time.
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
      local rho = mass_elc * n0 -- Electron mass density.
      local E = 1.5 * n0 * mass_elc * vt * vt -- Electron total energy density (at rest).

      return rho, 0.0, 0.0, 0.0, E
    end,

    evolve = true -- Evolve species?
  },

  -- Field.
  field = PKPM.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    externalFieldInit = function (t, xn)
      local Ez = math.cos(omega * t) -- External electric field (z-direction).

      return 0.0, 0.0, Ez, B0, 0.0, 0.0
    end,
    evolveExternalField = true,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()
