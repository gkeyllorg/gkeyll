-- Buneman instability with the PKPM system of equations (1x, normalized parallel velocity w =
-- v_par/xi), the counterpart of vlasov/creg/rt_vlasov_buneman_1x1v.c. Parameters follow entry JE33
-- of Ammar's Simulation Journal (https://ammar-hakim.org/sj/je/je33/je33-buneman.html): m_i/m_e =
-- 25, electron drift 0.159 so that k v_d = omega_pe for the box mode k = 2 pi, vte = 0.02, vti =
-- 0.001, density perturbation 1e-5 with the electric field from Gauss's law. The ghost current
-- cancels the net current of the drifting electrons, so that Ampere's law equals the periodic
-- Poisson solve. In the PKPM variables the drift is the electron fluid velocity and F is a
-- Maxwellian in w; the thermal-speed normalization follows the heating, so the w extents stay at
-- +-8. A uniform magnetic field along x sets the field direction.
-- The fluid-field coupling is advanced explicitly within the SSP-RK3 step.
-- Kinetic linear theory at k = 2 pi: gamma = 0.1976 omega_pe, omega_r = 0.1425 omega_pe (cold-beam
-- maximum 0.235).
-- Figures of merit (16 x 64 cells): the field energy grows at gamma = 0.1881 over t = 20-65
-- (Vlasov, serendipity p2 on the same grid: 0.1950). It saturates near t = 72 at about a quarter of
-- the Vlasov peak (1.6e-3 at t = 72.3 against 7.0e-3 at t = 72.1) and is 1.4e-3 at t = 80 (Vlasov
-- 1.1e-3). Electron and ion number are conserved to 1e-12.

local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 25.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.
Vx_drift_elc = 0.159 -- Electron drift velocity (x-direction).
Vx_drift_ion = 0.0 -- Ion drift velocity (x-direction).
vte = 0.02 -- Electron thermal velocity.
vti = 0.001 -- Ion thermal velocity.
alpha = 1.0e-5 -- Applied perturbation amplitude.
B0 = 1.0 -- Magnetic field along x (sets the field direction b).

Lx = 1.0 -- Domain size (configuration space: x-direction).
kx = 2.0 * pi / Lx -- Perturbed wave number (x-direction).
omega_pe = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_elc)) -- Electron plasma frequency.

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Nvx = 64 -- Cell count (velocity space: normalized parallel velocity w).
vx_max_elc = 8.0 -- Domain boundary (electron velocity space: w = v_par/xi).
vx_max_ion = 8.0 -- Domain boundary (ion velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 80.0 / omega_pe -- Final simulation time.
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
    lower = { -vx_max_elc },
    upper = { vx_max_elc },
    cells = { Nvx },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[2]

      local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * pi) -- Normalized Maxwellian in w = v_par/xi.

      return F, F
    end,

    initFluid = function (t, xn)
      local x = xn[1]

      local rho = mass_elc * n0 * (1.0 + alpha * math.cos(kx * x)) -- Electron mass density.
      local ux = Vx_drift_elc -- Electron drift.
      local E = 0.5 * rho * ux * ux + 1.5 * rho * vte * vte -- Electron total energy density.

      return rho, rho * ux, 0.0, 0.0, E
    end,

    evolve = true -- Evolve species?
  },

  -- Ions.
  ion = PKPM.Species.new {
    charge = charge_ion, mass = mass_ion,

    -- Velocity space grid.
    lower = { -vx_max_ion },
    upper = { vx_max_ion },
    cells = { Nvx },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[2]

      local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * pi) -- Normalized Maxwellian in w = v_par/xi.

      return F, F
    end,

    initFluid = function (t, xn)
      local rho = mass_ion * n0 -- Ion mass density.
      local ux = Vx_drift_ion -- Ion drift.
      local E = 0.5 * rho * ux * ux + 1.5 * rho * vti * vti -- Ion total energy density.

      return rho, rho * ux, 0.0, 0.0, E
    end,

    evolve = true -- Evolve species?
  },

  -- Field.
  field = PKPM.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local Ex = -alpha * math.sin(kx * x) / kx -- Electric field from Gauss's law.

      return Ex, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    externalFieldInit = function (t, xn)
      return 0.0, 0.0, 0.0, B0, 0.0, 0.0
    end,

    useGhostCurrent = true,

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()
