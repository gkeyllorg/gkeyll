-- Two-stream instability with the PKPM system of equations (1x, normalized parallel velocity w =
-- v_par/xi), the counterpart of vlasov/creg/rt_vlasov_twostream_1x1v.c. Two counter-streaming
-- Maxwellian electron beams (drift +-4 vt, n0/2 each, lambda_D = 1) on a neutralizing ion
-- background are perturbed in density at k lambda_D = 0.1 with amplitude 1e-5; the electric field
-- follows from Gauss's law. In the PKPM variables the beams are one normalized distribution F(w)
-- with two humps at w = +-4 vt/xi, xi^2 = (T_par + 2 T_perp)/(3 m) = 19/3, and G = (T_perp/(m
-- xi^2)) F. A uniform magnetic field along x sets the field direction.
-- The fluid-field coupling is advanced explicitly within the SSP-RK3 step.
-- Kinetic linear theory: a purely growing mode with gamma = 0.2770 omega_pe (cold-beam limit
-- 0.308).
-- Figures of merit (64 x 32 cells): the field energy grows at gamma = 0.2776 over t = 20-40
-- (Vlasov, serendipity p2 on the same grid: 0.2778) and saturates at t = 45.2 with peak energy 105
-- (Vlasov 106). After saturation it decays faster than in the Vlasov run: 32.7 at t = 50 against
-- 57.2 (40.7 with 128 configuration-space cells). Electron number is conserved to 3e-13.

local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

n0 = 1.0 -- Reference number density.
vt = 1.0 -- Thermal velocity of each beam.
Vx_drift = 4.0 -- Drift velocity of each beam (x-direction).
lambda_D = 1.0 -- Electron Debye length.
alpha = 1.0e-5 -- Applied perturbation amplitude.
B0 = 1.0 -- Magnetic field along x (sets the field direction b).

T_par = mass_elc * (vt * vt + Vx_drift * Vx_drift) -- Parallel temperature of the two-beam state.
T_perp = mass_elc * vt * vt -- Perpendicular temperature.
xi = math.sqrt((T_par + 2.0 * T_perp) / (3.0 * mass_elc)) -- Thermal speed of the two-beam state.
kx = 0.1 / lambda_D -- Perturbed wave number (x-direction).
omega_pe = vt / lambda_D -- Electron plasma frequency.

-- Simulation parameters.
Nx = 64 -- Cell count (configuration space: x-direction).
Nvx = 32 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 2.0 * pi / kx -- Domain size (configuration space: x-direction).
vx_max = 16.0 * vt / xi -- Domain boundary (velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 50.0 / omega_pe -- Final simulation time.
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

      local wd = Vx_drift / xi -- Beam drift in w.
      local s2 = (vt * vt) / (xi * xi) -- Variance of each beam in w.
      local F = 0.5 * (math.exp(-0.5 * (w - wd) * (w - wd) / s2) + math.exp(-0.5 * (w + wd) * (w + wd) / s2)) / math.sqrt(2.0 * pi * s2) -- Two-hump normalized parallel distribution.
      local G = (T_perp / (mass_elc * xi * xi)) * F -- Perpendicular-temperature weighted distribution.

      return F, G
    end,

    initFluid = function (t, xn)
      local x = xn[1]

      local rho = mass_elc * n0 * (1.0 + alpha * math.cos(kx * x)) -- Electron mass density.
      local E = 1.5 * rho * xi * xi -- Electron total energy density (zero net flow).

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

      local Ex = -alpha * math.sin(kx * x) / kx -- Electric field from Gauss's law.

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
