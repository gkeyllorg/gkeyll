-- Propagation into a plasma wave beach with a DG 5-moment (Euler) electron fluid and the Maxwell
-- field of the Vlasov app (1x), the counterpart of moments/creg/rt_5m_mom_beach.c. The setup
-- follows entry JE8 of Ammar's Simulation Journal
-- (https://ammar-hakim.org/sj/je/je8/je8-plasmabeach.html), adapted from Section III. A. of D. N.
-- Smithe (2007), "Finite-difference time-domain simulation of fusion plasmas at radiofrequency time
-- scales", Physics of Plasmas 14 (5): 056104: the plasma frequency falls as omega_p deltaT = 25 (1
-- - x)^5 with deltaT = 0.01 / c, and an antenna current J_y = -J_0 sin(omega_drive t) on 1 - 1/400
-- < x < 1 drives a wave with omega_drive deltaT = pi / 10 toward the cutoff omega_p = omega_drive
-- at x = 0.583. The upper boundary of the field, next to the antenna, is absorbing
-- (GKYL_FIELD_ABSORB). The fluid-field coupling is the staggered operator split with the implicit
-- fluid-Maxwell source step, at omega_p dt up to 3.75.
--
-- Figures of merit at t = 5e-9 (200 cells, serendipity p1, 1000 steps): on 0.3 < x < 0.95 the cell
-- averages of E_y are the moments finite-volume solution (400 cells) times a constant to 0.013 in
-- L2 (0.002 with 400 cells, 0.023 with 100 cells where omega_p dt reaches 7.5). The constant is -1
-- / 6.79: the moments test has a copy boundary next to the antenna, and its wave has a different
-- amplitude and the opposite sign. The maximum of |E_y| is 8.70e-13, 1.85 times the amplitude mu_0
-- c J_0 / 800 the antenna radiates to each side in free space, and |E_y| stays below 0.9% of the
-- maximum for x < 0.483, beyond the cutoff. With a copy boundary in place of the absorbing one
-- (central flux of the field) the maximum is 12.6 times larger, |E_y| reaches 4.9% of it beyond the
-- cutoff and the residual of the shape is 0.60. With 100 cells the antenna covers a quarter of the
-- last cell and the two-point projection of its current doubles the amplitude (3.79 times the
-- free-space one); the shape is not affected.

local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi -- Pi.

-- Physical constants (using non-normalized physical units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.
epsilon0 = 8.854187817620389850536563031710750260608e-12 -- Permittivity of free space.
mu0 = 12.56637061435917295385057353311801153679e-7 -- Permeability of free space.
mass_elc = 9.10938215e-31 -- Electron mass.
charge_elc = -1.602176487e-19 -- Electron charge.

J0 = 1.0e-12 -- Reference current density (Amps / m^3).

-- Derived physical quantities (using non-normalized physical units).
light_speed = 1.0 / math.sqrt(mu0 * epsilon0) -- Speed of light.

-- Simulation parameters.
Nx = 200 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
Lx100 = Lx / 100.0 -- Domain size over 100 (x-direction).
x_last_edge = Lx - Lx / 400.0 -- Lower edge of the antenna, which reaches the upper boundary (low density side).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

deltaT = Lx100 / light_speed -- Arbitrary constant, with units of time.
factor = deltaT * deltaT * charge_elc * charge_elc / (mass_elc * epsilon0) -- Numerical factor for calculation of electron number density.
omega_drive = pi / 10.0 / deltaT -- Drive current angular frequency.

t_end = 5.0e-9 -- Final simulation time.
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
  periodicDirs = { }, -- Periodic directions.

  elc = Vlasov.FluidSpecies.new {
    charge = charge_elc, mass = mass_elc,
    equation = Euler.new { gasGamma = gas_gamma },

    init = function (t, xn)
      local x = xn[1]

      -- Plasma frequency profile (times deltaT).
      local omegaPdt = 25.0 * (1.0 - x) * (1.0 - x) * (1.0 - x) * (1.0 - x) * (1.0 - x)
      local ne = omegaPdt * omegaPdt / factor -- Electron number density.

      return mass_elc * ne, 0.0, 0.0, 0.0, ne * (-charge_elc) / (gas_gamma - 1.0)
    end,

    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy },
  },

  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    init = function (t, xn)
      return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    appliedCurrent = function (t, xn)
      local x = xn[1]

      -- Applied current (y-direction): the antenna next to the upper boundary.
      local app_y = 0.0
      if x > x_last_edge then
        app_y = -J0 * math.sin(omega_drive * t)
      end

      return 0.0, app_y, 0.0
    end,
    evolveAppliedCurrent = true,

    bcx = { G0.FieldBc.bcCopy, G0.FieldBc.bcAbsorb },
  },
}

vlasovApp:run()
