-- Traveling electromagnetic plane wave with the Maxwell solver (1x). A unit-amplitude plane wave
-- with two wavelengths per box length along x is advanced for five periods on a periodic domain,
-- after which it must coincide with its initial state. After five periods the field returns to its
-- initial state to 8e-5 in the relative L2 norm, with the energy conserved to 2e-4.
local Vlasov = G0.Vlasov

pi = math.pi -- Pi.
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
E0 = 1.0 / math.sqrt(3.0) -- Electric field amplitude (|E| = 1).

k_wave_x = 2.0 -- Wave number (x-direction, wavelengths per box length).
omega = 2.0 * math.pi * math.sqrt(k_wave_x * k_wave_x) -- Wave frequency (c = 1).
t_period = 2.0 * math.pi / omega -- Wave period.

Nx = 64 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
poly_order = 3 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.8 -- CFL coefficient.

t_end = 5.0 * t_period -- Final simulation time.
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
  periodicDirs = { 1 }, -- Periodic directions.

  -- Field (Maxwell).
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0,

    -- Initial conditions function (E, B, phi, psi).
    init = function (t, xn)
      local x = xn[1]

      local phase = (2.0 * pi / Lx) * (k_wave_x * x)
      local amp = 1.0 * math.cos(phase) -- Traveling plane wave.

      local Ex = 0.0 * amp
      local Ey = 1.0 / math.sqrt(2.0) * amp
      local Ez = 1.0 / math.sqrt(2.0) * amp
      local Bx = 0.0 * amp
      local By = -1.0 / math.sqrt(2.0) * amp
      local Bz = 1.0 / math.sqrt(2.0) * amp

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

  }
}

vlasovApp:run()
