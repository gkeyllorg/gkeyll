-- Advection in specified electromagnetic fields with the DG 5-moment (Euler) fluid of the Vlasov
-- app (1x), the counterpart of moments/creg/rt_5m_em_advect_resonant.c and
-- vlasov/creg/rt_vlasov_em_advect_resonant_1x3v.c. The initial conditions follow entry JE32 of
-- Ammar's Simulation Journal (https://ammar-hakim.org/sj/je/je32/je32-vlasov-test-ptcl.html),
-- rotated so that the oscillating electric field E_z = cos(omega t) is along z and the background
-- magnetic field B_0 = 1 along x; omega = Omega_c with Omega_c = q B_0 / m the cyclotron frequency,
-- the resonant case. The fields are external (the field itself is static) and the fluid does not
-- feed back on them. The fluid-field coupling is explicit (use_explicit_sources): the Lorentz force
-- is in the fluid RHS of every SSP-RK3 stage, and the time-step resolves the rate of that source,
-- the larger of the cyclotron frequency and twice the acceleration over the thermal speed sqrt(2 e
-- / rho).
--
-- The uniform fluid follows the exact single-particle solution u_y = (q / m) Omega_c t sin(omega t)
-- / (2 omega), u_z = (q / m) (sin(omega t) / omega + t cos(omega t)) / 2, whose amplitude grows
-- linearly in time: at t = 100, u_y = -25.318282 and u_z = -42.862761. Figures of merit (2 cells,
-- serendipity p1, cfl_frac = 0.1, 28798 steps): u_y = -25.318174 and u_z = -42.862721, errors
-- 1.1e-4 and 4.0e-5 (2.3e-6 of the speed, 49.8). The time-step falls from 0.087 to 0.002 as the
-- gyrating fluid speeds up (its acceleration grows with its speed while its thermal speed does
-- not). The internal energy rises by 1.2e-4 of its initial value (the implicit variant,
-- rt_dg_5m_em_advect_resonant_1x, keeps it to 2.1e-10), and the state stays uniform.

local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi -- Pi.

-- Physical constants (using normalized code units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.
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
Lx = 4.0 * pi -- Domain size (configuration space: x-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.1 -- CFL coefficient.

t_end = 100.0 -- Final simulation time.
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

  -- Fluid-field coupling in the SSP-RK3 stages instead of the implicit operator split.
  useExplicitSources = true,

  elc = Vlasov.FluidSpecies.new {
    charge = charge_elc, mass = mass_elc,
    equation = Euler.new { gasGamma = gas_gamma },

    init = function (t, xn)
      local T = vt * vt * mass_elc
      local rho_elc = n0 * mass_elc -- Electron mass density.
      local E_elc = n0 * T / (gas_gamma - 1.0) -- Electron total energy density.

      return rho_elc, 0.0, 0.0, 0.0, E_elc
    end,
  },

  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    init = function (t, xn)
      return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    externalFieldInit = function (t, xn)
      return 0.0, 0.0, math.cos(omega * t), B0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    evolveExternalField = true,
  },
}

vlasovApp:run()
