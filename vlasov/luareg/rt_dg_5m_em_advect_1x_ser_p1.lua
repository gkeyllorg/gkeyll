-- Advection in specified electromagnetic fields with the DG 5-moment (Euler) fluid of the Vlasov
-- app (1x), the counterpart of moments/creg/rt_5m_em_advect.c and
-- vlasov/creg/rt_vlasov_em_advect_1x3v.c. The initial conditions follow entry JE32 of Ammar's
-- Simulation Journal (https://ammar-hakim.org/sj/je/je32/je32-vlasov-test-ptcl.html), rotated so
-- that the oscillating electric field E_z = cos(omega t) is along z and the background magnetic
-- field B_0 = 1 along x; omega = 0.5 Omega_c with Omega_c = q B_0 / m the cyclotron frequency, the
-- non-resonant case. The fields are external (the field itself is static) and the fluid does not
-- feed back on them. The fluid-field coupling is the implicit source step of the staggered operator
-- split (implicit midpoint at the nodes, the external field evaluated at the midpoint of each
-- source interval).
--
-- The uniform fluid follows the exact single-particle solution u_y = A (cos(omega t) - cos(Omega_c
-- t)), u_z = A (sin(Omega_c t) - (omega / Omega_c) sin(omega t)) with A = Omega_c (q / m) E_0 /
-- (Omega_c^2 - omega^2): at t = 100, u_y = 0.136863 and u_z = 0.500238. Figures of merit (2 cells,
-- serendipity p1, 3945 steps, Omega_c dt = 0.025): u_y = 0.140498 and u_z = 0.506395, errors 3.6e-3
-- and 6.2e-3. The internal energy and the density stay at their initial values to round-off
-- (1e-13), and the state stays uniform.

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
omega = 0.5 -- Oscillating field frequency normalized to cyclotron frequency.

-- Simulation parameters.
Nx = 2 -- Cell count (configuration space: x-direction).
Lx = 4.0 * pi -- Domain size (configuration space: x-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 / 64.0 -- CFL coefficient.

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
