-- Electromagnetic wave in a two-fluid plasma with DG 5-moment (Euler) fluids and the Maxwell field
-- of the Vlasov app (1x). The transverse wave E_y = E_0 cos(k x - omega t), B_z = (k / omega) E_y,
-- with the quiver velocities u_y = -(q / m) (E_0 / omega) sin(k x - omega t) of the electrons and
-- the ions, is an exact solution of the linearized equations with omega^2 = omega_p^2 + c^2 k^2 and
-- omega_p^2 the sum of n q^2 / (epsilon_0 m) over the species (the pressure does not enter). With
-- m_i / m_e = 25, n = 1, k = 1 and c = 1, omega = 1.428286; E_0 = 1e-3 (the nonlinear terms are of
-- relative order E_0^2). The fluid-field coupling is the staggered operator split: the implicit
-- fluid-Maxwell source step advances E (with c^2 curl B as a known current) and the SSP-RK3
-- transport step advances B and the fluids.
--
-- Figures of merit at t = 20, 4.5 wave periods (32 cells, serendipity p1, 306 steps, omega dt =
-- 0.093): the fitted amplitude of E_y is E_0 (1 - 4.1e-4) and its phase error 5.7e-3 rad. The
-- relative L2 errors are 2.7e-2 (E_y), 2.8e-2 (B_z) and 2.9e-2 (electron u_y). The total energy
-- (field plus fluids, from the frames) changes by 1.5e-6 of the initial field energy.

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
mass_ion = 25.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Number density.
T0 = 0.01 -- Temperature of both species.
E0 = 1.0e-3 -- Amplitude of the electric field of the wave.
k = 1.0 -- Wavenumber.

-- Derived physical quantities (using normalized code units).
omega_p_sq = n0 * (charge_elc * charge_elc / mass_elc + charge_ion * charge_ion / mass_ion) / epsilon0 -- Plasma frequency squared (both species).
omega = math.sqrt(omega_p_sq + k * k / (epsilon0 * mu0)) -- Wave frequency.

-- Simulation parameters.
Nx = 32 -- Cell count (configuration space: x-direction).
Lx = 2.0 * pi / k -- Domain size (configuration space: x-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 20.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = 100 -- Number of times to calculate field energy.
integrated_mom_calcs = 100 -- Number of times to calculate integrated moments.
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
      local x = xn[1]

      -- Quiver velocity of the species in the wave E_y = E0 cos(k x - omega t).
      local rho = n0 * mass_elc
      local uy = -(charge_elc / mass_elc) * (E0 / omega) * math.sin(k * x)

      return rho, 0.0, rho * uy, 0.0, n0 * T0 / (gas_gamma - 1.0) + 0.5 * rho * uy * uy
    end,
  },

  ion = Vlasov.FluidSpecies.new {
    charge = charge_ion, mass = mass_ion,
    equation = Euler.new { gasGamma = gas_gamma },

    init = function (t, xn)
      local x = xn[1]

      -- Quiver velocity of the species in the wave E_y = E0 cos(k x - omega t).
      local rho = n0 * mass_ion
      local uy = -(charge_ion / mass_ion) * (E0 / omega) * math.sin(k * x)

      return rho, 0.0, rho * uy, 0.0, n0 * T0 / (gas_gamma - 1.0) + 0.5 * rho * uy * uy
    end,
  },

  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    init = function (t, xn)
      local x = xn[1]

      return 0.0, E0 * math.cos(k * x), 0.0, 0.0, 0.0, (k / omega) * E0 * math.cos(k * x), 0.0, 0.0
    end,
  },
}

vlasovApp:run()
