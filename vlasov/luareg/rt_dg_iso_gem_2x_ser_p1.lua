-- Geospace Environment Modeling (GEM) magnetic reconnection challenge with two DG isothermal Euler
-- fluids and the Maxwell field of the Vlasov app (2x), the counterpart of
-- moments/creg/rt_iso_gem.c. The initial conditions follow J. Birn et al. (2001), "Geospace
-- Environmental Modeling (GEM) Magnetic Reconnection Challenge", Journal of Geophysical Research:
-- Space Physics 106 (A3): 3715-3719: a Harris sheet of half-width 0.5 d_i with background density
-- 0.2, T_i / T_e = 5, m_i / m_e = 25, c = 10 v_A, and a flux perturbation of 0.1 B_0 d_i; periodic
-- in x, with reflecting walls for the fluids and perfectly conducting walls for the field in y.
-- Each fluid keeps the temperature of its initial condition (p = v_t^2 rho). The fluid-field
-- coupling is the staggered operator split with the implicit fluid-Maxwell source step, in which
-- the isothermal fluids carry no energy. The fluids use the Lax flux of the isothermal wave
-- equation object (rp_type) with the characteristic slope limiter (limit_fluid) and the vertex
-- positivity limiter of the density, and the field its slope limiter (limit_em) and divergence
-- cleaning of the magnetic field (mgnErrorSpeedFactor = 1).
--
-- Figures of merit at t = 25 / Omega_ci (64 x 32 cells, serendipity p1, 4167 steps): the
-- reconnected flux int |B_y(x, y = 0)| dx / 2 grows from 0.0200 to 0.1566 (the moments finite-
-- volume solution on the same cells: 0.0200 to 0.1265). The out-of-plane field reaches |B_z| =
-- 0.027 and the electron speed 0.20 (cell averages); the cell averages of the correction potential
-- stay below 2.4e-3, the densities between 0.11 and 2.4, and the field energy falls from 1.511 to
-- 1.052. The run amplifies round-off, so the late-time figures agree between runs only to a few
-- percent. B_y changes sign once along the midplane, as in the moments solution. The mirror
-- symmetry about the midplane y = 0, kept to round-off by the moments solution, is lost to 1.6% of
-- the maximum of the electron z-momentum.
local Vlasov = G0.Vlasov
local IsoEuler = G0.Vlasov.Eq.IsoEuler

-- Mathematical constants (dimensionless).
pi = math.pi -- Pi.

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_ion = 1.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.
mass_elc = 1.0 / 25.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

di = 1.0 -- Ion skin depth.
n0 = 1.0 -- Reference in-plane number density.
nb_over_n0 = 0.2 -- Background number density / reference number density.
B0 = 0.1 -- Reference magnetic field strength.
Ti_over_Te = 5.0 -- Ion temperature / electron temperature.
beta = 1.0 -- Total plasma beta.

lambda = 0.5 * di -- Current sheet width.
psi0 = 0.1 * B0 * di -- Reference magnetic scalar potential.

-- Derived physical quantities (using normalized code units).
Ti_frac = Ti_over_Te / (1.0 + Ti_over_Te) -- Fraction of total temperature from ions.
Te_frac = 1.0 / (1.0 + Ti_over_Te) -- Fraction of total temperature from electrons.
T_tot = beta * (B0 * B0) / 2.0 / n0 -- Total temperature.
omega_ci = math.abs(charge_ion * B0 / mass_ion) -- Ion cyclotron frequency.

-- Simulation parameters.
Nx = 64 -- Cell count (x-direction).
Ny = 32 -- Cell count (y-direction).
Lx = 25.6 * di -- Domain size (x-direction).
Ly = 12.8 * di -- Domain size (y-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 25.0 / omega_ci -- Final simulation time.
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
  lower = { -0.5 * Lx, -0.5 * Ly },
  upper = { 0.5 * Lx, 0.5 * Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions.

  elc = Vlasov.FluidSpecies.new {
    charge = charge_elc, mass = mass_elc,
    limitFluid = true, -- Characteristic slope limiter (before the positivity limiter).
    equation = IsoEuler.new {
      vThermal = math.sqrt(T_tot * Te_frac / mass_elc),
      rpType = G0.IsoEulerRP.Lax
    },

    init = function (t, xn)
      local y = xn[2]

      local sech_sq = (1.0 / math.cosh(y / lambda)) * (1.0 / math.cosh(y / lambda)) -- Hyperbolic secant squared.

      local n = n0 * (sech_sq + nb_over_n0) -- Total number density.
      local Jz = -(B0 / lambda) * sech_sq -- Total current density (z-direction).

      local rho = n * mass_elc -- Mass density.
      local mom_z = (mass_elc / charge_elc) * Jz * Te_frac -- Momentum density (z-direction).

      return rho, 0.0, 0.0, mom_z
    end,

    bcy = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcReflect },
  },

  ion = Vlasov.FluidSpecies.new {
    charge = charge_ion, mass = mass_ion,
    limitFluid = true, -- Characteristic slope limiter (before the positivity limiter).
    equation = IsoEuler.new {
      vThermal = math.sqrt(T_tot * Ti_frac / mass_ion),
      rpType = G0.IsoEulerRP.Lax
    },

    init = function (t, xn)
      local y = xn[2]

      local sech_sq = (1.0 / math.cosh(y / lambda)) * (1.0 / math.cosh(y / lambda)) -- Hyperbolic secant squared.

      local n = n0 * (sech_sq + nb_over_n0) -- Total number density.
      local Jz = -(B0 / lambda) * sech_sq -- Total current density (z-direction).

      local rho = n * mass_ion -- Mass density.
      local mom_z = (mass_ion / charge_ion) * Jz * Ti_frac -- Momentum density (z-direction).

      return rho, 0.0, 0.0, mom_z
    end,

    bcy = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcReflect },
  },

  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,
    mgnErrorSpeedFactor = 1.0,
    limitField = true,

    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local Bxb = B0 * math.tanh(y / lambda) -- Magnetic field of the current sheet.

      -- Current sheet plus the perturbation.
      local Bx = Bxb - psi0 * (pi / Ly) * math.cos(2.0 * pi * x / Lx) * math.sin(pi * y / Ly)
      local By = psi0 * (2.0 * pi / Lx) * math.sin(2.0 * pi * x / Lx) * math.cos(pi * y / Ly)

      return 0.0, 0.0, 0.0, Bx, By, 0.0, 0.0, 0.0
    end,

    bcy = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall },
  },
}

vlasovApp:run()
