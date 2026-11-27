-- Langmuir oscillation about a held, nonuniform electron background with the DG 5-moment (Euler)
-- equations and the Poisson field of the Vlasov app (1x). Electrons (vt = 0.01, lambda_D = 0.01,
-- no ions) with the density n0 (1 + 0.5 cos(k x)), k = 1, are held about that state at rest: the
-- potential is solved from the departure of the electrons from their background, so the
-- background, which is not neutral, exerts no force (its own field would be 5e3 times that of the
-- oscillation), and the app derives the body force that balances its pressure gradient. The
-- initial velocity u0 sin(k x), u0 = 1e-4, starts an oscillation of every fluid element at its
-- local plasma frequency: for a cold plasma u_x = u0 sin(k x) cos(omega_p(x) t) and E_x = (m /
-- |q|) u0 sin(k x) omega_p(x) sin(omega_p(x) t), with omega_p(x)^2 = omega_pe^2 (1 + 0.5 cos(k
-- x)). The oscillation phase-mixes along x: after two periods of omega_pe, the end of the run, the
-- phase differs by 6.5 radians across the domain. Thermal corrections are of the order of (k
-- lambda_D)^2. The coupling is the kick-drift-kick split, whose kicks take E = -grad(phi); the
-- step resolves the largest plasma frequency (cfl_frac = 0.2).
-- rt_vp_langmuir_background_1x1v_ser_p2.lua and rt_pkpm_poisson_langmuir_background_1x.lua are the
-- same problem with kinetic and PKPM electrons.
--
-- Figures of merit at t = 4 pi / omega_pe (64 cells, serendipity p1, 77 steps, omega_pe dt =
-- 0.163): the cell averages of E_x differ from the cold local solution by 1.1e-2 in L2, and by
-- 7.3e-3 with the leapfrog frequency (2 / dt) asin(omega_p dt / 2) of the split in it (1.7e-2 at
-- 32 cells and 3.3e-3 at 128); those of u_x differ by 1.7e-2. The density changes by at most
-- 3.1e-4 of itself, the compression of the oscillation.

local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi -- Pi.

-- Physical constants (using normalized code units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.
epsilon0 = 1.0 -- Permittivity of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

n0 = 1.0 -- Number density.
T0 = 1.0e-4 -- Electron temperature.
a_bg = 0.5 -- Relative amplitude of the density variation of the background.
u0 = 1.0e-4 -- Amplitude of the initial velocity.
k = 1.0 -- Wavenumber.

-- Derived physical quantities (using normalized code units).
omega_pe = math.sqrt(n0 * charge_elc * charge_elc / (epsilon0 * mass_elc)) -- Plasma frequency.
vt = math.sqrt(T0 / mass_elc) -- Electron thermal velocity.

-- Simulation parameters.
Nx = 64 -- Cell count (configuration space: x-direction).
Lx = 2.0 * pi / k -- Domain size (configuration space: x-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.2 -- CFL coefficient.

t_end = 4.0 * pi / omega_pe -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Background of the electrons: a density that varies along x, at rest, with a uniform temperature.
function elcBackground (t, xn)
  local x = xn[1]

  local n = n0 * (1.0 + a_bg * math.cos(k * x)) -- Number density.

  return n * mass_elc, 0.0, 0.0, 0.0, n * T0 / (gas_gamma - 1.0)
end

-- Initial condition: the background with the velocity u0 sin(k x).
function elcInit (t, xn)
  local x = xn[1]

  local rho, _, _, _, E = elcBackground(t, xn)
  local ux = u0 * math.sin(k * x) -- Flow velocity.

  return rho, rho * ux, 0.0, 0.0, E + 0.5 * rho * ux * ux
end

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

  -- Electrons.
  elc = Vlasov.FluidSpecies.new {
    charge = charge_elc, mass = mass_elc,
    equation = Euler.new { gasGamma = gas_gamma },

    init = elcInit,
    -- The electrons are held about their background, evaluated at the nodes (continuous).
    background = elcBackground,
    evalOnNodes = true,
  },

  -- Field: the electrostatic potential of the departure of the electrons from their background.
  field = Vlasov.Field.new {
    fieldID = G0.FieldModel.Phi,
    epsilon0 = epsilon0,

    poissonBcs = {
      lowerType = {
        G0.PoissonBc.bcPeriodic
      },
      upperType = {
        G0.PoissonBc.bcPeriodic
      }
    }
  }
}

vlasovApp:run()
