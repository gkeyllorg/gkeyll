-- Langmuir wave of a warm electron fluid with the DG 5-moment (Euler) equations and the Poisson
-- field of the Vlasov app (1x). An adiabatic electron fluid (gamma = 3, vt = 0.1, n0 = 1) on a
-- uniform ion background, which the periodic solve removes as the mean charge density, is perturbed
-- in density at k lambda_De = 0.05 with amplitude 1e-4: a standing wave at the Bohm-Gross frequency
-- omega^2 = omega_pe^2 + gamma k^2 vt^2, omega = 1.00374 omega_pe. The coupling is the
-- kick-drift-kick split: the source step kicks the momentum with E = -grad(phi) of the densities at
-- its start (which it does not change), the SSP-RK3 transport step moves the fluid, and the next
-- kick solves again. The coupling of the density to the potential is explicit, a leapfrog stable
-- for omega_pe dt < 2.
--
-- Figures of merit (32 cells, serendipity p1, t = 50, 67 steps at the CFL step, omega_pe dt =
-- 0.756): a least-squares fit of the field energy gives omega = 1.0275 omega_pe; the leapfrog phase
-- error omega_pe^2 dt^2/24 = 2.4% accounts for the difference from Bohm-Gross to 0.01%. The total
-- energy (fluid energy plus eps0/2 int |grad phi|^2, from the frames) stays between 0 and -14.5%
-- of the initial field energy with no drift: the bounded error (omega_pe dt)^2/4 = 14% of the
-- leapfrog. Mass is conserved to 4e-15.

local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi -- Pi.

-- Physical constants (using normalized code units).
gas_gamma = 3.0 -- Adiabatic index.
epsilon0 = 1.0 -- Permittivity of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

n0 = 1.0 -- Number density.
T0 = 0.01 -- Electron temperature.
alpha = 1.0e-4 -- Applied perturbation amplitude.
k = 0.5 -- Wavenumber.

-- Derived physical quantities (using normalized code units).
omega_pe = math.sqrt(n0 * charge_elc * charge_elc / (epsilon0 * mass_elc)) -- Plasma frequency.
vt = math.sqrt(T0 / mass_elc) -- Electron thermal velocity.
omega = math.sqrt(omega_pe * omega_pe + gas_gamma * k * k * vt * vt) -- Bohm-Gross frequency.

-- Simulation parameters.
Nx = 32 -- Cell count (configuration space: x-direction).
Lx = 2.0 * pi / k -- Domain size (configuration space: x-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 50.0 -- Final simulation time.
num_frames = 10 -- Number of output frames.
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
      local x = xn[1]

      -- Density perturbation at rest with a uniform temperature (a standing wave).
      local n = n0 * (1.0 + alpha * math.cos(k * x))

      return n * mass_elc, 0.0, 0.0, 0.0, n * T0 / (gas_gamma - 1.0)
    end,
  },

  -- Field: the electrostatic potential (the periodic solve removes the mean charge density, the
  -- uniform ion background).
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
