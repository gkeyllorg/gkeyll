-- Jeans instability of a self-gravitating fluid with the DG Euler equations and the Poisson field
-- of the Vlasov app (1x), the fluid counterpart of rt_vp_jeans_1x1v. A neutral adiabatic gas
-- (gamma = 3, cs = 1, n0 = 1, m = 1) with gravitational coupling alpha_g = 4 pi G = 1, so that
-- omega_J = 1 and k_J = omega_J/cs = 1, starts in the growing mode at k/k_J = 0.5 with density
-- amplitude 1e-4:
-- adiabatic pressure and the flow that continuity gives the mode. Only nabla^2 phi_g = alpha_g m n
-- is solved (the periodic solve removes the mean density, the Jeans swindle). Linear theory:
-- omega^2 = cs^2 k^2 - omega_J^2, a growth rate gamma = 0.8660 omega_J, so the gravitational field
-- energy grows at 2 gamma = 1.7321 omega_J. The fluid takes -grad(phi_g) in the kicks of the
-- kick-drift-kick split between the SSP-RK3 transport steps.
--
-- Figures of merit (32 cells, serendipity p1, t = 10, 87 steps): a fit of ln int |grad phi_g|^2
-- over t = 1-8 gives 2 gamma = 1.7268 (-0.30%). The total energy (fluid energy minus half the
-- field-grav-energy file, int |grad phi_g|^2 / alpha_g) changes by 8.2e-4 of the final
-- gravitational energy; mass is conserved to 5e-15.

local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi -- Pi.

-- Physical constants (using normalized code units).
gas_gamma = 3.0 -- Adiabatic index.
alpha_g = 1.0 -- Gravitational coupling 4 pi G.
mass = 1.0 -- Particle mass.
charge = 0.0 -- Particle charge (neutral).

n0 = 1.0 -- Reference number density.
T0 = 1.0 / 3.0 -- Temperature.
alpha = 1.0e-4 -- Applied perturbation amplitude.

-- Derived physical quantities (using normalized code units).
omega_J = math.sqrt(alpha_g * mass * n0) -- Jeans frequency.
cs = math.sqrt(gas_gamma * T0 / mass) -- Sound speed.
k_J = omega_J / cs -- Jeans wave number.
k = 0.5 * k_J -- Perturbed wave number.
gamma_J = math.sqrt(omega_J * omega_J - cs * cs * k * k) -- Growth rate.

-- Simulation parameters.
Nx = 32 -- Cell count (configuration space: x-direction).
Lx = 2.0 * pi / k -- Domain size (configuration space: x-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 10.0 -- Final simulation time.
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

  neut = Vlasov.FluidSpecies.new {
    charge = charge, mass = mass,
    equation = Euler.new { gasGamma = gas_gamma },

    init = function (t, xn)
      local x = xn[1]

      -- The growing mode of linear theory: density and adiabatic pressure perturbations with the
      -- flow that continuity gives them.
      local n = n0 * (1.0 + alpha * math.cos(k * x))
      local p = n0 * T0 * (1.0 + gas_gamma * alpha * math.cos(k * x))
      local ux = -(gamma_J * alpha / k) * math.sin(k * x)

      return n * mass, n * mass * ux, 0.0, 0.0, p / (gas_gamma - 1.0) + 0.5 * n * mass * ux * ux
    end,
  },

  -- Field: self-gravity only (the species is neutral, so no electrostatic solve; the periodic
  -- solve removes the mean density, the Jeans swindle).
  field = Vlasov.Field.new {
    fieldID = G0.FieldModel.Phi,
    alphaG = alpha_g,

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
