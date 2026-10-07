-- Acoustic oscillation of a pressure perturbation with the DG Euler solver (1x). A gas at rest with
-- density 1 and pressure 1 (gamma = 1.4) carries the pressure perturbation 0.01 sin(x) on a
-- periodic domain of length 2 pi. The perturbation splits into two sound waves and, after one
-- acoustic period 2 pi / c_s with c_s = sqrt(1.4), the state must return to its initial one up to
-- the 1e-4 nonlinear correction. Figures of merit (serendipity p2, 128 cells, 715 steps): after one
-- period the pressure differs from its initial profile by at most 4.6e-5 and the density by 3.3e-5,
-- both at the 1e-4 level of the second-order correction to the linear wave, the residual velocity
-- is below 1.2e-4, and mass and energy are conserved to 4e-14.
local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

pi = math.pi -- Pi.
gas_gamma = 1.4 -- Adiabatic index.
rho0 = 1.0 -- Background mass density.
p0 = 1.0 -- Background pressure.
dp = 0.01 -- Pressure perturbation amplitude.
c_s = math.sqrt(gas_gamma * p0 / rho0) -- Sound speed.
t_period = 2.0 * math.pi / c_s -- Acoustic period of the k = 1 mode.

Nx = 128 -- Cell count (x-direction).
Lx = 2.0 * math.pi -- Domain size (x-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 1.0 * t_period -- Final simulation time.
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

  -- Fluid.
  fluid = Vlasov.FluidSpecies.new {
    equation = Euler.new { gasGamma = gas_gamma },

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      -- Standing pressure perturbation on a uniform gas at rest.
      local rho = rho0
      local p = p0 + dp * math.sin(x)

      return rho, 0.0, 0.0, 0.0, p / (gas_gamma - 1.0)
    end,
  },

  skipField = true,
}

vlasovApp:run()
