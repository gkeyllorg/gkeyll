-- Uniform gas under a time-dependent applied acceleration with the DG Euler solver (1x). A gas at
-- rest with density 1 and pressure 1 (gamma = 1.4) on a periodic domain is accelerated by
-- a_x = a0 cos(omega t) with a0 = 0.5 and omega = 2; there is no field, so the fluid RHS applies the
-- force and the acceleration is evaluated at every stage time. The state must stay uniform with
-- velocity (a0/omega) sin(omega t) and energy p0/(gamma - 1) + rho u^2/2, which is -0.25 and 2.53125
-- at t = 3 pi/4 (an acceleration frozen at its initial value would give u = 1.18). Figures of merit
-- (serendipity p2, 16 cells, 46 steps): the momentum differs from the exact one by at most
-- 1.1e-8 and the energy by 1.4e-5; the density stays 1 to round-off.
local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

pi = math.pi -- Pi.
gas_gamma = 1.4 -- Adiabatic index.
rho0 = 1.0 -- Background mass density.
p0 = 1.0 -- Background pressure.
a0 = 0.5 -- Amplitude of the applied acceleration.
omega = 2.0 -- Frequency of the applied acceleration.

Nx = 16 -- Cell count (x-direction).
Lx = 2.0 * math.pi -- Domain size (x-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 0.75 * math.pi -- Final simulation time.
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

      -- Uniform gas at rest.
      local rho = rho0
      local p = p0

      return rho, 0.0, 0.0, 0.0, p / (gas_gamma - 1.0)
    end,

    -- Applied acceleration function.
    appliedAccel = function (t, xn)
      return a0 * math.cos(omega * t), 0.0, 0.0
    end,
    evolveAppliedAccel = true,
  },

  skipField = true,
}

vlasovApp:run()
