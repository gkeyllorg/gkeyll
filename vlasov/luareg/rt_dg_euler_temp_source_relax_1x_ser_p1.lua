-- Relaxation of the temperature of an Euler fluid under the temperature source, the moments of the
-- BGK-form source nu_Q (f_M(n, u, T_Q) - f), with the DG Euler solver (1x), explicit update.
-- Two homogeneous drifting gases (gamma = 1.4) start at T = 1 on a periodic domain. The source of
-- the first has T_Q = 10 (heating) and that of the second T_Q = 0.1 (cooling). Source rate 1, run
-- for 5 relaxation times.
-- Figures of merit (serendipity p1, 2 cells, to t = 5): density and momentum are conserved to
-- round-off. The temperature follows T_Q + (T_0 - T_Q) exp(-nu_Q t): 9.9394 and 0.10606 at
-- t = 5, with T - T_Q decaying at rate 1.0000 (largest departure from the exponential 3.4e-5 of
-- T_0 - T_Q, with nu_Q dt up to 0.14).
local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

pi = math.pi -- Pi.
gas_gamma = 1.4 -- Adiabatic index.
mass = 1.0 -- Species mass.
n0 = 1.0 -- Number density.
T0 = 1.0 -- Initial temperature.
ux0 = 0.5 -- Flow velocity (x-direction).
uy0 = 0.3 -- Flow velocity (y-direction).
nu_Q = 1.0 -- Rate of the temperature source.
T_heat = 10.0 -- Temperature of the source of the heated species.
T_cool = 0.1 -- Temperature of the source of the cooled species.

Nx = 2 -- Cell count (x-direction).
Lx = 2.5 -- Domain size (x-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 5.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Initial conditions function of both species.
local function init_fluid(t, xn)
  local rho = mass * n0 -- Mass density.
  local p = n0 * T0 -- Pressure.

  return rho, rho * ux0, rho * uy0, 0.0, 0.5 * rho * (ux0 * ux0 + uy0 * uy0) + p / (gas_gamma - 1.0)
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

  -- Heated species.
  heat = Vlasov.FluidSpecies.new {
    mass = mass,
    equation = Euler.new { gasGamma = gas_gamma },

    -- Initial conditions function.
    init = init_fluid,

    -- Temperature source.
    sourceBGK = {
      rateProfile = function (t, xn)
        return nu_Q
      end,
      tempShape = function (t, xn)
        return T_heat
      end,
    },
  },

  -- Cooled species.
  cool = Vlasov.FluidSpecies.new {
    mass = mass,
    equation = Euler.new { gasGamma = gas_gamma },

    -- Initial conditions function.
    init = init_fluid,

    -- Temperature source.
    sourceBGK = {
      rateProfile = function (t, xn)
        return nu_Q
      end,
      tempShape = function (t, xn)
        return T_cool
      end,
    },
  },

  skipField = true,
}

vlasovApp:run()
