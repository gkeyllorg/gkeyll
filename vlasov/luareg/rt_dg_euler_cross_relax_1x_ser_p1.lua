-- Cross-species relaxation of two DG Euler fluids (1x, no field): the momentum and energy exchange
-- between fluid species.
-- Homogeneous electrons (T = 1, drifting at half the ion thermal speed) and ions (mass 100, T =
-- 0.5, at rest), both with gamma = 5/3, with constant collision frequencies nu_ei = 5 and nu_ie =
-- (m_e / m_i) nu_ei, so that m_e n_e nu_ei = m_i n_i nu_ie. Run for 10 electron-ion collision times
-- in 20 steps, with (nu_ei + nu_ie) dt between 0.29 and 0.54.
-- The pair exchanges momentum and energy through the fluid states, in one backward-Euler step in
-- each SSP-RK3 stage (the exchange of the PKPM species, see rt_pkpm_hybrid_cross_relax_1x.lua).
-- Figures of merit (2 cells, to t = 2): mass is conserved exactly and the total momentum and total
-- energy to 1e-15. Each step multiplies the drift difference by R((nu_ei + nu_ie) dt) to round-off
-- (1e-14) and the temperature difference by R(2 (nu_ei m_e + nu_ie m_i) / (m_e + m_i) dt) to 5e-4
-- (the frictional heating), with R(z) = [1/3 + (2/3) (3/4 + (1/4)/(1 + z))/(1 + z/4)]/(1 + 2 z/3),
-- the factor of the three backward-Euler stage steps of a step with z = nu dt. The temperatures end
-- at 0.91914 (electrons) and 0.58127 (ions).
local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Physical constants (using normalized code units).
mass_elc = 1.0 -- Electron mass.
mass_ion = 100.0 -- Ion mass.

n0 = 1.0 -- Reference number density.
Te = 1.0 -- Electron temperature.
Ti = 0.5 -- Ion temperature.
gas_gamma = 5.0 / 3.0 -- Adiabatic index of the fluid species.
nu_ei = 5.0 -- Electron-ion collision frequency.
nu_ie = (mass_elc / mass_ion) * nu_ei -- Ion-electron collision frequency.

-- Derived physical quantities (using normalized code units).
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.
ue = 0.5 * vti -- Electron drift velocity.
ui = 0.0 -- Ion drift velocity.

-- Simulation parameters.
Nx = 2 -- Cell count (configuration space: x-direction).
Lx = 0.8 -- Domain size (configuration space: x-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 2.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate integrated L2 norm of distribution function.
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
  decompCuts = { 1 }, -- Cuts in each coordinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Electrons (adiabatic fluid, HLLC fluxes by default): they collide with the ions through the
  -- exchange.
  elc = Vlasov.FluidSpecies.new {
    charge = 0.0, mass = mass_elc,
    equation = Euler.new { gasGamma = gas_gamma },

    -- Initial conditions.
    init = function (t, xn)
      local rho = mass_elc * n0 -- Electron mass density.
      local E = 0.5 * rho * ue * ue + n0 * Te / (gas_gamma - 1.0) -- Electron total energy density.

      return rho, rho * ue, 0.0, 0.0, E
    end,

    collisions = {
      numCrossCollisions = 1,
      collideWith = { "ion" },
      collideWithCrossNu = {
        {
          crossNu = function (t, xn)
            return nu_ei
          end
        }
      }
    }
  },

  -- Ions (adiabatic fluid, HLLC fluxes by default): they collide with the electrons through the
  -- exchange.
  ion = Vlasov.FluidSpecies.new {
    charge = 0.0, mass = mass_ion,
    equation = Euler.new { gasGamma = gas_gamma },

    -- Initial conditions.
    init = function (t, xn)
      local rho = mass_ion * n0 -- Ion mass density.
      local E = 0.5 * rho * ui * ui + n0 * Ti / (gas_gamma - 1.0) -- Ion total energy density.

      return rho, rho * ui, 0.0, 0.0, E
    end,

    collisions = {
      numCrossCollisions = 1,
      collideWith = { "elc" },
      collideWithCrossNu = {
        {
          crossNu = function (t, xn)
            return nu_ie
          end
        }
      }
    }
  },

  skipField = true,
}

vlasovApp:run()
