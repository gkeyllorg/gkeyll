-- Turbulence with the Hasegawa-Wakatani system. 
-- Input parameters match the initial conditions found in entry JE17 of Ammar's Simulation Journal 
-- (https://ammar-hakim.org/sj/je/je17/je17-hasegawa-wakatani.html)

local Vlasov = G0.Vlasov
local HasegawaWakatani = G0.Vlasov.Eq.HasegawaWakatani

alpha = 0.1 -- Adiabatic coupling constant
s = 2.0 -- width of the initial Gaussian density

-- Simulation parameters.
Nx = 32 -- Cell count (x-direction).
Ny = 32 -- Cell count (y-direction).
Lx = 40.0 -- Domain size (x-direction).
Ly = 40.0 -- Domain size (y-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 200.0 -- Final simulation time.
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
  lower = { -Lx/2.0, -Ly/2.0 },
  upper = { Lx/2.0, Ly/2.0 },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,
    
  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1, 2 }, -- Periodic directions.
  
  -- Fluid.
  fluid = Vlasov.FluidSpecies.new {
    equation = HasegawaWakatani.new { alpha = alpha, is_modified = false },

    -- Background (linear) density gradient for driving turbulence. 
    n0 = function(t, xn)
      local x, y = xn[1], xn[2]
      return x
    end, 
    
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]
      local r = x^2 + y^2
      local phi = math.exp(-r/s^2) -- initial potential, same as density
      local zeta = 4.0*(r-s^2)*math.exp(-r/s^2)/s^4 -- grad^2 phi 
      
      return zeta, phi
    end,
  },

  skipField = true
}

-- Run application.
vlasovApp:run()