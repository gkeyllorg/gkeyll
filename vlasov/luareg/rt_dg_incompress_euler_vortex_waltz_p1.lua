-- Vortez-Waltz with incompressible Euler equations. 
-- Input parameters match the initial conditions found in entry JE13 of Ammar's Simulation Journal 
-- (https://ammar-hakim.org/sj/je/je13/je13-incomp-euler-2d.html)

local Vlasov = G0.Vlasov
local IncompressEuler = G0.Vlasov.Eq.IncompressEuler

x1 = 3.5 -- x location of first Gaussian
y1 = 5.0 -- y location of first Gaussian
x2 = 6.5 -- x location of second Gaussian
y2 = 5.0 -- y location of second Gaussian
w = 0.8 -- width of Gaussians

-- Simulation parameters.
Nx = 64 -- Cell count (x-direction).
Ny = 64 -- Cell count (y-direction).
Lx = 10.0 -- Domain size (x-direction).
Ly = 10.0 -- Domain size (y-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 100.0 -- Final simulation time.
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
  lower = { 0.0, 0.0 },
  upper = { Lx, Ly },
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
    equation = IncompressEuler.new { },
    
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local r1 = (x-x1)^2 + (y-y1)^2
      local r2 = (x-x2)^2 + (y-y2)^2
      local chi = math.exp(-r1/w^2) + math.exp(-r2/w^2)
      
      return chi
    end,
  },

  skipField = true
}

-- Run application.
vlasovApp:run()