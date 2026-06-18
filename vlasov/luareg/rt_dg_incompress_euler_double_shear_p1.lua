-- Double shear flow with incompressible Euler equations. 
-- Input parameters match the initial conditions found in entry JE13 of Ammar's Simulation Journal 
-- (https://ammar-hakim.org/sj/je/je13/je13-incomp-euler-2d.html)

local Vlasov = G0.Vlasov
local IncompressEuler = G0.Vlasov.Eq.IncompressEuler

rho = math.pi/15.0 -- Shear layer width. 
delta = 0.05 -- Perturbation. 

-- Simulation parameters.
Nx = 64 -- Cell count (x-direction).
Ny = 64 -- Cell count (y-direction).
Lx = 2.0*math.pi -- Domain size (x-direction).
Ly = 2.0*math.pi -- Domain size (y-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 8.0 -- Final simulation time.
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

      local chi = 0.0
      if y > math.pi then
        chi = delta*math.cos(x) + 1.0/rho*(1.0 / math.cosh((3.0*math.pi/2.0 - y) / rho))^2 -- Right shear layer
      else
        chi = delta*math.cos(x) - 1.0/rho*(1.0 / math.cosh((y - math.pi/2.0) / rho))^2 -- Left shear layer
      end
      
      return chi
    end,
  },

  skipField = true
}

-- Run application.
vlasovApp:run()