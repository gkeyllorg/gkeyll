-- Constant advection in 2x using a p1 DG discretization of the advection equation. 

local Vlasov = G0.Vlasov
local LinearAdvection = G0.Vlasov.Eq.LinearAdvection

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
v_advect = 1.0 -- Advection velocity.

r0 = 0.2 -- Distribution radius.
x0 = 1.0 / 4.0 -- Distribution center (x-coordinate).
y0 = 1.0 / 2.0 -- Distribution center (y-coordinate).

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Ny = 16 -- Cell count (configuration space: y-direction).
Lx = 1.0 -- Domain size (configuration space: x-direction).
Ly = 1.0 -- Domain size (configuration space: y-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 2.0 * pi -- Final simulation time.
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
  decompCuts = { 1, 1 }, -- Cuts in each coordinate direction (x- and y-directions only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1, 2 }, -- Periodic directions (x- and y-directions only).
  
  -- Fluid.
  fluid = Vlasov.FluidSpecies.new {
    equation = LinearAdvection.new { },

    -- Constant advection function.
    appAdvect = function (t, xn)
      local x, y = xn[1], xn[2]

      local ux = -y + 0.5 -- Advection velocity (x-direction).
      local uy = x - 0.5 -- Advection velocity (y-direction).
      local uz = 0.0 -- Advection velocity (z-direction).
      
      return ux, uy, uz
    end,
    
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local r = math.min(math.sqrt(((x - x0) * (x - x0)) + ((y - y0) * (y - y0))), r0) / r0
      local f = 0.25 * (1.0 + math.cos(pi * r)) -- Advected quantity.
      
      return f
    end,

    evolve = true -- Evolve species?
  },

  skipField = true
}

-- Run application.
vlasovApp:run()