-- Constant advection of a 3D sine wave, with a different speed in each direction, using a tensor p2
-- DG discretization of the advection equation.
-- The wave crosses the periodic box once in x, twice in y (backwards) and three times in z,
-- returning to its initial state at t = 2 pi.
-- Figures of merit (16 x 16 x 16 cells): relative L2 error against the initial state 5.7e-04,
-- maximum cell-average error 2.9e-04 of the peak, integral of f^2 retained to 0.9994.

local Vlasov = G0.Vlasov
local LinearAdvection = G0.Vlasov.Eq.LinearAdvection

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
v_advect = 1.0 -- Advection velocity.

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Ny = 16 -- Cell count (configuration space: y-direction).
Nz = 16 -- Cell count (configuration space: z-direction).
Lx = 2.0 * pi -- Domain size (configuration space: x-direction).
Ly = 2.0 * pi -- Domain size (configuration space: y-direction).
Lz = 2.0 * pi -- Domain size (configuration space: z-direction).
poly_order = 2 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
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
  lower = { 0.0, 0.0, 0.0 },
  upper = { Lx, Ly, Lz },
  cells = { Nx, Ny, Nz },
  cflFrac = cfl_frac,
    
  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1, 1 }, -- Cuts in each coordinate direction (x-, y- and z-directions).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1, 2, 3 }, -- Periodic directions (x-, y- and z-directions).
  
  -- Fluid.
  fluid = Vlasov.FluidSpecies.new {
    equation = LinearAdvection.new { },

    -- Constant advection function.
    appAdvect = function (t, xn)
      local ux = v_advect -- Advection velocity (x-direction).
      local uy = -2.0 * v_advect -- Advection velocity (y-direction).
      local uz = 3.0 * v_advect -- Advection velocity (z-direction).

      return ux, uy, uz
    end,
    
    -- Initial conditions function.
    init = function (t, xn)
      local x, y, z = xn[1], xn[2], xn[3]

      local f = math.sin(x) * math.sin(y) * math.sin(z) -- Advected quantity.

      return f
    end,

    evolve = true -- Evolve species?
  },

  skipField = true
}

-- Run application.
vlasovApp:run()