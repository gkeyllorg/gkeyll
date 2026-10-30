-- Constant 4th-order diffusion of a 3D sine wave using a serendipity p1 DG discretization of the
-- advection-diffusion equation.
-- The wave sin(x) sin(y) sin(z) in the periodic 2 pi box decays at the rate D k^4 = 1 per direction
-- (k = 1), so the integral of f^2 decays as exp(-2 x 3 t).
-- Figures of merit (8 x 8 x 8 cells): measured decay rate of the integral of f^2 2.88 against 3;
-- its ratio at t = 0.1 is 0.5621 against 0.5488.

local Vlasov = G0.Vlasov
local LinearAdvection = G0.Vlasov.Eq.LinearAdvection

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
v_advect = 1.0 -- Advection velocity.
diffusion_coeff = 1.0 -- Diffusion coefficient.
diffusion_order = 4 -- Diffusion order.

-- Simulation parameters.
Nx = 8 -- Cell count (configuration space: x-direction).
Ny = 8 -- Cell count (configuration space: y-direction).
Nz = 8 -- Cell count (configuration space: z-direction).
Lx = 2.0 * pi -- Domain size (configuration space: x-direction).
Ly = 2.0 * pi -- Domain size (configuration space: y-direction).
Lz = 2.0 * pi -- Domain size (configuration space: z-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.1 -- Final simulation time.
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
  decompCuts = { 1, 1, 1 }, -- Cuts in each coordinate direction (x-, y- and z-directions only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1, 2, 3 }, -- Periodic directions (x-, y- and z-directions only).
  
  -- Fluid.
  fluid = Vlasov.FluidSpecies.new {
    equation = LinearAdvection.new { },

    -- Constant advection function.
    appAdvect = function (t, xn)
      local ux = 0.0 -- Advection velocity (x-direction).
      local uy = 0.0 -- Advection velocity (y-direction).
      local uz = 0.0 -- Advection velocity (z-direction).

      return ux, uy, uz
    end,
    
    -- Initial conditions function.
    init = function (t, xn)
      local x, y, z = xn[1], xn[2], xn[3]

      local f = math.sin(x) * math.sin(y) * math.sin(z) -- Advected quantity.

      return f
    end,

    -- Diffusion.
    diffusion = {
      diffusionCoefficient = diffusion_coeff,
      diffusionOrder = diffusion_order
    },

    evolve = true -- Evolve species?
  },

  skipField = true
}

-- Run application.
vlasovApp:run()