-- General diffusion (with constant, positive-definite diffusion tensor) of a 3D plane wave using a
-- serendipity p1 DG discretization of the advection-diffusion equation.
-- The wave decays at the rate sum_ij ki kj Dij = 7.2 for Dxx = 1, Dyy = 0.7, Dzz = 1.2, Dxy = 0.5,
-- Dxz = -0.3, Dyz = 0.2 with (kx, ky, kz) = (1, 2, 1), so the integral of f^2 decays as exp(-2 x
-- 7.2 t).
-- Figures of merit (8 x 8 x 8 cells): measured decay rate of the integral of f^2 7.145 against 7.2;
-- its ratio at t = 0.1 is 0.2396 against 0.2369.

local Vlasov = G0.Vlasov
local LinearAdvection = G0.Vlasov.Eq.LinearAdvection

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
v_advect = 1.0 -- Advection velocity.
diffusion_xx = 1.0 -- Diffusion tensor (xx-component).
diffusion_xy = 0.5 -- Diffusion tensor (xy-component).
diffusion_xz = -0.3 -- Diffusion tensor (xz-component).
diffusion_yy = 0.7 -- Diffusion tensor (yy-component).
diffusion_yz = 0.2 -- Diffusion tensor (yz-component).
diffusion_zz = 1.2 -- Diffusion tensor (zz-component).

kx = 1.0 -- Wavenumber (x-direction).
ky = 2.0 -- Wavenumber (y-direction).
kz = 1.0 -- Wavenumber (z-direction).

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
  lower = { -0.5 * Lx, -0.5 * Ly, -0.5 * Lz },
  upper = { 0.5 * Lx, 0.5 * Ly, 0.5 * Lz },
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

      local f = math.sin(kx * x + ky * y + kz * z) -- Advected quantity.

      return f
    end,

    -- Diffusion.
    diffusion = {
      diffusionTensor = function (t, xn)
        return diffusion_xx, diffusion_xy, diffusion_xz, diffusion_yy, diffusion_yz, diffusion_zz
      end
    },

    evolve = true -- Evolve species?
  },

  skipField = true
}

-- Run application.
vlasovApp:run()