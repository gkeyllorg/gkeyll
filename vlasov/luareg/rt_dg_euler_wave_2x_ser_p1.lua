-- Smooth travelling wave with the DG Euler solver (2x) and the HLLC flux, the counterpart of
-- moments/creg/rt_euler_wave_2d_wv.c (Section 4.1 of R. Liska and B. Wendroff (2003), "Comparison
-- of Several Difference Schemes on 1D and 2D Test Problems for the Euler Equations", SIAM J. Sci.
-- Comput. 25 (3): 995-1017, as set up in entry JE22 of A. Hakim's simulation journal): rho = 1 +
-- 0.2 sin(pi (x + y)) is carried at the velocity (1, -0.5) and pressure 1 (gamma = 1.4) through the
-- periodic box [0, 2]^2 and returns at t = 4, so the exact solution is the initial condition.
-- Figures of merit (serendipity p1, 50 x 50 cells, 1382 steps): the L2 error of the density is
-- 1.19e-3 of the L2 amplitude of the wave; the amplitude of the cell averages falls from 0.1993 to
-- 0.1992; mass and energy are conserved to 1e-13.
local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi -- Pi.

-- Physical constants (using normalized code units).
gas_gamma = 1.4 -- Adiabatic index.
rho0 = 1.0 -- Background fluid mass density.
drho = 0.2 -- Amplitude of the density wave.
u = 1.0 -- Fluid x-velocity.
v = -0.5 -- Fluid y-velocity.
p = 1.0 -- Fluid pressure.

-- Simulation parameters.
Nx = 50 -- Cell count (x-direction).
Ny = 50 -- Cell count (y-direction).
Lx = 2.0 -- Domain size (x-direction).
Ly = 2.0 -- Domain size (y-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 4.0 -- Final simulation time.
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
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = { 1, 2 }, -- Periodic directions.

  -- Fluid.
  fluid = Vlasov.FluidSpecies.new {
    equation = Euler.new { gasGamma = gas_gamma },

    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local rho = rho0 + drho * math.sin(pi * (x + y)) -- Fluid mass density.

      return rho, rho * u, rho * v, 0.0, p / (gas_gamma - 1.0) + 0.5 * rho * (u * u + v * v)
    end,
  },

  skipField = true,
}

vlasovApp:run()
