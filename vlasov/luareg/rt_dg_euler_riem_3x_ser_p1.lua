-- Spherical Riemann problem between two walls with the DG Euler solver (3x) and the HLLC flux, the
-- counterpart of moments/creg/rt_euler_riem_3d.c: a gas at rest with rho = 1 and pressure 5 inside
-- the sphere of radius 0.2 centred at (0, 0, 0.4) and pressure 1 outside, gamma = 1.4, on [0, 1.5]
-- x [0, 1.5] x [0, 1] with copy boundaries in x and y and reflecting walls in z, run to t = 0.7.
-- The problem is adapted from Section 3.2 of J. O. Langseth and R. J. LeVeque (2000), "A Wave
-- Propagation Method for Three-Dimensional Hyperbolic Conservation Laws", Journal of Computational
-- Physics 165 (1): 126-166, as set up in entry JE23 of A. Hakim's simulation journal.
--
-- Figures of merit (serendipity p1, 37 x 37 x 25 cells, 492 steps): against the moments solutions
-- on the same cells (ending at t = 0.7023) and on 74 x 74 x 50 cells, block-averaged, the L1
-- differences of the cell averages are 1.2e-2 and 7.9e-3 of the mean in density and 1.1e-2 and
-- 8.6e-3 in pressure. The density spans [0.358, 1.574] and the largest pressure is 1.901 (moments
-- on the same cells: [0.390, 1.520] and 1.811). The cell averages are symmetric under the exchange
-- of x and y to 7.3e-15. The copy boundaries at x = 0 and y = 0 pass through the centre of the
-- sphere and do not conserve mass: it grows by 6.1% (moments: 6.3%).
local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Physical constants (using normalized code units).
gas_gamma = 1.4 -- Adiabatic index.
rhol = 1.0 -- Inner fluid mass density.
pl = 5.0 -- Inner fluid pressure.
rhor = 1.0 -- Outer fluid mass density.
pr = 1.0 -- Outer fluid pressure.

-- Simulation parameters.
Nx = 37 -- Cell count (x-direction).
Ny = 37 -- Cell count (y-direction).
Nz = 25 -- Cell count (z-direction).
Lx = 1.5 -- Domain size (x-direction).
Ly = 1.5 -- Domain size (y-direction).
Lz = 1.0 -- Domain size (z-direction).
rloc = 0.2 -- Fluid boundary (radial coordinate).
zloc = 0.4 -- Height of the centre of the sphere.
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 0.7 -- Final simulation time.
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
  decompCuts = { 1, 1, 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = {  }, -- Periodic directions.

  -- Fluid.
  fluid = Vlasov.FluidSpecies.new {
    equation = Euler.new { gasGamma = gas_gamma },

    -- Initial conditions function.
    init = function (t, xn)
      local x, y, z = xn[1], xn[2], xn[3]

      local r = math.sqrt(x * x + y * y + (z - zloc) * (z - zloc))

      -- Inner and outer states (fluid at rest).
      local rho, p = rhor, pr
      if r < rloc then
        rho, p = rhol, pl
      end

      return rho, 0.0, 0.0, 0.0, p / (gas_gamma - 1.0)
    end,

    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy },
    bcy = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy },
    bcz = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcReflect },
  },

  skipField = true,
}

vlasovApp:run()
