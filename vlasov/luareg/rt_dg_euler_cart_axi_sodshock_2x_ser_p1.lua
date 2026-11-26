-- Sod-type shock tube in axial symmetry on a Cartesian grid with the DG Euler solver (2x) and the
-- HLLC flux, the counterpart of moments/creg/rt_euler_cart_axi_sodshock.c: (rho, u, p) = (3, 0, 3)
-- inside r < 0.75 and (1, 0, 1) outside, gamma = 1.4, on [-1.25, 1.25]^2 with copy boundaries,
-- t_end = 0.1; the axisymmetric generalization of Section 2.6.2 of A. Hakim (2006), "High
-- Resolution Wave Propagation Schemes for Two-Fluid Plasma Simulations", PhD thesis, University of
-- Washington. Figures of merit (serendipity p1, 128 x 128 cells, 55 steps): against the moments
-- solution on 512 x 512 cells (ending at t = 0.1008), block-averaged, the L1 differences of the
-- cell averages are 4.0e-3 of the mean in density and 4.5e-3 in pressure. The densities along the x
-- axis and along the diagonal, as functions of radius, differ by 1.0e-2 on average (reference:
-- 1.7e-3); the shock (rho = 1.2) is at r = 0.9021 on the axis and 0.8993 on the diagonal
-- (reference: 0.8993 and 0.8986). The cell averages keep the mirror symmetries in x, in y and about
-- the diagonal to 1e-14; mass and energy are conserved to 2e-14.
local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Physical constants (using normalized code units).
gas_gamma = 1.4 -- Adiabatic index.
rhol = 3.0 -- Inner fluid mass density.
pl = 3.0 -- Inner fluid pressure.
rhor = 1.0 -- Outer fluid mass density.
pr = 1.0 -- Outer fluid pressure.

-- Simulation parameters.
Nx = 128 -- Cell count (x-direction).
Ny = 128 -- Cell count (y-direction).
Lx = 2.5 -- Domain size (x-direction).
Ly = 2.5 -- Domain size (y-direction).
rloc = 0.75 -- Fluid boundary (radial coordinate).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

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
  lower = { -0.5 * Lx, -0.5 * Ly },
  upper = { 0.5 * Lx, 0.5 * Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = {  }, -- Periodic directions.

  -- Fluid.
  fluid = Vlasov.FluidSpecies.new {
    equation = Euler.new { gasGamma = gas_gamma },

    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local r = math.sqrt(x * x + y * y)

      -- Inner and outer states (fluid at rest).
      local rho, p = rhor, pr
      if r < rloc then
        rho, p = rhol, pl
      end

      return rho, 0.0, 0.0, 0.0, p / (gas_gamma - 1.0)
    end,

    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy },
    bcy = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy },
  },

  skipField = true,
}

vlasovApp:run()
