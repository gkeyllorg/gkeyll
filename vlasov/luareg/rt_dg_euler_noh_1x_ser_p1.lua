-- Noh problem with the DG Euler solver (1x), the counterpart of moments/creg/rt_euler_noh_1d.c and
-- rt_pkpm_neut_bgk_noh_1x.c: Section 4 of W. F. Noh (1987), "Errors for calculations of
-- strong shocks using an artificial viscosity and an artificial heat flux", Journal of
-- Computational Physics 72 (1): 78-120, as in problem 3 of entry JE2 of Ammar's Simulation Journal
-- (https://ammar-hakim.org/sj/je/je2/je2-euler-shock.html#problem-3). Density 1 and pressure 1e-6
-- with gamma = 5/3, moving at +1 left of x = 0.5 and at -1 right of it, to t = 1: two shocks leave
-- the centre at speed 1/3 and the gas between them is at rest with density 4 and pressure 4/3. The
-- reference is the exact Riemann solution.
--
-- Figures of merit (serendipity p1, 100 cells, 371 steps): L1 errors of the cell averages of 0.45%
-- (density, of the mean density), 0.051% (velocity, of the largest |u| + c_s) and 0.14% (pressure,
-- of the mean pressure); shocks at x = 0.1662 and 0.8338 (exact 0.1667 and 0.8333). More than three
-- cells from the shocks and the centre the density is within 0.002 of 4 and the pressure within
-- 7e-4 of 4/3; the two central cells have density 3.732 at pressure 1.3334 (Noh's wall heating).
local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Physical constants (using normalized code units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.

rhol = 1.0 -- Left fluid mass density.
ul = 1.0 -- Left fluid velocity.
pl = 1.0e-6 -- Left fluid pressure.

rhor = 1.0 -- Right fluid mass density.
ur = -1.0 -- Right fluid velocity.
pr = 1.0e-6 -- Right fluid pressure.

-- Simulation parameters.
Nx = 100 -- Cell count (x-direction).
x_min = 0.0 -- Lower domain boundary (x-direction).
Lx = 1.0 -- Domain size (x-direction).
x_loc = 0.5 -- Interface of the two states.
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 1.0 -- Final simulation time.
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
  lower = { x_min },
  upper = { x_min + Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = {  }, -- Periodic directions.

  -- Fluid.
  fluid = Vlasov.FluidSpecies.new {
    equation = Euler.new { gasGamma = gas_gamma },

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local rho = 0.0
      local u = 0.0
      local p = 0.0

      if x < x_loc then
        rho = rhol -- Fluid mass density (left).
        u = ul -- Fluid velocity (left).
        p = pl -- Fluid pressure (left).
      else
        rho = rhor -- Fluid mass density (right).
        u = ur -- Fluid velocity (right).
        p = pr -- Fluid pressure (right).
      end

      -- Set fluid mass density, momentum density and total energy density.
      return rho, rho * u, 0.0, 0.0, p / (gas_gamma - 1.0) + 0.5 * rho * u * u
    end,

    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }
  },

  skipField = true,
}

vlasovApp:run()
