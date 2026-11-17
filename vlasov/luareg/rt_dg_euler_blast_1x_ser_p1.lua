-- Interacting blast waves with the DG Euler solver (1x), the fluid counterpart of
-- pkpm/creg/rt_pkpm_neut_bgk_blast_1x.c: the problem of P. Woodward and P. Colella (1984), "The
-- numerical simulation of two-dimensional fluid flow with strong shocks", Journal of Computational
-- Physics 54 (1): 115-173. Density 1 at rest with gamma = 5/3 between reflecting walls at x = 0 and
-- 1; pressure 1000 left of x = 0.1, 0.01 between 0.1 and 0.9, and 100 right of 0.9, to t = 0.038.
-- The reference is a second-order finite-volume solution (MUSCL, HLLE) on 6400 cells.
--
-- Figures of merit (serendipity p1, 160 cells, 757 steps): L1 errors of the cell averages against
-- the reference of 7.3% (density, of the mean density), 0.54% (velocity, of the largest |u| + c_s)
-- and 2.6% (pressure, of the mean pressure); largest density 3.67 at x = 0.778 (reference 4.15 at x
-- = 0.784 on the same cells), largest pressure 311.3 (reference 318.5), largest velocity 13.26
-- (reference 13.17). The mass and the energy change by 3e-14 and 4e-14 of their initial values.
local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Physical constants (using normalized code units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.

rhol = 1.0 -- Left fluid mass density.
ul = 0.0 -- Left fluid velocity.
pl = 1000.0 -- Left fluid pressure.

rhom = 1.0 -- Middle fluid mass density.
um = 0.0 -- Middle fluid velocity.
pm = 0.01 -- Middle fluid pressure.

rhor = 1.0 -- Right fluid mass density.
ur = 0.0 -- Right fluid velocity.
pr = 100.0 -- Right fluid pressure.

-- Simulation parameters.
Nx = 160 -- Cell count (x-direction).
x_min = 0.0 -- Lower domain boundary (x-direction).
Lx = 1.0 -- Domain size (x-direction).
x_lm = 0.1 -- Interface of the left and middle states.
x_mr = 0.9 -- Interface of the middle and right states.
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 0.038 -- Final simulation time.
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

      if x < x_lm then
        rho = rhol -- Fluid mass density (left).
        u = ul -- Fluid velocity (left).
        p = pl -- Fluid pressure (left).
      elseif x < x_mr then
        rho = rhom -- Fluid mass density (middle).
        u = um -- Fluid velocity (middle).
        p = pm -- Fluid pressure (middle).
      else
        rho = rhor -- Fluid mass density (right).
        u = ur -- Fluid velocity (right).
        p = pr -- Fluid pressure (right).
      end

      -- Set fluid mass density, momentum density and total energy density.
      return rho, rho * u, 0.0, 0.0, p / (gas_gamma - 1.0) + 0.5 * rho * u * u
    end,

    bcx = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcReflect }
  },

  skipField = true,
}

vlasovApp:run()
