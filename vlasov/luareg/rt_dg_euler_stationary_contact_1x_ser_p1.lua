-- Strong shock with a stationary contact with the DG Euler solver (1x), the fluid counterpart of
-- rt_pkpm_neut_bgk_stationary_contact_1x.c. Density 1 moving at -18.281207 with gamma =
-- 5/3, with pressure 1000 left of x = 0.8 and 0.01 right of it, to t = 0.012 (the left half of the
-- Woodward-Colella blast wave in the frame in which its contact is at rest). A rarefaction moves to
-- the left, a Mach 189 shock to the right and the contact stays at x = 0.8; the reference is the
-- exact Riemann solution.
--
-- Figures of merit (serendipity p1, 100 cells, 239 steps): L1 errors of the cell averages of 1.5%
-- (density, of the mean density), 0.28% (velocity, of the largest |u| + c_s) and 0.86% (pressure,
-- of the mean pressure); contact at x = 0.8006 (exact 0.8) with its density jump inside two cells,
-- shock at x = 0.8730 (exact 0.8731); density 0.6161 left of the contact (exact 0.6157) and 4.0026
-- right of it (exact 3.9997), pressure 444.2 between the contact and the shock (exact 445.6).
local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Physical constants (using normalized code units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.

rhol = 1.0 -- Left fluid mass density.
ul = -18.281207 -- Left fluid velocity.
pl = 1000.0 -- Left fluid pressure.

rhor = 1.0 -- Right fluid mass density.
ur = -18.281207 -- Right fluid velocity.
pr = 0.01 -- Right fluid pressure.

-- Simulation parameters.
Nx = 100 -- Cell count (x-direction).
x_min = 0.0 -- Lower domain boundary (x-direction).
Lx = 1.0 -- Domain size (x-direction).
x_loc = 0.8 -- Interface of the two states.
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 0.012 -- Final simulation time.
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
