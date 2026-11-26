-- Collision of two strong shocks with the DG Euler solver (1x), the fluid counterpart of
-- rt_pkpm_neut_bgk_two_shock_1x.c. The states of test 4 of E. F. Toro, Riemann Solvers
-- and Numerical Methods for Fluid Dynamics, with gamma = 5/3: (rho, u, p) = (5.99924, 19.5975,
-- 460.894) left of x = 0.4 and (5.99242, -6.19633, 46.0895) right of it, to t = 0.035. Two shocks
-- move apart with a contact between them; the reference is the exact Riemann solution.
--
-- Figures of merit (serendipity p1, 100 cells, 361 steps): L1 errors of the cell averages of 1.2%
-- (density, of the mean density), 0.056% (velocity, of the largest |u| + c_s) and 0.25% (pressure,
-- of the mean pressure); shocks at x = 0.3570 and 0.8951 (exact 0.3560 and 0.8939), contact at x =
-- 0.7029 (exact 0.6996); between the shocks pressure 1841.5 (exact 1841.8) and velocity 8.559
-- (exact 8.560), with density 12.725 left of the contact (exact 12.743) and 21.908 right of it
-- (exact 21.925).
local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Physical constants (using normalized code units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.

rhol = 5.99924 -- Left fluid mass density.
ul = 19.5975 -- Left fluid velocity.
pl = 460.894 -- Left fluid pressure.

rhor = 5.99242 -- Right fluid mass density.
ur = -6.19633 -- Right fluid velocity.
pr = 46.0895 -- Right fluid pressure.

-- Simulation parameters.
Nx = 100 -- Cell count (x-direction).
x_min = 0.0 -- Lower domain boundary (x-direction).
Lx = 1.0 -- Domain size (x-direction).
x_loc = 0.4 -- Interface of the two states.
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 0.035 -- Final simulation time.
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
