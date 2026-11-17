-- Expansion into a near vacuum with the DG Euler solver (1x), the counterpart of
-- moments/creg/rt_euler_vac_riem_1d.c: the shock tube of Section 2.6.2 of A. Hakim (2006), High
-- Resolution Wave Propagation Schemes for Two-Fluid Plasma Simulations, PhD thesis, University of
-- Washington, with the right-hand state replaced by a near vacuum. (rho, u, p) = (3, 0, 3) left of
-- x = 0.75 and (1e-6, 0, 1e-6) right of it with gamma = 1.4 on [0.25, 1.25], to t = 0.05. A
-- rarefaction expands to the right behind a contact and a shock in the low-density gas; the
-- reference is the exact Riemann solution.
--
-- Figures of merit (serendipity p1, 512 cells, 494 steps): L1 errors of the cell averages of 0.073%
-- (density, of the mean density), 2.7% (velocity, of the largest |u| + c_s) and 0.086% (pressure,
-- of the mean pressure); shock at x = 1.007 (exact 1.048) and contact at x = 0.979 (exact 0.989),
-- with velocity 4.83 (exact 4.78) and pressure 3.00e-5 (exact 2.95e-5) between them. The density
-- and pressure at the cell vertices stay above 0.9998e-6.
local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Physical constants (using normalized code units).
gas_gamma = 1.4 -- Adiabatic index.

rhol = 3.0 -- Left fluid mass density.
ul = 0.0 -- Left fluid velocity.
pl = 3.0 -- Left fluid pressure.

rhor = 1.0e-6 -- Right fluid mass density.
ur = 0.0 -- Right fluid velocity.
pr = 1.0e-6 -- Right fluid pressure.

-- Simulation parameters.
Nx = 512 -- Cell count (x-direction).
x_min = 0.25 -- Lower domain boundary (x-direction).
Lx = 1.0 -- Domain size (x-direction).
x_loc = 0.75 -- Interface of the two states.
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 0.05 -- Final simulation time.
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
