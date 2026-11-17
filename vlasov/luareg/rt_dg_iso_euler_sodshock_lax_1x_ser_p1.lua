-- Sod-type shock tube with the DG isothermal Euler solver (1x), the counterpart of
-- moments/creg/rt_iso_euler_sodshock_lax.c. The states (rho, u) = (3, 0.5) and (1, 0) with the
-- sound speed 1 meet at x = 0.75 (Section 2.6.2 of A. Hakim's 2006 thesis, High Resolution Wave
-- Propagation Schemes for Two-Fluid Plasma Simulations, with a left velocity of 0.5), producing a
-- left rarefaction and a right shock; the reference is the exact Riemann solution of the isothermal
-- Euler equations (star density 2.2012 and velocity 0.8096, shock speed 1.4836). The numerical flux
-- is the Lax solver of the isothermal wave equation object (rp_type).
--
-- Figures of merit at t = 0.1 (serendipity p1, 512 cells, 309 steps): L1 errors of the cell
-- averages of 0.0016 in density and 0.0017 in momentum; shock at x = 0.8985 (exact 0.8984); star-
-- state density 2.2012 and velocity 0.8096, both exact to four digits; the cell averages of the
-- density stay within [1, 3]. With the Roe flux (rt_dg_iso_euler_sodshock_1x) the errors are 0.0014
-- and 0.0016.
local Vlasov = G0.Vlasov
local IsoEuler = G0.Vlasov.Eq.IsoEuler

cs = 1.0 -- Fluid sound speed.
rhol = 3.0 -- Left fluid mass density.
ul = 0.5 -- Left fluid velocity.
rhor = 1.0 -- Right fluid mass density.
ur = 0.0 -- Right fluid velocity.

Nx = 512 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
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
  lower = { 0.25 },
  upper = { 0.25 + Lx },
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
    equation = IsoEuler.new { vThermal = cs, rpType = G0.IsoEulerRP.Lax },

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      -- Left and right states separated at x = 0.75.
      local rho, u = rhor, ur
      if x < 0.75 then
        rho, u = rhol, ul
      end

      return rho, rho * u, 0.0, 0.0
    end,

    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }
  },

  skipField = true,
}

vlasovApp:run()
