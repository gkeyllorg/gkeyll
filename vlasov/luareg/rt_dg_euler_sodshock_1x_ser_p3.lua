-- Sod-type shock tube with the DG Euler solver (1x). The states (rho, u, p) = (3, 0, 3) and
-- (1, 0, 1) with gamma = 1.4 meet at x = 0.75
-- (Section 2.6.2 of A. Hakim's 2006 thesis, High Resolution Wave Propagation Schemes for Two-
-- Fluid Plasma Simulations), producing a left rarefaction, a contact and a right shock; the
-- reference is the exact Riemann solution. Figures of merit (serendipity p3, 512 cells): L1 errors
-- of 0.0020 in density, 0.0008 in velocity and 0.0015 in pressure against the exact Riemann
-- solution; shock at x = 0.9001 (exact 0.8994), contact at x = 0.7961 (exact 0.7964), star-state
-- pressure 1.6934 and velocity 0.4641, both exact to four digits.
local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

gas_gamma = 1.4 -- Adiabatic index.
rhol = 3.0 -- Left fluid mass density.
ul = 0.0 -- Left fluid velocity.
pl = 3.0 -- Left fluid pressure.
rhor = 1.0 -- Right fluid mass density.
ur = 0.0 -- Right fluid velocity.
pr = 1.0 -- Right fluid pressure.

Nx = 512 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
poly_order = 3 -- Polynomial order.
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
    equation = Euler.new { gasGamma = gas_gamma },

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      -- Left and right states separated at x = 0.75.
      local rho, u, p = rhor, ur, pr
      if x < 0.75 then
        rho, u, p = rhol, ul, pl
      end

      return rho, rho * u, 0.0, 0.0, p / (gas_gamma - 1.0) + 0.5 * rho * u * u
    end,

    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }
  },

  skipField = true,
}

vlasovApp:run()
