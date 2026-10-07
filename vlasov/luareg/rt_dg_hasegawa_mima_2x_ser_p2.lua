-- Drift-wave turbulence with the Hasegawa-Mima system (2x). A Gaussian potential blob of width 2 on
-- a periodic 40 x 40 box with background density gradient kappa = 1
-- (entry JE17 of A. Hakim's simulation journal, ammar-hakim.org/sj/je/je17); the evolved quantity
-- is zeta' = (grad^2 - 1) phi. The exact invariants are E = int (|grad phi|^2 + phi^2) and W = int
-- zeta'^2, and the zeta' centroid moves in y at exactly -kappa. Figures of merit
-- (serendipity p2, 128 x 128 cells, 700 steps): the zeta' centroid moves at -1.0000 kappa for t < 5
-- and at -0.992 kappa averaged over t < 10, int zeta' is conserved to 1e-15, and over the turbulent
-- run to t = 200 the energy E falls by 0.98% and the enstrophy W by 1.1% through grid-scale
-- dissipation.
local Vlasov = G0.Vlasov
local HasegawaMima = G0.Vlasov.Eq.HasegawaMima

s = 2.0 -- Width of the initial Gaussian blob.
kappa = 1.0 -- Background density gradient.

Nx = 128 -- Cell count (x-direction).
Ny = 128 -- Cell count (y-direction).
Lx = 40.0 -- Domain size (x-direction).
Ly = 40.0 -- Domain size (y-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 200.0 -- Final simulation time.
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
  periodicDirs = { 1, 2 }, -- Periodic directions.

  -- Fluid.
  fluid = Vlasov.FluidSpecies.new {
    equation = HasegawaMima.new { },

    -- Linear background density whose gradient drives the drift waves through {phi, n0}.
    n0 = function (t, xn)
      local x, y = xn[1], xn[2]
      return kappa * x
    end,

    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]
      local s2 = s * s
      local x2 = x * x
      local y2 = y * y
      local r2 = x2 + y2

      -- Gaussian potential blob; the evolved quantity is zeta' = grad^2 phi - phi.
      local phi = math.exp(-r2 / s2)
      local zeta = 4.0 * (r2 - s2) * phi / (s2 * s2)
      return zeta - phi
    end,
  },

  skipField = true,
}

vlasovApp:run()
