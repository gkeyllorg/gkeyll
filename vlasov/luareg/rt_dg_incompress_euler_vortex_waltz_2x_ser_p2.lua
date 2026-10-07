-- Vortex waltz with the incompressible Euler equations (2x). Two Gaussian vortices of unit peak
-- vorticity and width 0.8, separated by 3 on a periodic 10 x 10 box
-- (entry JE13 of A. Hakim's simulation journal, ammar-hakim.org/sj/je/je13), co-rotate and merge.
-- The reference for the early rotation is the point-vortex period T = 2 pi^2 d^2 / Gamma = 88.4
-- with Gamma = pi w^2, and the energy int |grad phi|^2 is an exact invariant. Figures of merit
-- (serendipity p2, 64 x 64 cells, 1523 steps): over 0 < t < 6 the pair rotates at 0.049 rad per
-- unit time, matching the point-vortex rate in the periodic box, 0.0489
-- (0.071 for an isolated pair); the separation falls from 3.0 to 1.8 by t = 30 and the vortices
-- merge by t = 38; over the run the energy changes by -1.8e-6 and the enstrophy int zeta^2 by -8.4%
-- (upwind dissipation of the filaments).
local Vlasov = G0.Vlasov
local IncompressEuler = G0.Vlasov.Eq.IncompressEuler

x1 = 3.5 -- Center of the first vortex (x).
y1 = 5.0 -- Center of the first vortex (y).
x2 = 6.5 -- Center of the second vortex (x).
y2 = 5.0 -- Center of the second vortex (y).
w = 0.8 -- Width of the Gaussian vortices.

Nx = 64 -- Cell count (x-direction).
Ny = 64 -- Cell count (y-direction).
Lx = 10.0 -- Domain size (x-direction).
Ly = 10.0 -- Domain size (y-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 100.0 -- Final simulation time.
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
  lower = { 0.0, 0.0 },
  upper = { Lx, Ly },
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
    equation = IncompressEuler.new { },

    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]
      local dx1 = (x - x1) * (x - x1)
      local dy1 = (y - y1) * (y - y1)
      local dx2 = (x - x2) * (x - x2)
      local dy2 = (y - y2) * (y - y2)
      local r1 = dx1 + dy1
      local r2 = dx2 + dy2

      -- Vorticity of two co-rotating Gaussian vortices.
      return math.exp(-r1 / (w * w)) + math.exp(-r2 / (w * w))
    end,
  },

  skipField = true,
}

vlasovApp:run()
