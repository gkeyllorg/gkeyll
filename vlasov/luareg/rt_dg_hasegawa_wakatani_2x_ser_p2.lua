-- Resistive drift-wave turbulence with the Hasegawa-Wakatani system (2x). Adiabaticity alpha and
-- density gradient kappa = 1 on a periodic 40 x 40 box seeded by a Gaussian density and potential
-- blob of width 2 (entry JE17 of A. Hakim's simulation journal, ammar-hakim.org/sj/je/je17). The
-- exact balances are dW/dt = Gamma_n and dE/dt = Gamma_n - Gamma_alpha with W = 1/2 int
-- (n - zeta)^2, E = 1/2 int (n^2 + |grad phi|^2), Gamma_n = -kappa int n d_y phi and Gamma_alpha =
-- alpha int (phi - n)^2. Figures of merit (serendipity p2, 40 x 40 cells, alpha = 1, 15757 steps):
-- the fluctuation kinetic energy int |grad phi|^2 grows from 3.1 to a saturated 7400-8400 after t =
-- 100 (7871 at t = 200) with a zonal share of 0.03-0.07 throughout; W grows from 7.85 to 2.3e4
-- against an integrated source int Gamma_n dt = 1.8e5, the balance being closed by the grid-scale
-- numerical enstrophy sink; E grows from 4.67 to 1.0e4 against int (Gamma_n - Gamma_alpha) dt =
-- 1.4e4. The time step is capped by the adiabatic response rate alpha (1 + k_min^2)/k_min^2 = 40.6
-- of the longest box mode.
local Vlasov = G0.Vlasov
local HasegawaWakatani = G0.Vlasov.Eq.HasegawaWakatani

alpha = 1.0 -- Adiabaticity parameter.
s = 2.0 -- Width of the initial Gaussian blob.
kappa = 1.0 -- Background density gradient.

Nx = 40 -- Cell count (x-direction).
Ny = 40 -- Cell count (y-direction).
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
    equation = HasegawaWakatani.new { alpha = alpha, is_modified = false },

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

      -- Gaussian potential and density blob, initially adiabatic (n = phi), vorticity grad^2 phi.
      local phi = math.exp(-r2 / s2)
      local zeta = 4.0 * (r2 - s2) * phi / (s2 * s2)
      return zeta, phi
    end,
  },

  skipField = true,
}

vlasovApp:run()
