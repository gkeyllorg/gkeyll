-- Double shear layer with the incompressible Euler equations (2x). Two tanh shear layers of width
-- pi/15 at y = pi/2 and y = 3 pi/2 on a periodic 2 pi box
-- (entry JE13 of A. Hakim's simulation journal, ammar-hakim.org/sj/je/je13) are seeded with a k = 1
-- vorticity perturbation of amplitude 1e-3 and roll up. The reference is the Rayleigh-equation
-- growth rate 0.7007 of the k = 1 mode, and the energy int |grad phi|^2 is an exact invariant.
-- Figures of merit (serendipity p2, 64 x 64 cells, 786 steps): the k = 1 mode grows at 0.68-0.70
-- per unit time over 2 < t < 5 against the Rayleigh rate 0.7007 and the layers have rolled up by t
-- = 12; the energy changes by -1.2e-7 and the enstrophy int zeta^2 by -1.6%.
local Vlasov = G0.Vlasov
local IncompressEuler = G0.Vlasov.Eq.IncompressEuler

pi = math.pi -- Pi.
rho = math.pi / 15.0 -- Shear layer width.
delta = 0.001 -- Amplitude of the k = 1 vorticity perturbation.

Nx = 64 -- Cell count (x-direction).
Ny = 64 -- Cell count (y-direction).
Lx = 2.0 * math.pi -- Domain size (x-direction).
Ly = 2.0 * math.pi -- Domain size (y-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 12.0 -- Final simulation time.
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

      -- Vorticity of two tanh shear layers at y = pi/2 and y = 3 pi/2.
      local y_up = 1.5 * pi
      local y_lo = 0.5 * pi
      local chi = 0.0
      if y > pi then
        local arg = (y_up - y) / rho
        chi = 1.0 / (rho * math.cosh(arg) * math.cosh(arg))
      else
        local arg = (y - y_lo) / rho
        chi = -1.0 / (rho * math.cosh(arg) * math.cosh(arg))
      end

      -- Add the k = 1 perturbation.
      local pert = delta * math.cos(x)
      return chi + pert
    end,
  },

  skipField = true,
}

vlasovApp:run()
