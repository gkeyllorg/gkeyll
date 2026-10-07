-- Kelvin-Helmholtz instability with the DG Euler solver (2x). A dense stream
-- (density 2, velocity -0.5) between two lighter streams (density 1, velocity 0.5) at pressure 2.5
-- with gamma = 1.4, after Section 4 of M. Sementilli, R. Zangeneh and J. Chen, Fluids 9 (2024) 52,
-- on a periodic unit box; a deterministic multi-mode velocity perturbation of amplitude 0.01 seeds
-- the shear layers and the run follows the instability into its nonlinear roll-up. Figures of merit
-- (serendipity p1, 64 x 64 cells): the y-directed kinetic energy density grows from 5.4e-4 to a
-- peak of 0.117 at t = 5.4 and ends at 0.086 (the incompressible vortex-sheet rate of the m = 1
-- mode, k du sqrt(rho_in rho_out)/(rho_in + rho_out) = 2.96, sets the time scale; the multi-mode
-- seed has no single-mode linear phase to fit); kinetic energy 0.1886 -> 0.1564 is converted into
-- internal energy 6.2500 -> 6.2822 while mass, x-momentum and total energy are conserved to 6e-13.
local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

pi = math.pi -- Pi.
gas_gamma = 1.4 -- Adiabatic index.
rho_in = 2.0 -- Inner fluid mass density.
u_in = -0.5 -- Inner fluid velocity (x-direction).
rho_out = 1.0 -- Outer fluid mass density.
u_out = 0.5 -- Outer fluid velocity (x-direction).
p0 = 2.5 -- Pressure.
yloc = 0.25 -- Half width of the inner fluid (|y| < yloc).
alpha = 0.01 -- Velocity perturbation amplitude.

Nx = 64 -- Cell count (x-direction).
Ny = 64 -- Cell count (y-direction).
Lx = 1.0 -- Domain size (x-direction).
Ly = 1.0 -- Domain size (y-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 10.0 -- Final simulation time.
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
    equation = Euler.new { gasGamma = gas_gamma },

    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      -- Shear layers at y = +/- yloc between a dense inner stream and a lighter outer stream.
      local rho, vx, vy, p = rho_out, u_out, 0.0, p0
      if math.abs(y) < yloc then
        rho, vx = rho_in, u_in
      end

      -- Deterministic multi-mode velocity perturbation (modes 1..4 in x and y, fixed phases).
      local k = 2.0 * pi
      for i = 1, 4 do
        for j = 1, 4 do
          local phase = 2.0 * pi * ((i * j) % 7) / 7.0
          vx = vx + alpha * math.sin(i * k * x + j * k * y + phase)
          vy = vy + alpha * math.sin(i * k * x + j * k * y - phase)
        end
      end

      return rho, rho * vx, rho * vy, 0.0, p / (gas_gamma - 1.0) + 0.5 * rho * (vx * vx + vy * vy)
    end,
  },

  skipField = true,
}

vlasovApp:run()
