-- Pressure pulse in a box of reflecting walls with the DG Euler solver (2x) and the HLLC flux, the
-- counterpart of moments/creg/rt_euler_reflect_2d.c: a gas at rest with rho = 1 and p = 1 + 0.1
-- exp(-50 r^2), r measured from the centre of the unit square, gamma = 1.4, reflecting walls on all
-- four sides, run to t = 0.5.
--
-- Figures of merit (serendipity p1, 128 x 128 cells, 523 steps): the integrated mass and total
-- energy are conserved to 5e-14 throughout the run. Against the moments solution on 512 x 512 cells
-- (ending at t = 0.5007), block-averaged, the L2 difference of the pressure is 1.1% of the L2
-- deviation of the reference pressure from its mean; the pressure spans [0.9935, 1.0477] (reference
-- [0.9933, 1.0480]). The density cell averages keep the mirror symmetries in x, in y and about the
-- diagonal to 3e-7.
local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Physical constants (using normalized code units).
gas_gamma = 1.4 -- Adiabatic index.
rho0 = 1.0 -- Reference fluid mass density.
p0 = 1.0 -- Reference fluid pressure.
dp = 0.1 -- Amplitude of the pressure pulse.
beta = 50.0 -- Beta parameter in exponential pressure distribution.

-- Simulation parameters.
Nx = 128 -- Cell count (x-direction).
Ny = 128 -- Cell count (y-direction).
Lx = 1.0 -- Domain size (x-direction).
Ly = 1.0 -- Domain size (y-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 0.5 -- Final simulation time.
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
  periodicDirs = {  }, -- Periodic directions.

  -- Fluid.
  fluid = Vlasov.FluidSpecies.new {
    equation = Euler.new { gasGamma = gas_gamma },

    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local xc = 0.5 * Lx
      local yc = 0.5 * Ly

      local rho = rho0 -- Fluid mass density.
      local p = p0 + dp * math.exp(-beta * ((x - xc) * (x - xc) + (y - yc) * (y - yc))) -- Fluid pressure.

      return rho, 0.0, 0.0, 0.0, p / (gas_gamma - 1.0)
    end,

    bcx = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcReflect },
    bcy = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcReflect },
  },

  skipField = true,
}

vlasovApp:run()
