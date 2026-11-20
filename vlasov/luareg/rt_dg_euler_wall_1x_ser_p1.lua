-- Gas flowing between reflecting walls with the DG Euler solver (1x), the fluid counterpart of
-- rt_pkpm_neut_bgk_wall_1x.c. Density 1 and pressure 1/4 with gamma = 5/3 moving at
-- velocity 1 (Mach 1.55) between reflecting walls at x = 0 and 1, to t = 0.3. A shock leaves the
-- wall the gas runs into and a rarefaction the wall it leaves; the reference is the exact solution
-- by the method of images, one Riemann problem per wall between the gas and its mirror image.
--
-- Figures of merit (serendipity p1, 64 cells, 106 steps): L1 errors of the cell averages of 1.25%
-- (density, of the mean density), 0.56% (velocity, of the largest |u| + c_s) and 0.63% (pressure,
-- of the mean pressure); shock at x = 0.8194 (exact 0.8216) with density 2.6820 (exact 2.6817) and
-- pressure 1.84464 (exact 1.84463) behind it. Between the wall the gas leaves and the tail of the
-- rarefaction the density is 0.0975 (exact 0.1131), and 0.066 in the cell at the wall. The mass and
-- the energy change by 3e-15 and 4e-15 of their initial values.
local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Physical constants (using normalized code units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.

rho0 = 1.0 -- Fluid mass density.
u0 = 1.0 -- Fluid velocity.
p0 = 0.25 -- Fluid pressure.

-- Simulation parameters.
Nx = 64 -- Cell count (x-direction).
x_min = 0.0 -- Lower domain boundary (x-direction).
Lx = 1.0 -- Domain size (x-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 0.3 -- Final simulation time.
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
      local rho = rho0 -- Fluid mass density.
      local u = u0 -- Fluid velocity.
      local p = p0 -- Fluid pressure.

      -- Set fluid mass density, momentum density and total energy density.
      return rho, rho * u, 0.0, 0.0, p / (gas_gamma - 1.0) + 0.5 * rho * u * u
    end,

    bcx = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcReflect }
  },

  skipField = true,
}

vlasovApp:run()
