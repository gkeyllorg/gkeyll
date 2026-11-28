-- Rayleigh-Taylor instability with the DG Euler solver (2x) and the HLLC flux, the counterpart of
-- moments/creg/rt_euler_rt.c: heavy fluid (rho = 2) over light fluid (rho = 1) in [0, 1/6] x
-- [0, 1], interface y = 0.5 + 0.01 cos(2 pi x / Lx), gamma = 1.4, applied acceleration
-- (0, -0.1, 0), reflecting walls on all sides, t_end = 8.5. The fluid starts at rest with the
-- hydrostatic pressure p = 0.2 (1 - y) above the interface, continuous across it, and the total
-- energy (p + 0.01) / (gamma - 1).
--
-- Figures of merit (serendipity p1, 50 x 200 cells, 6883 steps; in brackets the moments solution on
-- 50 x 200 cells): mass is conserved to 4e-13. The kinetic energy grows as exp(2 gamma t) with
-- gamma = 0.920 while it rises from 1e-6 to 1e-5 (0.949) and is 6.65e-4 at t = 8.5 (4.15e-4), when
-- the heavy fluid (rho > 1.5) reaches down to y = 0.0575 and the light fluid up to y = 0.7975
-- (0.1625 and 0.7125). The total energy rises by 1.69e-4 while the potential energy g int rho y dV
-- falls by 1.69e-4. The departure from the mirror symmetry about x = Lx / 2 is 3e-8 at the end; it
-- measures round-off, which the run amplifies by about 1e8 (the C twin, whose initial state differs
-- in the last bits, ends 5e-8 from this run and 2e-8 from the symmetry).
local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi -- Pi.

-- Physical constants (using normalized code units).
gas_gamma = 1.4 -- Adiabatic index.
grav = 0.1 -- Gravitational acceleration.
rho_top = 2.0 -- Top fluid mass density.
rho_bot = 1.0 -- Bottom fluid mass density.
p_ref = 0.01 -- Reference fluid pressure.
pert_max = 0.01 -- Maximum amplitude of initial perturbation.

-- Simulation parameters.
Nx = 50 -- Cell count (x-direction).
Ny = 200 -- Cell count (y-direction).
Lx = 1.0 / 6.0 -- Domain size (x-direction).
Ly = 1.0 -- Domain size (y-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 8.5 -- Final simulation time.
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

      -- Height of the interface between the two fluids.
      local yloc = 0.5 + (pert_max * math.cos(2.0 * pi * (x / Lx)))

      local rho, p = 0.0, 0.0

      if y > yloc then
        rho = rho_top -- Fluid mass density (top).
        p = rho_top * grav * (1.0 - y) -- Fluid pressure (top).
      else
        rho = rho_bot -- Fluid mass density (bottom).
        p = (rho_top * grav * (1.0 - yloc)) + (rho_bot * grav * (yloc - y)) -- Fluid pressure (bottom).
      end

      return rho, 0.0, 0.0, 0.0, (p + p_ref) / (gas_gamma - 1.0)
    end,

    -- Applied acceleration function.
    appliedAccel = function (t, xn)
      return 0.0, -grav, 0.0
    end,

    bcx = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcReflect },
    bcy = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcReflect },
  },

  skipField = true,
}

vlasovApp:run()
