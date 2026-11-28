-- Quadrant Riemann problem with the DG Euler solver (2x) and the Lax flux, the counterpart of
-- moments/creg/rt_euler_riem_2d_lax.c: Section 4.3, case 3 of R. Liska and B. Wendroff (2003),
-- "Comparison of Several Difference Schemes on 1D and 2D Test Problems for the Euler Equations",
-- SIAM J. Sci. Comput. 25 (3): 995-1017, with gamma = 1.4, the four states meeting at (0.8, 0.8) on
-- the unit square, copy boundaries and t_end = 0.8 (the paper runs to 0.3). Figures of merit
-- (serendipity p1, 100 x 100 cells, 1211 steps): the largest density is 1.763. The density cell
-- averages keep the mirror symmetry of the problem about the diagonal x = y to 2.1e-13. The mass in
-- the box grows from 0.3187 to 0.8946 through the copy boundaries.
local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Physical constants (using normalized code units).
gas_gamma = 1.4 -- Adiabatic index.

rho_ul = 0.5323 -- Upper-left fluid mass density.
u_ul = 1.206 -- Upper-left fluid x-velocity.
v_ul = 0.0 -- Upper-left fluid y-velocity.
p_ul = 0.3 -- Upper-left fluid pressure.

rho_ur = 1.5 -- Upper-right fluid mass density.
u_ur = 0.0 -- Upper-right fluid x-velocity.
v_ur = 0.0 -- Upper-right fluid y-velocity.
p_ur = 1.5 -- Upper-right fluid pressure.

rho_ll = 0.138 -- Lower-left fluid mass density.
u_ll = 1.206 -- Lower-left fluid x-velocity.
v_ll = 1.206 -- Lower-left fluid y-velocity.
p_ll = 0.029 -- Lower-left fluid pressure.

rho_lr = 0.5323 -- Lower-right fluid mass density.
u_lr = 0.0 -- Lower-right fluid x-velocity.
v_lr = 1.206 -- Lower-right fluid y-velocity.
p_lr = 0.3 -- Lower-right fluid pressure.

-- Simulation parameters.
Nx = 100 -- Cell count (x-direction).
Ny = 100 -- Cell count (y-direction).
Lx = 1.0 -- Domain size (x-direction).
Ly = 1.0 -- Domain size (y-direction).
loc = 0.8 -- Fluid boundaries (both x and y coordinates).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 0.8 -- Final simulation time.
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
    equation = Euler.new { gasGamma = gas_gamma, rpType = G0.EulerRP.Lax },

    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local rho, u, v, p = 0.0, 0.0, 0.0, 0.0

      if y > loc then
        if x < loc then
          rho, u, v, p = rho_ul, u_ul, v_ul, p_ul -- Upper-left state.
        else
          rho, u, v, p = rho_ur, u_ur, v_ur, p_ur -- Upper-right state.
        end
      else
        if x < loc then
          rho, u, v, p = rho_ll, u_ll, v_ll, p_ll -- Lower-left state.
        else
          rho, u, v, p = rho_lr, u_lr, v_lr, p_lr -- Lower-right state.
        end
      end

      return rho, rho * u, rho * v, 0.0, p / (gas_gamma - 1.0) + 0.5 * rho * (u * u + v * v)
    end,

    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy },
    bcy = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy },
  },

  skipField = true,
}

vlasovApp:run()
