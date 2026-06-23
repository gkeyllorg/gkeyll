-- 2D Riemann (quadrant) problem, using the Harten-Lax-van Leer (HLL) Riemann solver, for the 5-moment (Euler) equations.
-- Input parameters match the initial conditions in Section 4.3, Case 3, with final time set to t = 0.8 rather than t = 0.3, from the article:
-- R. Liska and B. Wendroff (2003), "Comparison of Several Difference Schemes on 1D and 2D Test Problems for the Euler Equations",
-- SIAM Journal on Scientific Computing, Volume 25 (3): 995-1017.
-- https://epubs.siam.org/doi/10.1137/S1064827502402120

local Moments = G0.Moments
local Euler = G0.Moments.Eq.Euler

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
Nx = 200 -- Cell count (x-direction).
Ny = 200 -- Cell count (y-direction).
Lx = 1.0 -- Domain size (x-direction).
Ly = 1.0 -- Domain size (y-direction).
cfl_frac = 0.95 -- CFL coefficient.

t_end = 0.8 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

loc = 0.8 -- Fluid boundaries (both x and y coordinates).

momentApp = Moments.App.new {
  
  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.0, 0.0 },
  upper = { Lx, Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction (x- and y-directions).
   
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).
  
  -- Fluid.
  fluid = Moments.Species.new {
    equation = Euler.new {
      gasGamma = gas_gamma,
      rpType = G0.EulerRP.HLL
    },
    
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local rho = 0.0
      local u = 0.0
      local v = 0.0
      local p = 0.0

      if y > loc then
        if x < loc then
          rho = rho_ul -- Fluid mass density (upper-left).
          u = u_ul -- Fluid x-velocity (upper-left).
          v = v_ul -- Fluid y-velocity (upper-left).
          p = p_ul -- Fluid pressure (upper-left).
        else
          rho = rho_ur -- Fluid mass density (upper-right).
          u = u_ur -- Fluid x-velocity (upper-right).
          v = v_ur -- Fluid y-velocity (upper-right).
          p = p_ur -- Fluid pressure (upper-right).
        end
      else
        if x < loc then
          rho = rho_ll -- Fluid mass density (lower-left).
          u = u_ll -- Fluid x-velocity (lower-left).
          v = v_ll -- Fluid y-velocity (lower-left).
          p = p_ll -- Fluid pressure (lower-left).
        else
          rho = rho_lr -- Fluid mass density (lower-right).
          u = u_lr -- Fluid x-velocity (lower-right).
          v = v_lr -- Fluid y-velocity (lower-right).
          p = p_lr -- Fluid pressure (lower-right).
        end
      end
  
      local mom_x = rho * u -- Fluid momentum density (x-direction).
      local mom_y = rho * v -- Fluid momentum density (y-direction).
      local mom_z = 0.0 -- Fluid momentum density (z-direction).
      local Etot = (p / (gas_gamma - 1.0)) + (0.5 * rho * (u * u + v * v)) -- Fluid total energy density.
      
      return rho, mom_x, mom_y, mom_z, Etot
    end,
  
    evolve = true -- Evolve species?
  }
}

-- Run application.
momentApp:run()
