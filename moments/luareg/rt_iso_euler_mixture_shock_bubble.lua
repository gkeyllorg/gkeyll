local Moments = G0.Moments
local IsoEulerMixture = G0.Moments.Eq.IsoEulerMixture

-- Physical constants (using normalized code units).
vt1 = 1.067947 -- First species thermal velocity.
vt2 = 2.544589 -- Second species thermal velocity.

rho_pre = 1.0 -- Pre-shock fluid mass density.
u_pre = 0.0 -- Pre-shock fluid velocity (x-direction).
alpha1_pre = 0.99999 -- Pre-shock volume fraction (first species).

rho_post = 1.3764 -- Post-shock fluid mass density.
u_post = -0.3336 -- Post-shock fluid velocity (x-direction).
alpha1_post = 0.99999 -- Post-shock volume fraction (first species).

rho_bub = 0.1818 -- Bubble fluid mass density.
u_bub = 0.0 -- Bubble fluid velocity (x-direction).
alpha1_bub = 0.00001 -- Bubble volume fraction (first species).

-- Simulation parameters.
Nx = 325 -- Cell count (x-direction).
Ny = 89 -- Cell count (y-direction).
Lx = 0.325 -- Domain size (x-direction).
Ly = 0.089 -- Domain size (y-direction).
cfl_frac = 0.95 -- CFL coefficient.

t_end = 0.4 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

x_loc = 0.225 -- Shock location (x-direction).
bub_loc_x = 0.175 -- Bubble location (x-direction).
bub_loc_y = 0.5 * Ly -- Bubble location (y-direction).
bub_rad = 0.025 -- Bubble radius.

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
    equation = IsoEulerMixture.new {
      numComponents = 2,
      vThermal = {
        vt1,
        vt2
      }
    },
    
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local rho1 = 0.0
      local rho2 = 0.0
      local alpha1 = 0.0
    
      local vx_total = 0.0
      local vy_total = 0.0
      local vz_total = 0.0
    
      local r = math.sqrt(((x - bub_loc_x) * (x - bub_loc_x)) + ((y - bub_loc_y) * (y - bub_loc_y)))
    
      if x > x_loc then
        rho1 = rho_post -- First species fluid mass density (post-shock).
        rho2 = rho_bub -- Second species fluid mass density (bubble).
        alpha1 = alpha1_post -- First species volume fraction (post-shock).
    
        vx_total = u_post -- Total mixture velocity (post-shock).
      else
        rho1 = rho_pre -- First species fluid mass density (pre-shock).
        rho2 = rho_bub -- Second species fluid mass density (bubble).
        alpha1 = alpha1_pre -- First species volume fraction (pre-shock).
    
        vx_total = u_pre -- Total mixture velocity (pre-shock).
      end
    
      if r < bub_rad then
        rho1 = rho_pre -- First species fluid mass density (pre-shock).
        rho2 = rho_bub -- Second species fluid mass density (bubble).
        alpha1 = alpha1_bub -- First species volume fraction (bubble).
    
        vx_total = u_bub -- Total mixture velocity (bubble).
      end
    
      local rho_total = (alpha1 * rho1) + ((1.0 - alpha1) * rho2) -- Total mixture density.

      local momx_total = rho_total * vx_total -- Total mixture momentum (x-direction).
      local momy_total = rho_total * vy_total -- Total mixture momentum (y-direction).
      local momz_total = rho_total * vz_total -- Total mixture momentum (z-direction).

      local vol_frac1 = rho_total * alpha1 -- Mixture weighted volume fraction (first species).
      local mass_frac1 = alpha1 * rho1 -- Mixture volume-weighted mass density (first species).
      local mass_frac2 = (1.0 - alpha1) * rho2 -- Mixture volume-weighted mass density (second species).
      
      return rho_total, momx_total, momy_total, momz_total, vol_frac1, mass_frac1, mass_frac2
    end,
  
    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }, -- Copy boundary conditions (x-direction).
    bcy = { G0.SpeciesBc.bcWall, G0.SpeciesBc.bcWall } -- Wall boundary conditions (y-direction).
  }
}

-- Run application.
momentApp:run()
