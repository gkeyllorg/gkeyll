local Moments = G0.Moments
local IsoEulerMixture = G0.Moments.Eq.IsoEulerMixture

-- Physical constants (using normalized code units).
vt1 = 396.867658 -- First species thermal velocity.
vt2 = 1100.464666 -- Second species thermal velocity.

rhol = 1.3333 -- Left fluid mass density.
ul = 0.3535 * math.sqrt(math.pow(10.0, 5.0)) -- Left fluid velocity.
alpha1_l = 0.99999 -- Left fluid volume fraction (first species).

rhoc = 1.0 -- Central fluid mass density.
uc = 0.0 -- Central fluid velocity.
alpha1_c = 0.99999 -- Central fluid volume fraction (first species).

rhor = 0.1379 -- Right fluid mass density.
ur = 0.0 -- Right fluid velocity.
alpha1_r = 0.00001 -- Right fluid volume fraction (first species).

-- Simulation parameters.
Nx = 2048 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
cfl_frac = 0.95 -- CFL coefficient.

t_end = 0.0012 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

momentApp = Moments.App.new {
  
  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.0 },
  upper = { Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).
    
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
      local x = xn[1]

      local rho1 = 0.0
      local rho2 = 0.0
      local alpha1 = 0.0
    
      local vx_total = 0.0
      local vy_total = 0.0
      local vz_total = 0.0
    
      if x < 0.05 then
        rho1 = rhol -- First species fluid mass density (left).
        rho2 = rhor -- Second species fluid mass density (right).
        alpha1 = alpha1_l -- First species volume fraction (left).
    
        vx_total = ul -- Total mixture velocity (left).
      elseif x < 0.5 then
        rho1 = rhoc -- First species fluid mass density (central).
        rho2 = rhor -- Second species fluid mass density (right).
        alpha1 = alpha1_c -- First species volume fraction (central).
    
      else
        rho1 = rhoc -- First species fluid mass density (central).
        rho2 = rhor -- Second species fluid mass density (right).
        alpha1 = alpha1_r -- First species volume fraction (right).
    
        vx_total = ur -- Total mixture velocity (right).
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
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (x-direction).
  }
}

-- Run application.
momentApp:run()
