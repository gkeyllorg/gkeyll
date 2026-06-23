local Moments = G0.Moments
local EulerMixture = G0.Moments.Eq.EulerMixture

-- Physical constants (using normalized code units).
gas_gamma1 = 1.4 -- First species adiabatic index.
gas_gamma2 = 1.67 -- Second species adiabatic index.

rhol = 1.3333 -- Left fluid mass density.
ul = 0.3535 * math.sqrt(math.pow(10.0, 5.0)) -- Left fluid velocity.
pl = 1.5 * math.pow(10.0, 5.0) -- Left fluid pressure.
alpha1_l = 0.99999 -- Left fluid volume fraction (first species).

rhoc = 1.0 -- Central fluid mass density.
uc = 0.0 -- Central fluid velocity.
pc = 1.0 * math.pow(10.0, 5.0) -- Central fluid pressure.
alpha1_c = 0.99999 -- Central fluid volume fraction (first species).

rhor = 0.1379 -- Right fluid mass density.
ur = 0.0 -- Right fluid velocity.
pr = 1.0 * math.pow(10.0, 5.0) -- Right fluid pressure.
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
    equation = EulerMixture.new {
      numComponents = 2,
      gasGamma = {
        gas_gamma1,
        gas_gamma2
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
      local p_total = 0.0
    
      if x < 0.05 then
        rho1 = rhol -- First species fluid mass density (left).
        rho2 = rhor -- Second species fluid mass density (right).
        alpha1 = alpha1_l -- First species volume fraction (left).
    
        vx_total = ul -- Total mixture velocity (left).
        p_total = pl -- Total mixture pressure (left).
      elseif x < 0.5 then
        rho1 = rhoc -- First species fluid mass density (central).
        rho2 = rhor -- Second species fluid mass density (right).
        alpha1 = alpha1_c -- First species volume fraction (central).
    
        vx_total = uc -- Total mixture velocity (central).
        p_total = pc -- Total mixture pressure (central).
      else
        rho1 = rhoc -- First species fluid mass density (central).
        rho2 = rhor -- Second species fluid mass density (right).
        alpha1 = alpha1_r -- First species volume fraction (right).
    
        vx_total = ur -- Total mixture velocity (right).
        p_total = pr -- Total mixture pressure (right).
      end

      local rho_total = (alpha1 * rho1) + ((1.0 - alpha1) * rho2) -- Total mixture density.

      local momx_total = rho_total * vx_total -- Total mixture momentum density (x-direction).
      local momy_total = rho_total * vy_total -- Total mixture momentum density (y-direction).
      local momz_total = rho_total * vz_total -- Total mixture momentum density (z-direction).
    
      local E1 = (p_total / (gas_gamma1 - 1.0)) + (0.5 * rho1 * (vx_total * vx_total)) -- First species total energy.
      local E2 = (p_total / (gas_gamma2 - 1.0)) + (0.5 * rho2 * (vx_total * vx_total)) -- Second species total energy.
      local E_total = (alpha1 * E1) + ((1.0 - alpha1) * E2) -- Total mixture energy.

      local vol_frac1 = rho_total * alpha1 -- Mixture weighted volume fraction (first species).
      local mass_frac1 = alpha1 * rho1 -- Mixture volume-weighted mass density (first species).
      local mass_frac2 = (1.0 - alpha1) * rho2 -- Mixture volume-weighted mass density (second species).
      
      return rho_total, momx_total, momy_total, momz_total, E_total, vol_frac1, mass_frac1, mass_frac2
    end,
  
    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (x-direction).
  }
}

-- Run application.
momentApp:run()
