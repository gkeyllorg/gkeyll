local Moments = G0.Moments
local EulerRGFM = G0.Moments.Eq.EulerRGFM

-- Physical constants (using normalized code units).
gas_gamma1 = 1.4 -- First species adiabatic index.
gas_gamma2 = 1.648 -- Second species adiabatic index.

rho_pre = 1.0 -- Pre-shock fluid mass density.
u_pre = 0.0 -- Pre-shock fluid velocity (x-direction).
phi1_pre = 0.99999 -- Pre-shock level set value (first species).

rho_post = 1.3764 -- Post-shock fluid mass density.
u_post = -0.3336 -- Post-shock fluid velocity (x-direction).
phi1_post = 0.99999 -- Post-shock level set value (first species).

rho_bub = 0.1818 -- Bubble fluid mass density.
u_bub = 0.0 -- Bubble fluid velocity (x-direction).
phi1_bub = 0.00001 -- Bubble level set value (first species).

-- Derived physical quantities (using normalized code units).
p_pre = 1.0 / gas_gamma1 -- Pre-shock fluid pressure.
p_post = 1.5698 / gas_gamma1 -- Post-shock fluid pressure.
p_bub = 1.0 / gas_gamma1 -- Bubble fluid pressure.

-- Simulation parameters.
Nx = 325 -- Cell count (x-direction).
Ny = 89 -- Cell count (y-direction).
Lx = 0.325 -- Domain size (x-direction).
Ly = 0.089 -- Domain size (y-direction).
cfl_frac = 0.95 -- CFL coefficient.
reinit_freq = 3 -- Reinitialization frequency (for level set).

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
    equation = EulerRGFM.new {
      numComponents = 2,
      gasGamma = {
        gas_gamma1,
        gas_gamma2
      },
      reinitFreq = reinit_freq
    },
    
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local rho1 = 0.0
      local rho2 = 0.0
      local phi1 = 0.0
    
      local vx_total = 0.0
      local vy_total = 0.0
      local vz_total = 0.0
      local p_total = 0.0
    
      local r = math.sqrt(((x - bub_loc_x) * (x - bub_loc_x)) + ((y - bub_loc_y) * (y - bub_loc_y)))
    
      if x > x_loc then
        rho1 = rho_post -- First species fluid mass density (post-shock).
        rho2 = rho_bub -- Second species fluid mass density (bubble).
        phi1 = phi1_post -- First species level set value (post-shock).
    
        vx_total = u_post -- Total fluid velocity (post-shock).
        p_total = p_post -- Total fluid pressure (post-shock).
      else
        rho1 = rho_pre -- First species fluid mass density (pre-shock).
        rho2 = rho_bub -- Second species fluid mass density (bubble).
        phi1 = phi1_pre -- First species level set value (pre-shock).
    
        vx_total = u_pre -- Total fluid velocity (pre-shock).
        p_total = p_pre -- Total fluid pressure (pre-shock).
      end
    
      if r < bub_rad then
        rho1 = rho_pre -- First species fluid mass density (pre-shock).
        rho2 = rho_bub -- Second species fluid mass density (bubble).
        phi1 = phi1_bub -- First species level set value (bubble).
    
        vx_total = u_bub -- Total fluid velocity (bubble).
        p_total = p_bub -- Total fluid pressure (bubble).
      end
    
      local rho_total = (phi1 * rho1) + ((1.0 - phi1) * rho2) -- Total fluid density.

      local momx_total = rho_total * vx_total -- Total fluid momentum density (x-direction).
      local momy_total = rho_total * vy_total -- Total fluid momentum density (y-direction).
      local momz_total = rho_total * vz_total -- Total fluid momentum density (z-direction).
    
      local E1 = (p_total / (gas_gamma1 - 1.0)) + (0.5 * rho1 * (vx_total * vx_total)) -- First species total energy.
      local E2 = (p_total / (gas_gamma2 - 1.0)) + (0.5 * rho2 * (vx_total * vx_total)) -- Second species total energy.
      local E_total = (phi1 * E1) + ((1.0 - phi1) * E2) -- Total fluid energy.

      local level_set1 = rho_total * phi1 -- Conserved level set value (first species).
      local mass_frac1 = phi1 * rho1 -- Conserved mass density (first species).
      local mass_frac2 = (1.0 - phi1) * rho2 -- Conserved mass density (second species).
      
      return rho_total, momx_total, momy_total, momz_total, E_total, level_set1, mass_frac1, mass_frac2, 0.0
    end,
  
    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }, -- Copy boundary conditions (x-direction).
    bcy = { G0.SpeciesBc.bcWall, G0.SpeciesBc.bcWall } -- Wall boundary conditions (y-direction).
  }
}

-- Run application.
momentApp:run()
