-- Shock tube test on a plane-symmetric Gowdy spacetime for the coupled fluid-Einstein equations, assuming an ultra-relativistic equation of state.
-- Input parameters taken from the initial conditions in Section 9 (Riemann problem 1), from the article:
-- A. P. Barnes, P. G. Lefloch, B. G. Schmidt and J. M. Stewart (2004), "The Glimm scheme for perfect fluids on plane-symmetric Gowdy spacetimes",
-- Classical and Quantum Gravity, Volume 21 (22): 5043.
-- https://iopscience.iop.org/article/10.1088/0264-9381/21/22/003

local Moments = G0.Moments
local GRMedium = G0.Moments.Eq.GRMedium

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
gas_gamma = 4.0 / 3.0 -- Adiabatic index.
kappa = 8.0 * pi -- Stress-energy prefactor in the Einstein field equations.

exp_2a = math.exp(-4.0) -- Exponential appearing in dt and dx metric terms.

-- Derived physical quantities (using normalized code units).
rhol = 100.0 / kappa -- Left fluid mass density.
rhor = 1.0 / kappa -- Right fluid mass density.

Etot_l = rhol -- Left fluid total energy density.
Etot_r = rhor -- Right fluid total energy density.

-- Simulation parameters.
Nx = 4096 -- Cell count (x-direction).
Lx = 2.0 -- Domain size (x-direction).
cfl_frac = 0.95 -- CFL coefficient.

t_end = 0.5 -- Final simulation time.
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
  lower = { -0.5 * Lx },
  upper = { 0.5 * Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).
  
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Fluid.
  fluid = Moments.Species.new {
    equation = GRMedium.new {
      gasGamma = gas_gamma,
      kappa = kappa
    },

    hasEinsteinMedium = true,
    mediumGasGamma = gas_gamma,
    mediumKappa = kappa,
  
    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]
      
      local Etot = 0.0

      if x < 0.0 then
        Etot = Etot_l -- Fluid total energy density (left).
      else
        Etot = Etot_r -- Fluid total energy density (right).
      end
    
      local a_dt = 0.0 -- Time derivative of metric term a.
      local a_dx = 0.0 -- Space derivative of metric term a.
      local b_dt = 0.0 -- Time derivative of metric term b.
      local b_dx = -math.sqrt((kappa * exp_2a * Etot) / 3.0) * math.tan((0.5 * x * math.sqrt(3.0 * kappa * exp_2a * Etot))) -- Space derivative of metric term b.
      local c_dt = 0.0 -- Time derivative of metric term c.
      local c_dx = 0.0 -- Space derivative of metric term c.
    
      local b_dx_plus = -math.sqrt((kappa * exp_2a * Etot) / 3.0) * math.tan((0.5 * (x + (0.5 * math.pow(10.0, -8.0))) * math.sqrt(3.0 * kappa * exp_2a * Etot)))
      local b_dx_minus = -math.sqrt((kappa * exp_2a * Etot) / 3.0) * math.tan((0.5 * (x - (0.5 * math.pow(10.0, -8.0))) * math.sqrt(3.0 * kappa * exp_2a * Etot)))
    
      local a_dt_dx = 0.0 -- Mixed space-time derivative of metric term a.
      local a_dx_dx = 0.0 -- Second space derivative of metric term a.
      local b_dt_dx = 0.0 -- Mixed space-time derivative of metric term b.
      local b_dx_dx = (b_dx_plus - b_dx_minus) / math.pow(10.0, -8.0) -- Second space derivative of metric term b.
      local c_dt_dx = 0.0 -- Mixed space-time derivative of metric term c.
      local c_dx_dx = 0.0 -- Second space derivative of metric term c.
    
      local mom_x = 0.0 -- Fluid momentum (x-direction).
    
      return exp_2a,
        a_dt, a_dx,
        b_dt, b_dx,
        c_dt, c_dx,
        a_dt_dx, a_dx_dx,
        b_dt_dx, b_dx_dx,
        c_dt_dx, c_dx_dx,
        Etot, mom_x
    end,

    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (x-direction).
  }
}

-- Run application.
momentApp:run()
