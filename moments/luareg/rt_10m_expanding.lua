local Moments = G0.Moments
local TenMoment = G0.Moments.Eq.TenMoment

-- Physical constants (using normalized code units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.
U0 = 1.0 -- (Initial) comoving plasma velocity.
R0 = 1.0 -- (Initial) radial distance from expansion/contraction center.

rho = 1.0 -- Fluid mass density.
u = 0.0 -- Fluid velocity.
p = 1.0 -- Fluid pressure.

-- Simulation parameters.
Nx = 512 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
k0 = 0.1 -- Closure parameter.
cfl_frac = 0.95 -- CFL coefficient.

t_end = 0.5 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Neural network parameters.
use_nn_closure = false -- Use neural network-based closure?
poly_order = 1 -- Polynomial order of learned DG coefficients.
nn_closure_file = "moments/data/neural_nets/pkpm_periodic_es_shock_p1_moms_nn_1" -- File path of neural network to use.

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
    equation = TenMoment.new {
      k0 = k0,
      hasNNClosure = use_nn_closure,
      polyOrder = poly_order,
      NNClosureFile = nn_closure_file,
      NNSpeciesName = "elc"
    },

    hasVolumeSources = true,
    volumeGasGamma = gas_gamma,
    volumeU0 = U0,
    volumeR0 = R0,
    
    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]
  
      local mom_x = rho * u -- Fluid momentum density (x-direction).
      local mom_y = 0.0 -- Fluid momentum density (y-direction).
      local mom_z = 0.0 -- Fluid momentum density (z-direction).

      local p_xx = p + (0.5 * rho * u * u)-- Fluid pressure tensor (xx-component).
      local p_xy = 0.0 -- Fluid pressure tensor (xy-component).
      local p_xz = 0.0 -- Fluid pressure tensor (xz-component).
      local p_yy = p -- Fluid pressure tensor (yy-component).
      local p_yz = 0.0 -- Fluid pressure tensor (yz-component).
      local p_zz = p -- Fluid pressure tensor (zz-component).
      
      return rho, mom_x, mom_y, mom_z, p_xx, p_xy, p_xz, p_yy, p_yz, p_zz
    end,
  
    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (x-direction).
  },
    
  -- Field.
  field = Moments.Field.new {
    epsilon0 = 1.0, mu0 = 1.0,

    -- Initial conditions function.
    init = function (t, xn)
      return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    bcx = { G0.FieldBc.bcCopy, G0.FieldBc.bcCopy } -- Copy boundary conditions (x-direction).
  }
}

-- Run application.
momentApp:run()
