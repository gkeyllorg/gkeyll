-- Parallel-propagating firehose instability test for the 10-moment equations.
-- Input parameters match the initial conditions in Section 4.4, from the article:
-- M. W. Kunz, J. M. Stone and X-N. Bai (2014), "Pegasus: A new hybrid-kinetic particle-in-cell code for astrophysical plasma dynamics",
-- Journal of Computational Physics, Volume 259: 154-174.
-- https://www.sciencedirect.com/science/article/pii/S0021999113007973

local Moments = G0.Moments
local TenMoment = G0.Moments.Eq.TenMoment

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_ion = 1836.0 -- Proton mass.
charge_ion = 1.0 -- Proton charge.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

vAe = 0.0125 -- Electron Alfven velocity.
n0 = 1.0 -- Reference number density.

-- Derived physical quantities (using normalized code units).
light_speed = 1.0 / math.sqrt(mu0 * epsilon0) -- Speed of light.
B0 = vAe * math.sqrt(mu0 * n0 * mass_elc) -- Reference magnetic field strength.
beta = 300.0 / pi -- Trace proton plasma beta.
dbeta = 100.0 -- Parallel proton plasma beta - perpendicular proton plasma beta.
beta_par = beta + 2.0 * dbeta / 3.0 -- Parallel proton plasma beta.
beta_perp = beta - dbeta / 3.0 -- Perpendicular proton plasma beta.

vte = vAe * math.sqrt(beta) -- Electron thermal velocity.
Te = vte * vte * mass_elc / 2.0 -- Electron temperature.

Ti_par = vAe * vAe * (beta_par * mass_elc / 2.0) -- Parallel ion temperature.
Ti_perp = vAe * vAe * (beta_perp * mass_elc / 2.0) -- Perpendicular ion temperature.

omega_ci = charge_ion * B0 / mass_ion -- Ion cyclotron frequency.
omega_pe = math.sqrt(n0 * charge_elc * charge_elc / (epsilon0 * mass_elc)) -- Electron plasma frequency.
de = light_speed / omega_pe -- Electron skin depth.
omega_pi = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_ion)) -- Ion plasma frequency.
di = light_speed / omega_pi -- Ion skin depth.
lambdaD = vte / omega_pe -- Electron Debye length.

noise_amp = 1.0e-6 * B0 -- Noise level for perturbation.
mode_init = 1 -- Initial wave mode to perturb with noise.
mode_final = 48 -- Final wave mode to perturb with noise.

-- Simulation parameters.
Nx = 560 -- Cell count (x-direction).
Lx = 300.0 * di -- Domain size (x-direction).
k0_elc = 0.1 / de -- Closure parameter for electrons.
k0_ion = 0.1 / di -- Closure parameter for ions.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 10.0 / omega_ci -- Final simulation time.
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
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Electrons.
  elc = Moments.Species.new {
    charge = charge_elc, mass = mass_elc,
    equation = TenMoment.new {
      k0 = k0_elc,
      hasNNClosure = use_nn_closure,
      polyOrder = poly_order,
      NNClosureFile = nn_closure_file,
      NNSpeciesName = "elc"
    },

    -- Initial conditions function.
    init = function (t, xn)
      local rhoe = mass_elc * n0 -- Electron mass density.
      local mome_x = 0.0 -- Electron momentum density (x-direction).
      local mome_y = 0.0 -- Electron momentum density (y-direction).
      local mome_z = 0.0 -- Electron momentum density (z-direction).
    
      local pre_xx = n0 * Te + mome_x * mome_x / rhoe -- Electron pressure tensor (xx-component).
      local pre_xy = mome_x * mome_y / rhoe -- Electron pressure tensor (xy-component).
      local pre_xz = mome_x * mome_z / rhoe -- Electron pressure tensor (xz-component).
      local pre_yy = n0 * Te + mome_y * mome_y / rhoe -- Electron pressure tensor (yy-component).
      local pre_yz = mome_y * mome_y / rhoe -- Electron pressure tensor (yz-component).
      local pre_zz = n0 * Te + mome_z * mome_z / rhoe -- Electron pressure tensor (zz-component).
	 
      return rhoe, mome_x, mome_y, mome_z, pre_xx, pre_xy, pre_xz, pre_yy, pre_yz, pre_zz
    end,

    evolve = true -- Evolve species?
  },

  -- Ions.
  ion = Moments.Species.new {
    charge = charge_ion, mass = mass_ion,
    equation = TenMoment.new {
      k0 = k0_ion,
      hasNNClosure = use_nn_closure,
      polyOrder = poly_order,
      NNClosureFile = nn_closure_file,
      NNSpeciesName = "ion"
    },

    -- Initial conditions function.
    init = function (t, xn)
      local rhoi = mass_ion * n0 -- Ion mass density.
      local momi_x = 0.0 -- Ion momentum density (x-direction).
      local momi_y = 0.0 -- Ion momentum density (y-direction).
      local momi_z = 0.0 -- Ion momentum density (z-direction).
    
      local pri_xx = n0 * Ti_par + momi_x * momi_x / rhoi -- Ion pressure tensor (xx-component).
      local pri_xy = momi_x * momi_y / rhoi -- Ion pressure tensor (xy-component).
      local pri_xz = momi_x * momi_z / rhoi -- Ion pressure tensor (xz-component).
      local pri_yy = n0 * Ti_perp + momi_y * momi_y / rhoi -- Ion pressure tensor (yy-component).
      local pri_yz = momi_y * momi_y / rhoi -- Ion pressure tensor (yz-component).
      local pri_zz = n0 * Ti_perp + momi_z * momi_z / rhoi -- Ion pressure tensor (zz-component).
	 
      return rhoi, momi_x, momi_y, momi_z, pri_xx, pri_xy, pri_xz, pri_yy, pri_yz, pri_zz
    end,

    evolve = true -- Evolve species?
  },

  -- Field.
  field = Moments.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,
    mgnErrorSpeedFactor = 1.0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).
    
      local Bx = B0 -- Total magnetic field (x-direction).
      local By = 0.0
      local Bz = 0.0
    
      local alpha = noise_amp * Bx -- Applied amplitude.
      local kx = 2.0 * pi / Lx -- Wave number (x-direction).
    
      math.randomseed(0)
    
      for i = mode_init, mode_final do
        By = By - alpha * math.random() * math.sin(i * kx * x + 2.0 * pi * math.random()) -- Total magnetic field (y-direction).
        Bz = Bz - alpha * math.random() * math.sin(i * kx * x + 2.0 * pi * math.random()) -- Total magnetic field (z-direction).
      end

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true -- Evolve field?
  }
}

-- Run application.
momentApp:run()
