-- Lower-Hybrid Drift Instability (LHDI) test for the 10-moment equations.
-- Input parameters match the initial conditions in Section 5.2, from the article:
-- J. Ng, A. Hakim, J. Juno and A. Bhattacharjee (2019), "Drift Instabilities in Thin Current Sheets Using a Two-Fluid Model With Pressure Tensor Effects",
-- Journal of Geophysical Research: Space Physics, Volume 124 (5): 3331-3346.
-- https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JA026313

local Moments = G0.Moments
local TenMoment = G0.Moments.Eq.TenMoment

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_ion = 36.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
Te_over_Ti = 0.1 -- Electron temperature / ion temperature.

n0 = 1.0 -- Reference number density.
nb_over_n0 = 0.001 -- Background number density / reference number density.
vte = 0.06 -- Electron thermal velocity.

beta = 1.0 / 11.0 -- Electron plasma beta.

noise_amp = 0.0001 -- Noise level for perturbation.
mode = 8.0 -- Wave mode to perturb with noise.

-- Derived physical quantities (using normalized code units).
Te = vte * vte * mass_elc / 2.0 -- Electron temperature.
Ti = Te / Te_over_Ti -- Ion temperature.
vti = math.sqrt(2.0 * Ti / mass_ion) -- Ion thermal velocity.

vAe = vte / math.sqrt(beta) -- Electron Alfven velocity.
B0 = vAe -- Reference magnetic field strength (derived from normalization of mass_elc and n0).
vAi = vAe / math.sqrt(mass_ion) -- Ion Alfven velocity.

omega_ci = charge_ion * B0 / mass_ion -- Ion cyclotron frequency.
omega_ce = charge_ion * B0 / mass_elc -- Electron cyclotron frequency.

larmor_ion = vti / omega_ci -- Ion Larmor radius.
larmor_elc = vte / omega_ce -- Electron Larmor radius.

l = larmor_ion -- Current sheet width.

nb = n0 * nb_over_n0 -- Background number density.
Te_frac = Te / (Te + Ti) -- Fraction of total temperature from electrons.
Ti_frac = 1.0 - Te_frac -- Fraction of total temperature from ions.
ix = mode -- Current (x-direction).
iy = 1.0 -- Current (y-direction).

-- Simulation parameters.
Nx = 64 -- Cell count (x-direction).
Ny = 128 -- Cell count (y-direction).
Lx = 6.4 * l -- Domain size (x-direction).
Ly = 12.8 * l -- Domain size (y-direction).
k0_elc = 1.0 -- Closure parameter for electrons.
k0_ion = 1.0 / 6.0 -- Closure parameter for ions.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 1100.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Neural network parameters.
use_nn_closure = false -- Use neural network-based closure?
poly_order = 1 -- Polynomial order of learned DG coefficients.
nn_closure_file = "moments/data/neural_nets/pkpm_ot_p1_moms_nn_1" -- File path of neural network to use.

momentApp = Moments.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { -0.5 * Lx, -0.5 * Ly },
  upper = { 0.5 * Lx, 0.5 * Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction (x- and y-directions).
  
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
      local x, y = xn[1], xn[2]

      local sech_sq = (1.0 / math.cosh(y / l)) * (1.0 / math.cosh(y / l)) -- Hyperbolic secant squared.

      local n = n0 * sech_sq -- Total number density.
      local Jx_noise = -noise_amp * (iy * pi / Ly) * math.sin(iy * pi * y / Ly) * math.sin(ix * 2.0 * pi* x / Lx) / mode -- Current density noise (x-direction).
      local Jy_noise = -noise_amp * (ix * 2.0 * pi / Lx) * math.cos(iy * pi * y / Ly) * math.cos(ix * 2.0 * pi * x / Lx) / mode -- Current density noise (y-direction).
    
      local Jx  = (B0 / l) * (-sech_sq) + Jx_noise -- Total current density, with noise (x-direction).
      local Jy  = Jy_noise -- Total current density, with noise (y-direction).
    
      local rhoe = n * mass_elc -- Electron mass density.
      local mome_x = (mass_elc / charge_elc) * Jx * Te_frac -- Electron momentum density (x-direction).
      local mome_y = (mass_elc / charge_elc) * Jy * Te_frac -- Electron momentum density (y-direction).
      local mome_z = 0.0 -- Electron momentum density (z-direction).
      local pre = n * Te -- Electron pressure (scalar).

      local pre_xx = pre + (mome_x * mome_x) / rhoe -- Electron pressure tensor (xx-component).
      local pre_xy = (mome_x * mome_y) / rhoe -- Electron pressure tensor (xy-component).
      local pre_xz = 0.0 -- Electron pressure tensor (xz-component).
      local pre_yy = pre + (mome_y * mome_y) / rhoe -- Electron pressure tensor (yy-component).
      local pre_yz = 0.0 -- Electron pressure tensor (yz-component).
      local pre_zz = pre -- Electron pressure tensor (zz-component).
	 
      return rhoe, mome_x, mome_y, mome_z, pre_xx, pre_xy, pre_xz, pre_yy, pre_yz, pre_zz
    end,

    evolve = true, -- Evolve species?
    bcy = { G0.SpeciesBc.bcWall, G0.SpeciesBc.bcWall } -- Wall boundary conditions (y-direction).
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
      local x, y = xn[1], xn[2]

      local sech_sq = (1.0 / math.cosh(y / l)) * (1.0 / math.cosh(y / l)) -- Hyperbolic secant squared.

      local n = n0 * sech_sq -- Total number density.
      local Jx_noise = -noise_amp * (iy * pi / Ly) * math.sin(iy * pi * y / Ly) * math.sin(ix * 2.0 * pi* x / Lx) / mode -- Current density noise (x-direction).
      local Jy_noise = -noise_amp * (ix * 2.0 * pi / Lx) * math.cos(iy * pi * y / Ly) * math.cos(ix * 2.0 * pi * x / Lx) / mode -- Current density noise (y-direction).
    
      local Jx  = (B0 / l) * (-sech_sq) + Jx_noise -- Total current density, with noise (x-direction).
      local Jy  = Jy_noise -- Total current density, with noise (y-direction).
    
      local rhoi = n * mass_ion -- Ion mass density.
      local momi_x = (mass_ion / charge_ion) * Jx * Ti_frac -- Ion momentum density (x-direction).
      local momi_y = (mass_ion / charge_ion) * Jy * Ti_frac -- Ion momentum density (y-direction).
      local momi_z = 0.0 -- Ion momentum density (z-direction).
      local pri = n * Ti -- Ion pressure (scalar).

      local pri_xx = pri + (momi_x * momi_x) / rhoi -- Ion pressure tensor (xx-component).
      local pri_xy = (momi_x * momi_y) / rhoi -- Ion pressure tensor (xy-component).
      local pri_xz = 0.0 -- Ion pressure tensor (xz-component).
      local pri_yy = pri + (momi_y * momi_y) / rhoi -- Ion pressure tensor (yy-component).
      local pri_yz = 0.0 -- Ion pressure tensor (yz-component).
      local pri_zz = pri -- Ion pressure tensor (zz-component).
	 
      return rhoi, momi_x, momi_y, momi_z, pri_xx, pri_xy, pri_xz, pri_yy, pri_yz, pri_zz
    end,

    evolve = true, -- Evolve species?
    bcy = { G0.SpeciesBc.bcWall, G0.SpeciesBc.bcWall } -- Wall boundary conditions (y-direction).
  },

  -- Field.
  field = Moments.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,
    mgnErrorSpeedFactor = 1.0,

    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local Bz_noise = noise_amp * math.cos(iy * pi * y / Ly) * math.sin(ix * 2.0 * pi * x / Lx) / mode

      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = -B0 * math.tanh(y / l) + Bz_noise; -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
    bcy = { G0.FieldBc.bcWall, G0.FieldBc.bcWall } -- Wall boundary conditions (y-direction).
  }
}

-- Run application.
momentApp:run()
