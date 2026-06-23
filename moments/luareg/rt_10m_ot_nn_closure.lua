local Moments = G0.Moments
local TenMoment = G0.Moments.Eq.TenMoment

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 25.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.
vAe = 0.5 -- Electron Alfven velocity.
beta = 0.08 -- Plasma beta.

-- Derived physical quantities (using normalized code units).
B0 = vAe * math.sqrt(mu0 * n0 * mass_elc) -- Reference magnetic field strength.
vAi = vAe / math.sqrt(mass_ion) -- Ion Alfven velocity.

vte = vAe * math.sqrt(beta / 2.0) -- Electron thermal velocity.
vti = vte / math.sqrt(mass_ion) -- Ion thermal velocity.

omega_ci = charge_ion * B0 / mass_ion -- Ion cyclotron frequency.
d_i = vAi / omega_ci -- Ion sound inertial length.

delta_B0 = 0.2 * B0 -- Reference magnetic field strength perturbation.
delta_u0 = 0.2 * vAi -- Reference fluid velocity perturbation.

-- Simulation parameters.
Nx = 128 -- Cell count (x-direction).
Ny = 128 -- Cell count (y-direction).
Lx = 20.48 * d_i -- Domain size (x-direction).
Ly = 20.48 * d_i -- Domain size (y-direction).
k0 = 5000.0 -- Closure parameter.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 75.0 / omega_ci -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Neural network parameters.
use_nn_closure = true -- Use neural network-based closure?
poly_order = 1 -- Polynomial order of learned DG coefficients.
nn_closure_file = "moments/data/neural_nets/pkpm_ot_p1_moms_nn_1" -- File path of neural network to use.

momentApp = Moments.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.0 },
  upper = { Lx, Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction (x- and y-directions).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1, 2 }, -- Periodic directions (x- and y-directions only).

  -- Electrons.
  elc = Moments.Species.new {
    charge = charge_elc, mass = mass_elc,
    equation = TenMoment.new {
      k0 = k0,
      hasNNClosure = use_nn_closure,
      polyOrder = poly_order,
      NNClosureFile = nn_closure_file,
      NNSpeciesName = "elc"
    },

    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local Jz = (delta_B0 * (4.0 * pi / Lx) * math.cos(4.0 * pi * x / Lx) + delta_B0 * (2.0 * pi / Ly) *  math.cos(2.0 * pi * y / Ly)) / mu0

      local vx_drift = -delta_u0 * math.sin(2.0 * pi * y / Ly) -- Electron drift velocity (x-direction).
      local vy_drift = delta_u0 * math.sin(2.0 * pi * x / Lx) -- Electron drift velocity (y-direction).
      local vz_drift = -Jz / charge_ion -- Electron drift velocity (z-direction).

      local rhoe = n0 * mass_elc -- Electron mass density.
      local mome_x = mass_elc * vx_drift -- Electron total momentum density (x-direction).
      local mome_y = mass_elc * vy_drift -- Electron total momentum density (y-direction).
      local mome_z = mass_elc * vz_drift -- Electron total momentum density (z-direction).
      local pre = vte * vte * rhoe -- Electron pressure (scalar).

      local pre_xx = pre + (mome_x * mome_x) / rhoe -- Electron pressure tensor (xx-component).
      local pre_xy = 0.0 -- Electron pressure tensor (xy-component).
      local pre_xz = 0.0 -- Electron pressure tensor (xz-component).
      local pre_yy = pre + (mome_y * mome_y) / rhoe -- Electron pressure tensor (yy-component).
      local pre_yz = 0.0 -- Electron pressure tensor (yz-component).
      local pre_zz = pre + (mome_z * mome_z) / rhoe -- Electron pressure tensor (zz-component).
	 
      return rhoe, mome_x, mome_y, mome_z, pre_xx, pre_xy, pre_xz, pre_yy, pre_yz, pre_zz
    end,

    evolve = true -- Evolve species?
  },

  -- Ions.
  ion = Moments.Species.new {
    charge = charge_ion, mass = mass_ion,
    equation = TenMoment.new {
      k0 = k0,
      hasNNClosure = use_nn_closure,
      polyOrder = poly_order,
      NNClosureFile = nn_closure_file,
      NNSpeciesName = "ion"
    },

    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local vx_drift = -delta_u0 * math.sin(2.0 * pi * y / Ly) -- Ion drift velocity (x-direction).
      local vy_drift = delta_u0 * math.sin(2.0 * pi * x / Lx) -- Ion drift velocity (y-direction).
      local vz_drift = 0.0 -- Ion drift velocity (z-direction).

      local rhoi = n0 * mass_ion -- Ion mass density.
      local momi_x = mass_ion * vx_drift -- Ion total momentum density (x-direction).
      local momi_y = mass_ion * vy_drift -- Ion total momentum density (y-direction).
      local momi_z = mass_ion * vz_drift -- Ion total momentum density (z-direction).
      local pri = vti * vti * rhoi -- Ion pressure (scalar).

      local pri_xx = pri + (momi_x * momi_x) / rhoi -- Ion pressure tensor (xx-component).
      local pri_xy = 0.0 -- Ion pressure tensor (xy-component).
      local pri_xz = 0.0 -- Ion pressure tensor (xz-component).
      local pri_yy = pri + (momi_y * momi_y) / rhoi -- Ion pressure tensor (yy-component).
      local pri_yz = 0.0 -- Ion pressure tensor (yz-component).
      local pri_zz = pri + (momi_z * momi_z) / rhoi -- Ion pressure tensor (zz-component).
	 
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
      local x, y = xn[1], xn[2]

      local Jz = (delta_B0 * (4.0 * pi / Lx) * math.cos(4.0 * pi * x / Lx) + delta_B0 * (2.0 * pi / Ly) * math.cos(2.0 * pi * y / Ly)) / mu0

      local Bx = -delta_B0 * math.sin(2.0 * pi * y / Ly) -- Total magnetic field (x-direction).
      local By = delta_B0 * math.sin(4.0 * pi * x / Lx) -- Total magnetic field (y-direction).
      local Bz = B0 -- Total magnetic field (z-direction).

      local vx_drift = -delta_u0 * math.sin(2.0 * pi * y / Ly) -- Electron drift velocity (x-direction).
      local vy_drift = delta_u0 * math.sin(2.0 * pi * x / Lx) -- Electron drift velocity (y-direction).
      local vz_drift = -Jz / charge_ion -- Electron drift velocity (z-direction).

      local Ex = -((vy_drift * Bz) - (vz_drift * By)) -- Total electric field (x-direction).
      local Ey = -((vz_drift * Bx) - (vx_drift * Bz)) -- Total electric field (y-direction).
      local Ez = -((vx_drift * By) - (vy_drift * Bx)) -- Total electric field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true -- Evolve field?
  }
}

-- Run application.
momentApp:run()