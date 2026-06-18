-- Burch et al. magnetic reconnection test for the 10-moment equations.
-- Input parameters match the initial conditions in Section 3, from the article:
-- J. M. TenBarge, J. Ng, J. Juno, L. Wang, A. M. Hakim and A. Bhattacharjee (2019), "An Extended MHD Study of the 16 October 2015 MMS Diffusion Region Crossing",
-- Journal of Geophysical Research: Space Physics, Volume 124 (11): 8474-8487.
-- https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JA026731

local Moments = G0.Moments
local TenMoment = G0.Moments.Eq.TenMoment

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_ion = 1.0 -- Proton mass.
charge_ion = 1.0 -- Proton charge.
mass_elc = mass_ion / 100.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
vAe = 0.2 -- Electron Alfven velocity.
n0 = 1.0 -- Reference number density.
beta2 = 2.748 -- Magnetosheath plasma beta.
Ti1_over_Ti2 = 7.73 / 1.374 -- Magnetospheric ion temperature / magnetosheath ion temperature.
Te1_over_Te2 = 1.288 / 1.374 -- Magnetospheric electron temperature / magnetosheath ion temperature.

-- Derived physical quantities (using normalized code units).
light_speed = 1.0 / math.sqrt(mu0 * epsilon0) -- Speed of light.
B0 = vAe * math.sqrt(n0 * mass_elc) -- Reference magnetic field strength.
omega_pi = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_ion)) -- Ion plasma frequency.
di = light_speed / omega_pi -- Ion skin depth.
omega0 = 1.0 * di -- Reference frequency.
psi0 = 0.1 * B0 * di -- Reference magnetic scalar potential.
guide1 = 0.099 * B0 -- Magnetospheric guide field strength.
guide2 = guide1 -- Magnetosheath guide field strength.
b1 = 1.696 * B0 -- Magnetospheric magnetic field strength.
b2 = 1.00 * B0 -- Magnetosheath magnetic field strength.
n1 = 0.062 * n0 -- Magnetospheric number density.
n2 = 1.0 * n0 -- Magnetosheath number density.
Ti2 = beta2 * (b2 * b2) / (2.0 * n2 * mu0) -- Magnetosheath ion temperature.

Te1 = Ti2 * Te1_over_Te2 -- Magnetospheric electron temperature.
Ti1 = Ti2 * Ti1_over_Ti2 -- Magnetospheric ion temperature.

Te2 = (0.5 * (b1 * b1 - b2 * b2) + 0.5 * (guide1 * guide1 - guide2 * guide2)
  + n1 * (Ti1 + Te1) - n2 * Ti2) / n2 -- Magnetosheath electron temperature (so that the system is in force balance).

Nx = 256 -- Cell count (x-direction).
Ny = 128 -- Cell count (y-direction).
Lx = 40.96 * di -- Domain size (x-direction).
Ly = 20.48 * di -- Domain size (y-direction).
k0_elc = 1.0 -- Closure parameter for electrons.
k0_ion = 0.1 -- Closure parameter for ions.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 250.0 -- Final simulation time.
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
  lower = { 0.0, 0.0 },
  upper = { Lx * di, Ly * di},
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
      k0 = k0_elc,
      hasNNClosure = use_nn_closure,
      polyOrder = poly_order,
      NNClosureFile = nn_closure_file,
      NNSpeciesName = "elc"
    },

    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local b1x = 0.5 * (b2 + b1) * (math.tanh((y - Ly * 0.25) / omega0) - math.tanh((y - Ly * 0.75) / omega0)
        + math.tanh((y - Ly * 1.25) / omega0) - math.tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (b2 - b1) -- Magnetospheric magnetic field (x-direction).
      local b1y = 0.0 -- Magnetospheric magnetic field (y-direction).
      local b1z = 0.5 * (guide2 - guide1) * (math.tanh((y - Ly * 0.25) / omega0) - math.tanh((y - Ly * 0.75) / omega0)
        + math.tanh((y - Ly * 1.25) / omega0) - math.tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (guide2 + guide1) -- Magnetospheric magnetic field (z-direction).
  
      local Ti_tot = 0.5 * (Ti2 - Ti1) * (math.tanh((y - Ly * 0.25) / omega0) - math.tanh((y - Ly * 0.75) / omega0)
        + math.tanh((y - Ly * 1.25) / omega0) - math.tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (Ti2 + Ti1) -- Total ion temperature.
      local Te_tot = 0.5 * (Te2 - Te1) * (math.tanh((y - Ly * 0.25) / omega0) - math.tanh((y - Ly * 0.75) / omega0)
        + math.tanh((y - Ly * 1.25) / omega0) - math.tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (Te2 + Te1) -- Total electron temperature.
      local n = (0.5 * (b1 * b1 - b1x * b1x) + 0.5 * (guide1 * guide1 - b1z * b1z) + n1 * (Ti1 + Te1)) / (Ti_tot + Te_tot) -- Total number density.
  
      local Bx = b1x - psi0 * 4.0 * pi / Ly * math.sin(2.0 * pi * x / Lx) * math.sin(4.0 * pi * y / Ly) -- Total magnetic field (x-direction).
      local By = b1y + psi0 * 2.0 * pi / Lx * math.cos(2.0 * pi * x / Lx) * (1.0 - math.cos(4.0 * pi * y / Ly)) -- Total magnetic field (y-direction).
      local Bz = b1z -- Total magnetic field (z-direction).
  
      local Te_frac = Te_tot / (Te_tot + Ti_tot) -- Fraction of total temperature from electrons.
      local Ti_frac = Ti_tot / (Te_tot + Ti_tot) -- Fraction of total temperature from ions
  
      local Jx = 0.5 * (guide2 - guide1) / omega0 * ((1.0 / math.cosh((y - Ly * 0.25) / omega0)) * (1.0 / math.cosh((y - Ly * 0.25) / omega0)) 
        - (1.0 / math.cosh((y - Ly * 0.75) / omega0)) * (1.0 / math.cosh((y - Ly * 0.75) / omega0)) 
        + (1.0 / math.cosh((y - Ly * 1.25) / omega0)) * (1.0 / math.cosh((y - Ly * 1.25) / omega0)) 
        - (1.0 / math.cosh((y + Ly * 0.25) / omega0)) * (1.0 / math.cosh((y + Ly * 0.25) / omega0))) -- Total current density (x-direction).
      local Jy = 0.0 -- Total current density (y-direction).
      local Jz  = -0.5 * (b2 + b1) / omega0 * ((1.0 / math.cosh((y - Ly * 0.25) / omega0)) * (1.0 / math.cosh((y - Ly * 0.25) / omega0)) 
        - (1.0 / math.cosh((y - Ly * 0.75) / omega0)) * (1.0 / math.cosh((y - Ly * 0.75) / omega0)) 
        + (1.0 / math.cosh((y - Ly * 1.25) / omega0)) * (1.0 / math.cosh((y - Ly * 1.25) / omega0)) 
        - (1.0 / math.cosh((y + Ly * 0.25) / omega0)) * (1.0 / math.cosh((y + Ly * 0.25) / omega0))) 
        - psi0 * math.sin(2.0 * pi * x / Lx) * ((2.0 * pi / Lx) * (2.0 * pi / Lx) * (1.0 - math.cos(4.0 * pi * y / Ly))
        + (4.0 * pi / Ly) * (4.0 * pi / Ly) * math.cos(4.0 * pi * y / Ly)) -- Total current density (z-direction).
     
      local Jxe = Jx * Te_frac -- Electron current density (x-direction).
      local Jye = Jy * Te_frac -- Electron current density (y-direction).
      local Jze = Jz * Te_frac -- Electron current density (z-direction).
  
      local rhoe = mass_elc * n -- Electron mass density.
      local mome_x = (mass_elc / charge_elc) * Jxe -- Electron momentum density (x-direction).
      local mome_y = (mass_elc / charge_elc) * Jye -- Electron momentum density (y-direction).
      local mome_z = (mass_elc / charge_elc) * Jze -- Electron momentum density (z-direction).

      local pre_xx = n * Te_tot + mome_x * mome_x / rhoe -- Electron pressure tensor (xx-component).
      local pre_xy = mome_x * mome_y / rhoe -- Electron pressure tensor (xy-component).
      local pre_xz = mome_x * mome_z / rhoe -- Electron pressure tensor (xz-component).
      local pre_yy = n * Te_tot + mome_y * mome_y / rhoe -- Electron pressure tensor (yy-component).
      local pre_yz = mome_y * mome_y / rhoe -- Electron pressure tensor (yz-component).
      local pre_zz = n * Te_tot + mome_z * mome_z / rhoe -- Electron pressure tensor (zz-component).
    
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
      local x, y = xn[1], xn[2]

      local b1x = 0.5 * (b2 + b1) * (math.tanh((y - Ly * 0.25) / omega0) - math.tanh((y - Ly * 0.75) / omega0)
        + math.tanh((y - Ly * 1.25) / omega0) - math.tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (b2 - b1) -- Magnetospheric magnetic field (x-direction).
      local b1y = 0.0 -- Magnetospheric magnetic field (y-direction).
      local b1z = 0.5 * (guide2 - guide1) * (math.tanh((y - Ly * 0.25) / omega0) - math.tanh((y - Ly * 0.75) / omega0)
        + math.tanh((y - Ly * 1.25) / omega0) - math.tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (guide2 + guide1) -- Magnetospheric magnetic field (z-direction).
  
      local Ti_tot = 0.5 * (Ti2 - Ti1) * (math.tanh((y - Ly * 0.25) / omega0) - math.tanh((y - Ly * 0.75) / omega0)
        + math.tanh((y - Ly * 1.25) / omega0) - math.tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (Ti2 + Ti1) -- Total ion temperature.
      local Te_tot = 0.5 * (Te2 - Te1) * (math.tanh((y - Ly * 0.25) / omega0) - math.tanh((y - Ly * 0.75) / omega0)
        + math.tanh((y - Ly * 1.25) / omega0) - math.tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (Te2 + Te1) -- Total electron temperature.
      local n = (0.5 * (b1 * b1 - b1x * b1x) + 0.5 * (guide1 * guide1 - b1z * b1z) + n1 * (Ti1 + Te1)) / (Ti_tot + Te_tot) -- Total number density.
  
      local Bx = b1x - psi0 * 4.0 * pi / Ly * math.sin(2.0 * pi * x / Lx) * math.sin(4.0 * pi * y / Ly) -- Total magnetic field (x-direction).
      local By = b1y + psi0 * 2.0 * pi / Lx * math.cos(2.0 * pi * x / Lx) * (1.0 - math.cos(4.0 * pi * y / Ly)) -- Total magnetic field (y-direction).
      local Bz = b1z -- Total magnetic field (z-direction).
  
      local Te_frac = Te_tot / (Te_tot + Ti_tot) -- Fraction of total temperature from electrons.
      local Ti_frac = Ti_tot / (Te_tot + Ti_tot) -- Fraction of total temperature from ions
  
      local Jx = 0.5 * (guide2 - guide1) / omega0 * ((1.0 / math.cosh((y - Ly * 0.25) / omega0)) * (1.0 / math.cosh((y - Ly * 0.25) / omega0)) 
        - (1.0 / math.cosh((y - Ly * 0.75) / omega0)) * (1.0 / math.cosh((y - Ly * 0.75) / omega0)) 
        + (1.0 / math.cosh((y - Ly * 1.25) / omega0)) * (1.0 / math.cosh((y - Ly * 1.25) / omega0)) 
        - (1.0 / math.cosh((y + Ly * 0.25) / omega0)) * (1.0 / math.cosh((y + Ly * 0.25) / omega0))) -- Total current density (x-direction).
      local Jy = 0.0 -- Total current density (y-direction).
      local Jz  = -0.5 * (b2 + b1) / omega0 * ((1.0 / math.cosh((y - Ly * 0.25) / omega0)) * (1.0 / math.cosh((y - Ly * 0.25) / omega0)) 
        - (1.0 / math.cosh((y - Ly * 0.75) / omega0)) * (1.0 / math.cosh((y - Ly * 0.75) / omega0)) 
        + (1.0 / math.cosh((y - Ly * 1.25) / omega0)) * (1.0 / math.cosh((y - Ly * 1.25) / omega0)) 
        - (1.0 / math.cosh((y + Ly * 0.25) / omega0)) * (1.0 / math.cosh((y + Ly * 0.25) / omega0))) 
        - psi0 * math.sin(2.0 * pi * x / Lx) * ((2.0 * pi / Lx) * (2.0 * pi / Lx) * (1.0 - math.cos(4.0 * pi * y / Ly))
        + (4.0 * pi / Ly) * (4.0 * pi / Ly) * math.cos(4.0 * pi * y / Ly)) -- Total current density (z-direction).
     
      local Jxi = Jx * Ti_frac -- Ion current density (x-direction).
      local Jyi = Jy * Ti_frac -- Ion current density (y-direction).
      local Jzi = Jz * Ti_frac -- Ion current density (z-direction).
  
      local rhoi = mass_ion * n -- Ion mass density.
      local momi_x = (mass_ion / charge_ion) * Jxi -- Ion momentum density (x-direction).
      local momi_y = (mass_ion / charge_ion) * Jyi -- Ion momentum density (y-direction).
      local momi_z = (mass_ion / charge_ion) * Jzi -- Ion momentum density (z-direction).

      local pri_xx = n * Ti_tot + momi_x * momi_x / rhoi -- Ion pressure tensor (xx-component).
      local pri_xy = momi_x * momi_y / rhoi -- Ion pressure tensor (xy-component).
      local pri_xz = momi_x * momi_z / rhoi -- Ion pressure tensor (xz-component).
      local pri_yy = n * Ti_tot + momi_y * momi_y / rhoi -- Ion pressure tensor (yy-component).
      local pri_yz = momi_y * momi_y / rhoi -- Ion pressure tensor (yz-component).
      local pri_zz = n * Ti_tot + momi_z * momi_z / rhoi -- Ion pressure tensor (zz-component).
    
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

      local b1x = 0.5 * (b2 + b1) * (math.tanh((y - Ly * 0.25) / omega0) - math.tanh((y - Ly * 0.75) / omega0)
        + math.tanh((y - Ly * 1.25) / omega0) - math.tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (b2 - b1) -- Magnetospheric magnetic field (x-direction).
      local b1y = 0.0 -- Magnetospheric magnetic field (y-direction).
      local b1z = 0.5 * (guide2 - guide1) * (math.tanh((y - Ly * 0.25) / omega0) - math.tanh((y - Ly * 0.75) / omega0)
        + math.tanh((y - Ly * 1.25) / omega0) - math.tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (guide2 + guide1) -- Magnetospheric magnetic field (z-direction).
  
      local Ti_tot = 0.5 * (Ti2 - Ti1) * (math.tanh((y - Ly * 0.25) / omega0) - math.tanh((y - Ly * 0.75) / omega0)
        + math.tanh((y - Ly * 1.25) / omega0) - math.tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (Ti2 + Ti1) -- Total ion temperature.
      local Te_tot = 0.5 * (Te2 - Te1) * (math.tanh((y - Ly * 0.25) / omega0) - math.tanh((y - Ly * 0.75) / omega0)
        + math.tanh((y - Ly * 1.25) / omega0) - math.tanh((y + Ly * 0.25) / omega0) + 1.0) + 0.5 * (Te2 + Te1) -- Total electron temperature.
      local n = (0.5 * (b1 * b1 - b1x * b1x) + 0.5 * (guide1 * guide1 - b1z * b1z) + n1 * (Ti1 + Te1)) / (Ti_tot + Te_tot) -- Total number density.

      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).
  
      local Bx = b1x - psi0 * 4.0 * pi / Ly * math.sin(2.0 * pi * x / Lx) * math.sin(4.0 * pi * y / Ly) -- Total magnetic field (x-direction).
      local By = b1y + psi0 * 2.0 * pi / Lx * math.cos(2.0 * pi * x / Lx) * (1.0 - math.cos(4.0 * pi * y / Ly)) -- Total magnetic field (y-direction).
      local Bz = b1z -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true -- Evolve field?
  }
}

-- Run application.
momentApp:run()
