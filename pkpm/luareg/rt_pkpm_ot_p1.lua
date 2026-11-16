local PKPM = G0.PKPM

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

nu_elc = 0.01 * omega_ci -- Electron collision frequency.
nu_ion = 0.01 * omega_ci / math.sqrt(mass_ion) -- Ion collision frequency.

-- Simulation parameters.
Nx = 128 -- Cell count (configuration space: x-direction).
Ny = 128 -- Cell count (configuration space: y-direction).
Nvx = 16 -- Cell count (velocity space: vx-direction).
Lx = 20.48 * d_i -- Domain size (configuration space: x-direction).
Ly = 20.48 * d_i -- Domain size (configuration space: y-direction).
vx_max_elc = 6.0 * vte -- Domain boundary (electron velocity space: vx-direction).
vx_max_ion = 6.0 * vti -- Domain boundary (ion velocity space: vx-direction).
poly_order = 1 -- Polynomial order.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 75.0 / omega_ci -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate integrated L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

pkpmApp = PKPM.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.0, 0.0 },
  upper = { Lx, Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  useExplicitSource = true,

  -- Decomposition for configuration space.
  decompCuts = { 4, 2 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = { 1, 2 }, -- Periodic directions (x- and y-directions only).

  -- Electrons.
  elc = PKPM.Species.new {
    modelID = G0.Model.Default,
    charge = charge_elc, mass = mass_elc,

    -- Velocity space grid.
    lower = { -vx_max_elc },
    upper = { vx_max_elc },
    cells = { Nvx },

    -- Initial conditions (distribution function).
    initDist = function (t, xn)
      local vx = xn[3]

      local F0 = (n0 / math.sqrt(2.0 * pi * vte * vte)) * (math.exp(-(vx * vx) / (2.0 * vte * vte))) -- Electron distribution function (F0).
      local G = (vte * vte) * F0 -- Electron distribution function (G).
    
      return F0, G
    end,

    -- Initial conditions (fluid).
    initFluid = function (t, xn)
      local x, y = xn[1], xn[2]

      local Jz = (delta_B0 * (4.0 * pi / Lx) * math.cos(4.0 * pi * x / Lx) + delta_B0 * (2.0 * pi / Ly) * math.cos(2.0 * pi * y / Ly)) / mu0

      local vx_drift = -delta_u0 * math.sin(2.0 * pi * y / Ly) -- Electron drift velocity (x-direction).
      local vy_drift = delta_u0 * math.sin(2.0 * pi * x / Lx) -- Electron drift velocity (y-direction).
      local vz_drift = -Jz / charge_ion -- Electron drift velocity (z-direction).
      
      local mom_x = mass_elc * vx_drift -- Electron total momentum density (x-direction).
      local mom_y = mass_elc * vy_drift -- Electron total momentum density (y-direction).
      local mom_z = mass_elc * vz_drift -- Electron total momentum density (z-direction).

      return mom_x, mom_y, mom_z
    end,

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_elc -- Electron collision frequency.
      end
    },

    diffusion = {
      D = 1.0e-4,
      order = 4
    },

    evolve = true -- Evolve species?
  },

  -- Ions.
  ion = PKPM.Species.new {
    modelID = G0.Model.Default,
    charge = charge_ion, mass = mass_ion,
    
    -- Velocity space grid.
    lower = { -vx_max_ion },
    upper = { vx_max_ion },
    cells = { Nvx },

    -- Initial conditions (distribution function).
    initDist = function (t, xn)
      local vx = xn[3]

      local F0 = (n0 / math.sqrt(2.0 * pi * vti * vti)) * (math.exp(-(vx * vx) / (2.0 * vti * vti))) -- Ion distribution function (F0).
      local G = (vti * vti) * F0 -- Ion distribution function (G).
      
      return F0, G
    end,

    -- Initial conditions (fluid).
    initFluid = function (t, xn)
      local x, y = xn[1], xn[2]

      local vx_drift = -delta_u0 * math.sin(2.0 * pi * y / Ly) -- Ion drift velocity (x-direction).
      local vy_drift = delta_u0 * math.sin(2.0 * pi * x / Lx) -- Ion drift velocity (y-direction).
      local vz_drift = 0.0 -- Ion drift velocity (z-direction).
      
      local mom_x = mass_ion * vx_drift -- Ion total momentum density (x-direction).
      local mom_y = mass_ion * vy_drift -- Ion total momentum density (y-direction).
      local mom_z = mass_ion * vz_drift -- Ion total momentum density (z-direction).

      return mom_x, mom_y, mom_z
    end,

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_ion -- Ion collision frequency.
      end
    },

    diffusion = {
      D = 1.0e-4,
      order = 4
    },

    evolve = true -- Evolve species?
  },

  -- Field.
  field = PKPM.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

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

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()