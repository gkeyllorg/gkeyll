-- Orszag-Tang vortex on a mesh that is nonuniform in x and y with the PKPM system of equations (2x,
-- normalized parallel velocity w = v_par/xi) and LBO collisions: rt_pkpm_ot_2x with position maps.
-- Electrons and ions (mass ratio 25, v_Ae/c = 0.5, beta 0.08 and equal temperatures) at uniform
-- density in a guide field B0 along z, with the in-plane field 0.2 B0 (-sin(2 pi y/L), sin(4 pi
-- x/L)) and flow 0.2 vA (-sin(2 pi y/L), sin(2 pi x/L)), the current carried by the electrons, on a
-- doubly periodic box of side L = 20.48 d_i. Both species have fourth-order momentum hyperdiffusion
-- (D = 1e-4). The position maps cluster the cells at x = 0 and y = 0 (cells from 0.5 to 1.5 of the
-- uniform size in both directions), periodically. Run to t = 25/Omega_ci, through the formation of
-- the first current sheets.
-- Figures of merit (16 x 16 x 8 cells), from the cell averages of the current J = sum_s (q_s/m_s)
-- (rho u)_s in units of n0 e c, with lengths in d_e, the uniform mesh of rt_pkpm_ot_2x in
-- parentheses: int |J|^2 dA is 0.929 (0.945) initially, 1.764 (1.883) at t = 20/Omega_ci and 1.873
-- (2.058) at 25/Omega_ci. The largest |J_z| rises from 0.0180 (0.0170) to a first maximum of 0.0385
-- (0.0372) at 21/Omega_ci, in the sheets at x = L/2, y = 0, and is 0.0360 (0.0348) at 25/Omega_ci.
-- The sheets form at several places, so the clustering lowers int |J|^2.

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
d_i = vAi / omega_ci -- Ion inertial length.

delta_B0 = 0.2 * B0 -- Reference magnetic field strength perturbation.
delta_u0 = 0.2 * vAi -- Reference fluid velocity perturbation.

nu_elc = 0.01 * omega_ci -- Electron collision frequency.
nu_ion = 0.01 * omega_ci / math.sqrt(mass_ion) -- Ion collision frequency.

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Ny = 16 -- Cell count (configuration space: y-direction).
Nvx = 8 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 20.48 * d_i -- Domain size (configuration space: x-direction).
Ly = 20.48 * d_i -- Domain size (configuration space: y-direction).
map_amp_x = 0.5 -- Amplitude of the position map (x-direction).
map_amp_y = 0.5 -- Amplitude of the position map (y-direction).
vx_max_elc = 6.0 -- Domain boundary (electron velocity space: w = v_par/xi).
vx_max_ion = 6.0 -- Domain boundary (ion velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 25.0 / omega_ci -- Final simulation time.
num_frames = 25 -- Number of output frames.
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

  -- Position maps x(eta) = eta - a (L/2 pi) sin(2 pi eta/L) clustering the cells at x = 0 and
  -- y = 0: periodic with the domain, with the Jacobian dx/deta = 1 - a cos(2 pi eta/L), 1 - a at
  -- the periodic boundary and 1 + a at the center of the domain.
  mapc2pPos = {
    {
      pmap = function (t, xn)
        local eta = xn[1]

        return eta - map_amp_x * (Lx / (2.0 * pi)) * math.sin(2.0 * pi * eta / Lx)
      end
    },
    {
      pmap = function (t, xn)
        local eta = xn[1]

        return eta - map_amp_y * (Ly / (2.0 * pi)) * math.sin(2.0 * pi * eta / Ly)
      end
    },
  },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coordinate direction (x-direction, y-direction).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1, 2 }, -- Periodic directions (x-direction, y-direction).

  -- Electrons.
  elc = PKPM.Species.new {
    charge = charge_elc, mass = mass_elc,

    -- Velocity space grid.
    lower = { -vx_max_elc },
    upper = { vx_max_elc },
    cells = { Nvx },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[3] -- Normalized parallel velocity w = v_par/xi.

      -- Normalized isotropic Maxwellian: F = exp(-w^2/2)/sqrt(2 pi) (M0 = 1, M2 = 1) and
      -- G = (T_perp/T) F = F, independent of position; the density, flow and temperature
      -- enter through the fluid state.
      local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * math.pi)

      return F, F
    end,

    initFluid = function (t, xn)
      local x, y = xn[1], xn[2]

      local Jz = (delta_B0 * (4.0 * pi / Lx) * math.cos(4.0 * pi * x / Lx) +
                  delta_B0 * (2.0 * pi / Ly) * math.cos(2.0 * pi * y / Ly)) /
                 mu0

      local vx_drift = -delta_u0 * math.sin(2.0 * pi * y / Ly) -- Electron drift velocity (x-direction).
      local vy_drift = delta_u0 * math.sin(2.0 * pi * x / Lx) -- Electron drift velocity (y-direction).
      local vz_drift = -Jz / charge_ion -- Electron drift velocity (z-direction).

      local mom_x = mass_elc * n0 * vx_drift -- Electron total momentum density (x-direction).
      local mom_y = mass_elc * n0 * vy_drift -- Electron total momentum density (y-direction).
      local mom_z = mass_elc * n0 * vz_drift -- Electron total momentum density (z-direction).

      local rho = mass_elc * n0 -- Electron mass density.
      local vtsq = vte * vte -- Electron T/m.
      local E = 0.5 * (mom_x * mom_x + mom_y * mom_y + mom_z * mom_z) / rho + 1.5 * rho * vtsq -- Electron total energy density.

      return rho, mom_x, mom_y, mom_z, E
    end,

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_elc
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
    charge = charge_ion, mass = mass_ion,

    -- Velocity space grid.
    lower = { -vx_max_ion },
    upper = { vx_max_ion },
    cells = { Nvx },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[3] -- Normalized parallel velocity w = v_par/xi.

      -- Normalized isotropic Maxwellian: F = exp(-w^2/2)/sqrt(2 pi) (M0 = 1, M2 = 1) and
      -- G = (T_perp/T) F = F, independent of position; the density, flow and temperature
      -- enter through the fluid state.
      local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * math.pi)

      return F, F
    end,

    initFluid = function (t, xn)
      local x, y = xn[1], xn[2]

      local vx_drift = -delta_u0 * math.sin(2.0 * pi * y / Ly) -- Ion drift velocity (x-direction).
      local vy_drift = delta_u0 * math.sin(2.0 * pi * x / Lx) -- Ion drift velocity (y-direction).
      local vz_drift = 0.0 -- Ion drift velocity (z-direction).

      local mom_x = mass_ion * n0 * vx_drift -- Ion total momentum density (x-direction).
      local mom_y = mass_ion * n0 * vy_drift -- Ion total momentum density (y-direction).
      local mom_z = mass_ion * n0 * vz_drift -- Ion total momentum density (z-direction).

      local rho = mass_ion * n0 -- Ion mass density.
      local vtsq = vti * vti -- Ion T/m.
      local E = 0.5 * (mom_x * mom_x + mom_y * mom_y + mom_z * mom_z) / rho + 1.5 * rho * vtsq -- Ion total energy density.

      return rho, mom_x, mom_y, mom_z, E
    end,

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_ion
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

      local Jz = (delta_B0 * (4.0 * pi / Lx) * math.cos(4.0 * pi * x / Lx) +
                  delta_B0 * (2.0 * pi / Ly) * math.cos(2.0 * pi * y / Ly)) /
                 mu0

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
