-- Magnetic reconnection of a Harris sheet with a guide field (the GEM reconnection challenge setup)
-- with the PKPM system of equations (2x, normalized parallel velocity w = v_par/xi) and LBO
-- collisions.
-- Electrons and ions (mass ratio 25, v_te/c = 1/4, Ti/Te = 5, electron beta 1/6 in the in-plane
-- field B0, guide field B0, background density 0.2 n0) in a sheet of half-width 0.5 d_i, perturbed
-- by the flux function psi0 cos(2 pi x/Lx) cos(pi y/Ly) with psi0 = 0.1 B0 d_i, on a domain 8 pi
-- d_i by 4 pi d_i, periodic in x, with conducting walls that reflect the species in y. Position
-- maps cluster the cells at the X point (x cells from 0.15 to 1.85 of the uniform size) and at the
-- sheet (y cells from 0.3 to 1.7). Run to t = 22/Omega_ci.
-- Figures of merit (12 x 12 x 8 cells), with the reconnected flux psi = (1/2) int |B_y(x, y = 0)|
-- dx and the rate d psi/dt between frames: psi = 0.200 B0 d_i initially and 1.177 at t =
-- 22/Omega_ci. The rate is 0.02 to 0.03 B0 vA from 7 to 15/Omega_ci, rises from 15/Omega_ci on
-- (onset), peaks at 0.148 B0 vA between 19 and 20/Omega_ci and is 0.113 between 21 and 22/Omega_ci.
-- The onset time depends on the mesh at this resolution.

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
nb_over_n0 = 0.2 -- Background number density / reference number density.
vte_over_c = 1.0 / 4.0 -- Electron thermal velocity sqrt(Te/me) / speed of light.
beta_elc = 1.0 / 6.0 -- Electron plasma beta (in-plane field).
Te_over_Ti = 0.2 -- Electron temperature / ion temperature.
guide = 1.0 -- Guide field / in-plane field.

-- Derived physical quantities (using normalized code units).
light_speed = 1.0 / math.sqrt(epsilon0 * mu0) -- Speed of light.
Te = mass_elc * (vte_over_c * light_speed) * (vte_over_c * light_speed) -- Electron temperature.
Ti = Te / Te_over_Ti -- Ion temperature.
vte = math.sqrt(Te / mass_elc) -- Electron thermal velocity.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.
B0 = math.sqrt(2.0 * mu0 * n0 * Te / beta_elc) -- Asymptotic in-plane magnetic field strength.
Bg = guide * B0 -- Guide magnetic field strength.

omega_ci = charge_ion * B0 / mass_ion -- Ion cyclotron frequency (in-plane field).
omega_pi = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_ion)) -- Ion plasma frequency.
d_i = light_speed / omega_pi -- Ion inertial length.

w0 = 0.5 * d_i -- Current sheet width.
psi0 = 0.1 * B0 * d_i -- Flux function perturbation amplitude.

nu_elc = 0.01 * omega_ci -- Electron collision frequency.
nu_ion = 0.01 * omega_ci / math.sqrt(mass_ion) -- Ion collision frequency.

-- Simulation parameters.
Nx = 12 -- Cell count (configuration space: x-direction).
Ny = 12 -- Cell count (configuration space: y-direction).
Nvx = 8 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 8.0 * pi * d_i -- Domain size (configuration space: x-direction).
Ly = 4.0 * pi * d_i -- Domain size (configuration space: y-direction).
map_amp_x = 0.85 -- Amplitude of the position map clustering the cells at the X point.
map_amp_y = 0.7 -- Amplitude of the position map clustering the cells at the current sheet.
vx_max_elc = 6.0 -- Domain boundary (electron velocity space: w = v_par/xi).
vx_max_ion = 6.0 -- Domain boundary (ion velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 22.0 / omega_ci -- Final simulation time.
num_frames = 22 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate integrated L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

local function sech2(x)
  return 1.0 / (math.cosh(x) * math.cosh(x))
end

-- Harris sheet with image sheets at y = -/+ Ly (periodicity of the in-plane field), the
-- number density of pressure balance with uniform temperatures, and the out-of-plane
-- current of the sheet and of the flux perturbation psi0 cos(2 pi x/Lx) cos(pi y/Ly).
local function harris_bx(y)
  return -B0 * (math.tanh(y / w0) - math.tanh((y + Ly) / w0) - math.tanh((y - Ly) / w0))
end

local function harris_n(y)
  local T = Te + Ti
  local bx = harris_bx(y)
  return 0.5 * (B0 * B0 - bx * bx) / (mu0 * T) + nb_over_n0 * n0
end

local function total_jz(x, y)
  local kx = 2.0 * pi / Lx
  local ky = pi / Ly
  return (B0 / (w0 * mu0)) * (sech2(y / w0) - sech2((y + Ly) / w0) - sech2((y - Ly) / w0))
    - (psi0 / mu0) * (kx * kx + ky * ky) * math.cos(kx * x) * math.cos(ky * y)
end

pkpmApp = PKPM.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { -0.5 * Lx, -0.5 * Ly },
  upper = { 0.5 * Lx, 0.5 * Ly },
  cells = { Nx, Ny },

  -- Position maps x(eta) = eta - a (L/2 pi) sin(2 pi eta/L) clustering the cells at the X point
  -- (x-direction) and at the current sheet (y-direction): they fix the end points and have the
  -- Jacobian 1 - a cos(2 pi eta/L), 1 - a at the center and 1 + a at the boundary.
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
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

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

      local Te_frac = Te / (Te + Ti) -- Electron share of the current.

      local n = harris_n(y) -- Electron number density.
      local Jz = total_jz(x, y) -- Out-of-plane current density.

      local mom_x = 0.0 -- Electron total momentum density (x-direction).
      local mom_y = 0.0 -- Electron total momentum density (y-direction).
      local mom_z = mass_elc * Jz * Te_frac / charge_elc -- Electron total momentum density (z-direction).

      local rho = mass_elc * n -- Electron mass density.
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

    bcy = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcReflect }, -- Reflecting walls (y-direction).

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

      local Ti_frac = Ti / (Te + Ti) -- Ion share of the current.

      local n = harris_n(y) -- Ion number density.
      local Jz = total_jz(x, y) -- Out-of-plane current density.

      local mom_x = 0.0 -- Ion total momentum density (x-direction).
      local mom_y = 0.0 -- Ion total momentum density (y-direction).
      local mom_z = mass_ion * Jz * Ti_frac / charge_ion -- Ion total momentum density (z-direction).

      local rho = mass_ion * n -- Ion mass density.
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

    bcy = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcReflect }, -- Reflecting walls (y-direction).

    evolve = true -- Evolve species?
  },

  -- Field.
  field = PKPM.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local kx = 2.0 * pi / Lx
      local ky = pi / Ly

      local Bx = harris_bx(y) + psi0 * ky * math.cos(kx * x) * math.sin(ky * y) -- Total magnetic field (x-direction).
      local By = -psi0 * kx * math.sin(kx * x) * math.cos(ky * y) -- Total magnetic field (y-direction).
      local Bz = Bg -- Total magnetic field (z-direction).

      return 0.0, 0.0, 0.0, Bx, By, Bz, 0.0, 0.0
    end,

    bcy = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall }, -- Conducting walls (y-direction).

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()
