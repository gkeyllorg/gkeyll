-- Magnetic reconnection of a force-free current sheet with the PKPM system of equations (2x,
-- normalized parallel velocity w = v_par/xi) and LBO collisions.
-- Electrons and ions (mass ratio 25, v_te/c = 1/4, Ti/Te = 5, electron beta 1/6 in the in-plane
-- field B0) at uniform density and temperatures in the field B = (-B0 tanh(y/w0), 0, B0 sqrt(1 +
-- sech^2(y/w0))) of uniform strength, w0 = 0.5 d_i, with the current carried by the electrons,
-- perturbed by the flux function psi0 cos(2 pi x/Lx) cos(pi y/Ly) with psi0 = 0.1 B0 d_i, on a
-- domain 8 pi d_i by 4 pi d_i, periodic in x, with conducting walls that reflect the species in y.
-- Position maps cluster the cells at the X point and at the sheet (cells from 0.15 to 1.85 of the
-- uniform size in both directions). Run to t = 22/Omega_ci.
-- Figures of merit (12 x 12 x 8 cells), with the reconnected flux psi = (1/2) int |B_y(x, y = 0)|
-- dx and the rate d psi/dt between frames: psi = 0.200 B0 d_i initially and 0.894 at t =
-- 22/Omega_ci. The rate is 0.02 to 0.03 B0 vA from 5 to 17/Omega_ci, rises from 17/Omega_ci on
-- (onset) and is 0.086 B0 vA between 20 and 21/Omega_ci and 0.088 between 21 and 22/Omega_ci. The
-- onset time depends on the mesh at this resolution.

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

n0 = 1.0 -- Number density (uniform).
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
map_amp_y = 0.85 -- Amplitude of the position map clustering the cells at the current sheet.
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

-- Force-free sheet: B = (-B0 tanh(y/w0), 0, B0 sqrt(guide^2 + sech^2(y/w0))) with |B| uniform,
-- uniform density and temperatures, the current J = curl B / mu0 carried by the electrons,
-- plus the flux perturbation psi0 cos(2 pi x/Lx) cos(pi y/Ly).
local function total_current(x, y)
  local kx = 2.0 * pi / Lx
  local ky = pi / Ly
  local s2 = sech2(y / w0)

  local Jx = -(B0 / (w0 * mu0)) * math.tanh(y / w0) * s2 / math.sqrt(s2 + guide * guide)
  local Jy = 0.0
  local Jz = (B0 / (w0 * mu0)) * s2 - (psi0 / mu0) * (kx * kx + ky * ky) * math.cos(kx * x) * math.cos(ky * y)

  return Jx, Jy, Jz
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

      local Jx, Jy, Jz = total_current(x, y) -- Current density, carried by the electrons.

      local mom_x = mass_elc * Jx / charge_elc -- Electron total momentum density (x-direction).
      local mom_y = mass_elc * Jy / charge_elc -- Electron total momentum density (y-direction).
      local mom_z = mass_elc * Jz / charge_elc -- Electron total momentum density (z-direction).

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
      local mom_x = 0.0 -- Ion total momentum density (x-direction).
      local mom_y = 0.0 -- Ion total momentum density (y-direction).
      local mom_z = 0.0 -- Ion total momentum density (z-direction).

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

      local Bx = -B0 * math.tanh(y / w0) + psi0 * ky * math.cos(kx * x) * math.sin(ky * y) -- Total magnetic field (x-direction).
      local By = -psi0 * kx * math.sin(kx * x) * math.cos(ky * y) -- Total magnetic field (y-direction).
      local Bz = B0 * math.sqrt(guide * guide + sech2(y / w0)) -- Total magnetic field (z-direction).

      return 0.0, 0.0, 0.0, Bx, By, Bz, 0.0, 0.0
    end,

    bcy = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall }, -- Conducting walls (y-direction).

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()
