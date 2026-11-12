-- Cross-species relaxation of an electron-ion plasma under BGK collisions with the PKPM system of
-- equations (1x, normalized parallel velocity w = v_par/xi), implicit update: the backward-Euler
-- twin of rt_pkpm_bgk_cross_relax_1x.c at 50 times its collision frequencies.
-- Homogeneous neutral electrons (T = 1, drifting along the magnetic field at half the ion thermal
-- speed) and ions (mass 100, T = 0.5, at rest) with constant collision frequencies
-- nu_ee = nu_ei = 50, nu_ii = 5 and nu_ie = (m_e / m_i) nu_ei. Both species start with
-- T_perp / T_par = 1.3. Run for 50 electron collision times in 18 steps, with (nu_ee + nu_ei) dt
-- between 5.2 and 5.8.
-- The exchange of the pair and the relaxation of the distributions take one backward-Euler step in
-- each SSP-RK3 stage, over the time weight of the stage, so with z = nu dt a step multiplies a
-- relaxing quantity by R(z) = [1/3 + (2/3) (3/4 + (1/4)/(1 + z))/(1 + z/4)]/(1 + 2 z/3).
-- Particle number is conserved exactly and the total momentum and total energy to 7e-15; the
-- constraints hold to 1e-15. Each step multiplies the drift difference by R((nu_ei + nu_ie) dt) to
-- the round-off of that difference, the temperature difference by
-- R(2 (nu_ei m_e + nu_ie m_i) / (m_e + m_i) dt) to 8e-4 (the frictional heating) and the pressure
-- anisotropy of each species by R((nu_ss + nu_sr) dt) to 5e-3 (the stages combine distributions
-- normalized by different temperatures). The temperatures end at 0.78782 (electrons) and 0.71260
-- (ions). The same holds with the drift across the magnetic field.
local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
mass_ion = 100.0 -- Ion mass.
charge = 0.0 -- Species charge (neutral: no coupling to the field).

n0 = 1.0 -- Reference number density.
Te = 1.0 -- Electron temperature.
Ti = 0.5 -- Ion temperature.
alpha = 1.3 -- Ratio of perpendicular to parallel temperatures.
nu_ee = 50.0 -- Electron-electron collision frequency.
nu_ei = 50.0 -- Electron-ion collision frequency.
nu_ii = 5.0 -- Ion-ion collision frequency.
nu_ie = (mass_elc / mass_ion) * nu_ei -- Ion-electron collision frequency.
B0 = 1.0 -- Magnetic field along x (sets the field direction b).

-- Derived physical quantities (using normalized code units).
Tpar_over_T = 3.0 / (1.0 + 2.0 * alpha) -- Parallel temperature / temperature.
Tperp_over_T = alpha * Tpar_over_T -- Perpendicular temperature / temperature.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.
ue = 0.5 * vti -- Electron drift velocity.
ui = 0.0 -- Ion drift velocity.

-- Simulation parameters.
Nx = 2 -- Cell count (configuration space: x-direction).
Nvx = 48 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 5.0 -- Domain size (configuration space: x-direction).
vx_max = 8.0 -- Domain boundary (velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 1.0 -- Final simulation time.
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
  lower = { 0.0 },
  upper = { Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coordinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Electrons.
  elc = PKPM.Species.new {
    charge = charge, mass = mass_elc,

    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[2]

      -- Bi-Maxwellian in w = v_par/xi: F has variance T_par/T and G = (T_perp/T) F.
      local F = math.exp(-w * w / (2.0 * Tpar_over_T)) / math.sqrt(2.0 * pi * Tpar_over_T)

      return F, Tperp_over_T * F
    end,

    initFluid = function (t, xn)
      local rho = mass_elc * n0 -- Electron mass density.
      local E = 0.5 * rho * ue * ue + 1.5 * n0 * Te -- Electron total energy density.

      return rho, rho * ue, 0.0, 0.0, E
    end,

    collisions = {
      collisionID = G0.Collisions.BGK,
      useImplicitCollisionScheme = true,

      selfNu = function (t, xn)
        return nu_ee
      end,

      numCrossCollisions = 1,
      collideWith = { "ion" },
      collideWithCrossNu = {
        {
          crossNu = function (t, xn)
            return nu_ei
          end
        }
      }
    },

    evolve = true -- Evolve species?
  },

  -- Ions.
  ion = PKPM.Species.new {
    charge = charge, mass = mass_ion,

    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[2]

      -- Bi-Maxwellian in w = v_par/xi: F has variance T_par/T and G = (T_perp/T) F.
      local F = math.exp(-w * w / (2.0 * Tpar_over_T)) / math.sqrt(2.0 * pi * Tpar_over_T)

      return F, Tperp_over_T * F
    end,

    initFluid = function (t, xn)
      local rho = mass_ion * n0 -- Ion mass density.
      local E = 0.5 * rho * ui * ui + 1.5 * n0 * Ti -- Ion total energy density.

      return rho, rho * ui, 0.0, 0.0, E
    end,

    collisions = {
      collisionID = G0.Collisions.BGK,
      useImplicitCollisionScheme = true,

      selfNu = function (t, xn)
        return nu_ii
      end,

      numCrossCollisions = 1,
      collideWith = { "elc" },
      collideWithCrossNu = {
        {
          crossNu = function (t, xn)
            return nu_ie
          end
        }
      }
    },

    evolve = true -- Evolve species?
  },

  -- Field.
  field = PKPM.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      return 0.0, 0.0, 0.0, B0, 0.0, 0.0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()
