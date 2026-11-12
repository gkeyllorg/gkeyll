-- Cross-species relaxation of an electron-ion plasma under BGK collisions with the PKPM system of
-- equations (1x, normalized parallel velocity w = v_par/xi), the counterpart of
-- vlasov/creg/rt_vlasov_bgk_cross_relax_1x3v.c.
-- Homogeneous neutral electrons (T = 1, drifting along the magnetic field at half the ion thermal
-- speed) and ions (mass 100, T = 0.5, at rest) with constant collision frequencies nu_ee = nu_ei = 1,
-- nu_ii = 0.1 and nu_ie = (m_e / m_i) nu_ei, so that m_e n_e nu_ei = m_i n_i nu_ie. Both species
-- start with T_perp / T_par = 1.3. Run for 5 electron collision times.
-- The pair exchanges momentum and energy through the fluid states, in one backward-Euler step after
-- each time step, as in rt_pkpm_lbo_cross_relax_1x.c, and the distributions relax to the unit
-- Maxwellian: no target temperature is formed, so the relative drift is not limited by the ion
-- thermal speed as it is for the Vlasov target.
-- Particle number is conserved exactly and the total momentum and total energy to 3e-14; the
-- constraints hold to 1e-14. Each step divides the drift difference by 1 + (nu_ei + nu_ie) dt, with
-- (nu_ei + nu_ie) dt up to 0.049, so its fitted decay rate is 0.9864 against nu_ei + nu_ie = 1.0100.
-- The temperature difference decays at 2 (nu_ei m_e + nu_ie m_i) / (m_e + m_i) = 0.0396 (measured
-- 0.0395). The pressure anisotropy decays at nu_ee + nu_ei = 2 for the electrons and at
-- nu_ii + nu_ie = 0.11 for the ions (measured 2.0001 and 0.1100), half the LBO rates. Vlasov 1x3v
-- measures 1.02, 0.0395, 2.09 and 0.110. The same figures hold with the drift across the magnetic
-- field.
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
nu_ee = 1.0 -- Electron-electron collision frequency.
nu_ei = 1.0 -- Electron-ion collision frequency.
nu_ii = 0.1 -- Ion-ion collision frequency.
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

t_end = 5.0 -- Final simulation time.
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
