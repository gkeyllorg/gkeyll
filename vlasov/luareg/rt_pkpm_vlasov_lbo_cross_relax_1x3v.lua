-- Cross-species relaxation of kinetic electrons (1x3v) and PKPM ions under LBO collisions, the
-- mixed counterpart of vlasov/creg/rt_pkpm_lbo_cross_relax_1x.c (both species PKPM) and
-- vlasov/creg/rt_vlasov_lbo_cross_relax_1x3v.c (both species kinetic).
-- Homogeneous neutral electrons (T = 1, drifting along the magnetic field at half the ion thermal
-- speed) and ions (mass 100, T = 0.5, at rest) with constant collision frequencies nu_ee = nu_ei =
-- 1, nu_ii = 0.1 and nu_ie = (m_e / m_i) nu_ei, so that m_e n_e nu_ei = m_i n_i nu_ie. Both species
-- start with T_perp / T_par = 1.3. Run for 5 electron collision times.
-- The cross operator of the kinetic species carries the exchange of the pair: it takes the moments
-- of the ions from the node values of their fluid state, and the fluid state of the ions takes the
-- opposite of the momentum and energy that operator gives the electrons, in the same stage. The
-- normalized distributions of the ions take their self operator at a rate that includes the
-- electrons.
-- Figures of merit (2 cells; 8 x 8 x 8 in the electron velocity on the hybrid basis, 48 in w for
-- the ions; 1226 steps, to t = 5): particle number is conserved to 1e-13 and the total momentum and
-- total energy to 4e-13 (both species kinetic on the same basis and velocity grid,
-- vlasov/luareg/rt_vlasov_lbo_cross_relax_1x3v_tensor_p1.lua: 2e-8 and 5e-10). The drift difference
-- decays at 1.0102 against nu_ei + nu_ie = 1.0100 and the temperature difference at 0.0395 against
-- 2 (nu_ei m_e + nu_ie m_i) / (m_e + m_i) = 0.0396. The temperatures end at 0.95544 (electrons) and
-- 0.54496 (ions); both PKPM: 0.95546 and 0.54495, both kinetic: 0.95545 and 0.54494.
local Vlasov = G0.Vlasov

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
vte = math.sqrt(Te / mass_elc) -- Electron thermal velocity.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.
ue = 0.5 * vti -- Electron drift velocity.
ui = 0.0 -- Ion drift velocity.

-- Simulation parameters.
Nx = 2 -- Cell count (configuration space: x-direction).
Nvx = 48 -- Cell count (ion velocity space: normalized parallel velocity w).
Nv_elc = 8 -- Cell count (electron velocity space: each direction).
Lx = 5.0 -- Domain size (configuration space: x-direction).
vx_max = 8.0 -- Domain boundary (ion velocity space: w = v_par/xi).
v_max_elc = 5.0 * vte -- Domain boundary (electron velocity space: each direction).
poly_order = 1 -- Polynomial order.
basis_type = "tensor" -- Basis function set (hybrid phase-space basis of the kinetic species).
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 5.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate integrated L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

vlasovApp = Vlasov.App.new {

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

  -- Electrons (kinetic).
  elc = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge, mass = mass_elc,

    -- Velocity space grid.
    lower = { -v_max_elc, -v_max_elc, -v_max_elc },
    upper = { v_max_elc, v_max_elc, v_max_elc },
    cells = { Nv_elc, Nv_elc, Nv_elc },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local vx, vy, vz = xn[2], xn[3], xn[4]

          local vt_par_sq = Tpar_over_T * vte * vte
          local vt_perp_sq = Tperp_over_T * vte * vte

          -- Bi-Maxwellian: parallel direction v_x (drift, T_par), perpendicular plane (v_y, v_z) (T_perp).
          local norm = n0 / (math.pow(2.0 * pi, 1.5) * vt_perp_sq * math.sqrt(vt_par_sq))

          return norm * math.exp(-(vx - ue) * (vx - ue) / (2.0 * vt_par_sq) - (vy * vy + vz * vz) / (2.0 * vt_perp_sq))
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

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

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Ions (PKPM).
  ion = Vlasov.PKPMSpecies.new {
    charge = charge, mass = mass_ion,

    kinetic = {
      -- Velocity space grid.
      lower = { -vx_max },
      upper = { vx_max },
      cells = { Nvx },

      -- Initial conditions.
      numInit = 1,
      projections = {
        {
          projectionID = G0.Projection.Func,

          init = function (t, xn)
            local w = xn[2]

            -- Bi-Maxwellian in w = v_par/xi: F has variance T_par/T and G = (T_perp/T) F.
            local F = math.exp(-w * w / (2.0 * Tpar_over_T)) / math.sqrt(2.0 * pi * Tpar_over_T)

            return F, Tperp_over_T * F, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
          end
        }
      },

      collisions = {
        collisionID = G0.Collisions.LBO,

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
      }
    },

    fluid = {
      -- Initial conditions.
      init = function (t, xn)
        local rho = mass_ion * n0 -- Ion mass density.
        local E = 0.5 * rho * ui * ui + 1.5 * n0 * Ti -- Ion total energy density.

        return rho, rho * ui, 0.0, 0.0, E
      end
    },

    evolve = true -- Evolve species?
  },

  -- Field.
  field = Vlasov.Field.new {
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

vlasovApp:run()
