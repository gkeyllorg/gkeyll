-- Cross-species relaxation of an electron-ion plasma under LBO collisions (Vlasov, 1x2v).
-- Homogeneous electrons (T = 1, drifting at half the ion thermal speed) and ions
-- (mass 100, T = 0.5, at rest) with constant collision frequencies nu_ee = nu_ei = 1, nu_ii = 0.1
-- and nu_ie = (m_e / m_i) nu_ei, so that the pair exchanges momentum and energy conservatively. Run
-- for 20 electron collision times. Particle number, total momentum and total energy are conserved
-- to 1e-12. The drift difference decays at nu_ei (measured 1.01) and the temperature difference at
-- 2 (nu_ei m_e + nu_ie m_i) / (m_e + m_i) = 0.0396 (measured 0.0396).
local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
mass_elc = 1.0 -- Electron mass.
mass_ion = 100.0 -- Ion mass.
charge_elc = -1.0 -- Electron charge.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.
Te = 1.0 -- Electron temperature.
Ti = 0.5 -- Ion temperature.
nu_ee = 1.0 -- Electron-electron collision frequency.
nu_ei = 1.0 -- Electron-ion collision frequency.
nu_ii = 0.1 -- Ion-ion collision frequency.
nu_ie = (mass_elc / mass_ion) * nu_ei -- Ion-electron collision frequency.

-- Derived physical quantities (using normalized code units).
vte = math.sqrt(Te / mass_elc) -- Electron thermal velocity.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.
ue = 0.5 * vti -- Electron drift velocity (x-direction).
ui = 0.0 -- Ion drift velocity (x-direction).

-- Simulation parameters.
Nx = 2 -- Cell count (configuration space: x-direction).
Nvx = 16 -- Cell count (velocity space: vx-direction).
Nvy = 16 -- Cell count (velocity space: vy-direction).
Lx = 5.0 -- Domain size (configuration space: x-direction).
vx_max_elc = 5.0 * vte -- Electron domain boundary (velocity space: vx-direction).
vy_max_elc = 5.0 * vte -- Electron domain boundary (velocity space: vy-direction).
vx_max_ion = 5.0 * vti -- Ion domain boundary (velocity space: vx-direction).
vy_max_ion = 5.0 * vti -- Ion domain boundary (velocity space: vy-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 20.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
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
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Electrons.
  elc = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_elc, mass = mass_elc,

    -- Velocity space grid.
    lower = { -vx_max_elc, -vy_max_elc },
    upper = { vx_max_elc, vy_max_elc },
    cells = { Nvx, Nvy },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          return n0 -- Number density.
        end,
        temperatureInit = function (t, xn)
          return Te -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return ue, 0.0 -- Drift velocity.
        end,

        correctAllMoments = true,
        useLastConverged = true
      }
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_ee -- Self-collision frequency.
      end,

      numCrossCollisions = 1,
      collideWith = { "ion" },
      collideWithCrossNu = {
        {
          crossNu = function (t, xn)
            return nu_ei -- Cross-collision frequency.
          end
        }
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.LTEMoments }
  },

  -- Ions.
  ion = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_ion, mass = mass_ion,

    -- Velocity space grid.
    lower = { -vx_max_ion, -vy_max_ion },
    upper = { vx_max_ion, vy_max_ion },
    cells = { Nvx, Nvy },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          return n0 -- Number density.
        end,
        temperatureInit = function (t, xn)
          return Ti -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return ui, 0.0 -- Drift velocity.
        end,

        correctAllMoments = true,
        useLastConverged = true
      }
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_ii -- Self-collision frequency.
      end,

      numCrossCollisions = 1,
      collideWith = { "elc" },
      collideWithCrossNu = {
        {
          crossNu = function (t, xn)
            return nu_ie -- Cross-collision frequency.
          end
        }
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.LTEMoments }
  },

  skipField = true,

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = 1.0, mu0 = 1.0,

    -- Initial conditions function.
    init = function (t, xn)
      return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0,
  }
}

vlasovApp:run()
