-- Cross-species relaxation of an electron-ion plasma under BGK collisions (Vlasov, 1x3v).
-- Homogeneous electrons (T = 1, drifting at half the ion thermal speed) and ions
-- (mass 100, T = 0.5, at rest) with constant collision frequencies nu_ee = nu_ei = 1, nu_ii = 0.1
-- and nu_ie = (m_e / m_i) nu_ei, so that the pair exchanges momentum and energy conservatively.
-- Both species start with T_perp / T_par = 1.3 (parallel along v_x). Run for 5 electron collision
-- times. The BGK target of each species moves with the flow of its partner, so the relative drift
-- must stay below the ion thermal speed for the ion target temperature to remain positive. Particle
-- number, total momentum and total energy are conserved to 1e-8. The drift difference decays at
-- nu_ei (measured 1.02) and the temperature difference at 2 (nu_ei m_e + nu_ie m_i) / (m_e + m_i) =
-- 0.0396 (measured 0.0395). The anisotropy decays at nu_ee + nu_ei = 2 for the electrons and at
-- nu_ii + nu_ie = 0.11 for the ions (measured 2.09 and 0.110), half the LBO rates.
local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
mass_elc = 1.0 -- Electron mass.
mass_ion = 100.0 -- Ion mass.
charge_elc = -1.0 -- Electron charge.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.
Te = 1.0 -- Electron temperature.
Ti = 0.5 -- Ion temperature.
alpha = 1.3 -- Ratio of perpendicular to parallel temperatures.
Te_par = 3.0 * Te / (1.0 + 2.0 * alpha) -- Parallel electron temperature.
Te_perp = alpha * Te_par -- Perpendicular electron temperature.
Ti_par = 3.0 * Ti / (1.0 + 2.0 * alpha) -- Parallel ion temperature.
Ti_perp = alpha * Ti_par -- Perpendicular ion temperature.
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
Nvx = 8 -- Cell count (velocity space: vx-direction).
Nvy = 8 -- Cell count (velocity space: vy-direction).
Nvz = 8 -- Cell count (velocity space: vz-direction).
Lx = 5.0 -- Domain size (configuration space: x-direction).
vx_max_elc = 5.0 * vte -- Electron domain boundary (velocity space: vx-direction).
vy_max_elc = 5.0 * vte -- Electron domain boundary (velocity space: vy-direction).
vz_max_elc = 5.0 * vte -- Electron domain boundary (velocity space: vz-direction).
vx_max_ion = 5.0 * vti -- Ion domain boundary (velocity space: vx-direction).
vy_max_ion = 5.0 * vti -- Ion domain boundary (velocity space: vy-direction).
vz_max_ion = 5.0 * vti -- Ion domain boundary (velocity space: vz-direction).
poly_order = 1 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 5.0 -- Final simulation time.
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
    lower = { -vx_max_elc, -vy_max_elc, -vz_max_elc },
    upper = { vx_max_elc, vy_max_elc, vz_max_elc },
    cells = { Nvx, Nvy, Nvz },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local vx = xn[2]
          local vy = xn[3]
          local vz = xn[4]

          local vt_par_sq = Te_par / mass_elc
          local vt_perp_sq = Te_perp / mass_elc

          -- Bi-Maxwellian: parallel direction v_x (drift, T_par), perpendicular plane (v_y, v_z) (T_perp).
          local norm = n0 / (math.pow(2.0 * math.pi, 3.0 / 2.0) * vt_perp_sq * math.sqrt(vt_par_sq))
          local f = norm * math.exp(-(vx - ue) * (vx - ue) / (2.0 * vt_par_sq)
            - (vy * vy + vz * vz) / (2.0 * vt_perp_sq))

          return f
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.BGK,

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

    correct = {
      correctAllMoments = true,
      iterationEpsilon = 1.0e-12,
      maxIterations = 100,
      useLastConverged = true
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2ij, G0.Moment.LTEMoments }
  },

  -- Ions.
  ion = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_ion, mass = mass_ion,

    -- Velocity space grid.
    lower = { -vx_max_ion, -vy_max_ion, -vz_max_ion },
    upper = { vx_max_ion, vy_max_ion, vz_max_ion },
    cells = { Nvx, Nvy, Nvz },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local vx = xn[2]
          local vy = xn[3]
          local vz = xn[4]

          local vt_par_sq = Ti_par / mass_ion
          local vt_perp_sq = Ti_perp / mass_ion

          -- Bi-Maxwellian: parallel direction v_x (drift, T_par), perpendicular plane (v_y, v_z) (T_perp).
          local norm = n0 / (math.pow(2.0 * math.pi, 3.0 / 2.0) * vt_perp_sq * math.sqrt(vt_par_sq))
          local f = norm * math.exp(-(vx - ui) * (vx - ui) / (2.0 * vt_par_sq)
            - (vy * vy + vz * vz) / (2.0 * vt_perp_sq))

          return f
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.BGK,

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

    correct = {
      correctAllMoments = true,
      iterationEpsilon = 1.0e-12,
      maxIterations = 100,
      useLastConverged = true
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2ij, G0.Moment.LTEMoments }
  },

  skipField = true,
}

vlasovApp:run()
