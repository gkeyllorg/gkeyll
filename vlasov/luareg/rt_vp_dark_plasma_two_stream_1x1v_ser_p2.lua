-- Two-stream instability of a charged dark-matter pair plasma with the Vlasov-Poisson system
-- including both electrostatics and self-gravity (1x1v). Electrons and positrons (m = 1, |q| = 1,
-- vt = 1) each consist of two counter-streaming Maxwellian beams at +-5 vt with density n0/2 each.
-- The counter-streaming beams drive the electrostatic two-stream instability on the
-- plasma-frequency time scale while the mass in all four beams self-gravitates with coupling
-- alpha_g = 4 pi G epsilon0 m^2/q^2 = 1e-2, so that lambda_J = lambda_D/sqrt(alpha_g) = 10 lambda_D
-- and omega_J = sqrt(alpha_g) omega_pe = 0.1 omega_pe. The electron beams are seeded with
-- random-amplitude, random-phase density noise of amplitude 1e-4 in modes 1-32 of a box of length
-- 20 pi lambda_J (fundamental k lambda_J = 0.1); the positrons are unperturbed. The run lasts 10
-- Jeans times.
-- Kinetic linear theory for the counter-streaming beams (electrons and positrons contribute
-- equally): box modes 1-30 are electrostatically unstable, the fastest being mode 18 (k lambda_D =
-- 0.18) with gamma = 0.473 omega_pe, so the electrostatic energy grows at up to 2 gamma = 0.946.
-- Figures of merit (serendipity p2, 64 x 128 cells; the Lua noise realization differs from the C
-- one): the electrostatic energy grows at 0.896 over t = 12-20 and peaks at 2.13e3 at t = 23.7; the
-- gravitational energy then grows at 0.146 over t = 40-85 and peaks at 3.28e5 at t = 91.3; electron
-- and positron numbers conserved to 3e-13 and total energy (kinetic plus epsilon0/2 int |grad
-- phi|^2 minus half the gravitational energy file, int |grad phi_g|^2/alpha_g) to 9e-6.

local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
alpha_g = 1.0e-2 -- Gravitational coupling 4 pi G epsilon0 m^2/q^2.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_pos = 1.0 -- Positron mass.
charge_pos = 1.0 -- Positron charge.

n0 = 1.0 -- Reference number density.
T = 1.0 -- Temperature.
Vx_drift = 5.0 -- Beam drift velocity (x-direction).

delta_n = 1.0e-4 -- Applied perturbation amplitude.
mode_init = 1 -- Initial wave mode to perturb with noise.
mode_final = 32 -- Final wave mode to perturb with noise.

-- Derived physical quantities (using normalized code units).
vte = math.sqrt(T / mass_elc) -- Electron thermal velocity.
omega_pe = math.sqrt((charge_elc * charge_elc) * n0 / (epsilon0 * mass_elc)) -- Electron plasma frequency.
lambda_D = vte / omega_pe -- Electron Debye length.
lambda_J = lambda_D / math.sqrt(alpha_g) -- Jeans length.
omega_J = vte / lambda_J -- Jeans frequency.

-- Simulation parameters.
Nx = 64 -- Cell count (configuration space: x-direction).
Nvx = 128 -- Cell count (velocity space: vx-direction).
kx = 0.1 / lambda_J -- Perturbed wave number (x-direction); the fundamental mode of the box.
Lx = 2.0 * pi / kx -- Domain size (configuration space: x-direction).
vx_max = 64.0 * vte -- Domain boundary (velocity space: vx-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 10.0 / omega_J -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Electron beam number density: random-amplitude, random-phase noise in modes
-- mode_init-mode_final. The generator is re-seeded at every point so every point
-- sees the same noise.
function elcBeamDensity(x)
  math.randomseed(0)
  local perturb = 0.0
  for i = mode_init, mode_final do
    perturb = perturb + delta_n * math.random() * math.cos(i * kx * x + 2.0 * pi * math.random())
  end
  return 0.5 * (1.0 + perturb) * n0
end

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { -0.5 * Lx },
  upper = { 0.5 * Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Electrons: two counter-streaming beams with density noise.
  elc = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_elc, mass = mass_elc,

    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    numInit = 2,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]
          return elcBeamDensity(x) -- Electron beam number density.
        end,
        temperatureInit = function (t, xn)
          return T -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift -- Right-going beam drift velocity.
        end,

        correctAllMoments = true
      },
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]
          return elcBeamDensity(x) -- Electron beam number density.
        end,
        temperatureInit = function (t, xn)
          return T -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return -Vx_drift -- Left-going beam drift velocity.
        end,

        correctAllMoments = true
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Positrons: two unperturbed counter-streaming beams.
  pos = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_pos, mass = mass_pos,

    -- Velocity space grid.
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    numInit = 2,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          return 0.5 * n0 -- Positron beam number density.
        end,
        temperatureInit = function (t, xn)
          return T -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift -- Right-going beam drift velocity.
        end,

        correctAllMoments = true
      },
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          return 0.5 * n0 -- Positron beam number density.
        end,
        temperatureInit = function (t, xn)
          return T -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return -Vx_drift -- Left-going beam drift velocity.
        end,

        correctAllMoments = true
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  isElectrostatic = true,

  -- Field: electrostatics and self-gravity.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0,
    alphaG = alpha_g,

    poissonBcs = {
      lowerType = {
        G0.PoissonBc.bcPeriodic
      },
      upperType = {
        G0.PoissonBc.bcPeriodic
      }
    }
  }
}

vlasovApp:run()
