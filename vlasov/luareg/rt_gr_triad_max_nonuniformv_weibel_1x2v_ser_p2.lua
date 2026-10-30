-- Relativistic Weibel instability of two cold counter-streaming electron beams through the GR triad
-- Vlasov solver and the GR-Maxwell field solver in the flat preset geometry (1x2v) on quadratic
-- velocity maps: the exact twin of the special-relativistic input rt_vlasov_sr_weibel_1x2v (beams
-- at +-0.9 c along y, T = 0.01 mc^2, filamentation mode k = 0.5 along x seeded by a magnetic field
-- perturbation, ions as a neutralizing background), so the general-relativistic current coupling
-- and field evolution must reduce to special relativity. Figures of merit (serendipity p2, 24 x 20
-- x 20 cells on quadratic velocity maps): magnetic-energy growth rate 0.5236 over t = 5-20 against
-- the kinetic theory value 0.525, saturation at t = 32.
local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

n0 = 1.0 -- Reference number density.
uy_drift = 0.9 -- Drift velocity of the beams (y-direction, units of c).
T_elc = 0.01 -- Electron temperature (units of mc^2).

alpha = 1e-05 -- Applied perturbation amplitude.
kx = 0.5 -- Perturbed wave number (x-direction).

-- Derived physical quantities (using normalized code units).
gamma_drift = 1.0 / math.sqrt(1.0 - (uy_drift * uy_drift)) -- Lorentz factor of the beams.
uy_drift_sr = gamma_drift * uy_drift -- Relativistic drift velocity of the beams (y-direction).

-- Simulation parameters.
Nx = 24 -- Cell count (configuration space: x-direction).
Nvx = 20 -- Cell count (velocity space: vx-direction).
Nvy = 20 -- Cell count (velocity space: vy-direction).
Lx = 2.0 * pi / kx -- Domain size (configuration space: x-direction).
vx_max = 6.0 -- Domain boundary (velocity space: vx-direction).
vy_max = 6.0 -- Domain boundary (velocity space: vy-direction).
vx_lin = 1.0 -- Velocity map: cell size at the origin relative to a uniform grid (vx-direction).
vy_lin = 1.0 -- Velocity map: cell size at the origin relative to a uniform grid (vy-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 35.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

vlasovApp = Vlasov.App.new {

  -- Flat preset geometry makes the GR triad/GR-Maxwell system reduce to the special-relativistic
  -- Weibel setup.
  geom = Vlasov.Geom.new {
    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.Flat,
    massBH = 0.0,
    spinBH = 0.0,
  },

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
    modelID = G0.Model.TriadGR,
    charge = charge_elc, mass = mass_elc,

    -- Velocity space grid.
    lower = { -1.0, -1.0 },
    upper = { 1.0, 1.0 },
    cells = { Nvx, Nvy },

    -- Quadratic velocity maps: finest cells at the origin, stretching to the domain boundary.
    mapc2pVel = {
      -- vx mapping
      {
        vmap = function (t, xn)
          local vc = xn[1]
          return vx_lin * vc + (vx_max - vx_lin) * vc * math.abs(vc)
        end
      },
      -- vy mapping
      {
        vmap = function (t, xn)
          local vc = xn[1]
          return vy_lin * vc + (vy_max - vy_lin) * vc * math.abs(vc)
        end
      }
    },

    -- Initial conditions.
    numInit = 2,
    projections = {
      -- Two counter-streaming Maxwell-Juttner distributions.
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          return 0.5 * n0 -- Beam number density.
        end,
        temperatureInit = function (t, xn)
          return T_elc -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return 0.0, uy_drift_sr -- Relativistic drift velocity of the first beam.
        end,

        correctAllMoments = true,
        useLastConverged = true
      },
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          return 0.5 * n0 -- Beam number density.
        end,
        temperatureInit = function (t, xn)
          return T_elc -- Isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return 0.0, -uy_drift_sr -- Relativistic drift velocity of the second beam.
        end,

        correctAllMoments = true,
        useLastConverged = true
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.LTEMoments }
  },

  -- Field.
  field = Vlasov.Field.new {
    -- GR-Maxwell field model in the flat limit.
    fieldID = G0.FieldModel.GR,
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = alpha * math.sin(kx * x) -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()
