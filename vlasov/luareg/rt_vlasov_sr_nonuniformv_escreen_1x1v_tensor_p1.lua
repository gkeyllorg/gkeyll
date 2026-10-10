-- Electric-field screening by pair creation (special-relativistic Vlasov-Maxwell, 1x1v): the
-- setup of Tolman, Philippov and Timokhin (2022, ApJL 933, L37) as adapted in the Gkeyll
-- relativistic Vlasov paper. An electron-positron plasma (n = 1 per species, T = mc^2) sits in a
-- uniform electric field E0 = 50 while cold pairs (T = 0.1 mc^2) drifting at u = 5 are injected at
-- 0.5 per species per unit time. The injected pairs screen the field, which then oscillates at the
-- relativistic plasma frequency with an amplitude damped by the fresh pairs spun up at every zero
-- crossing. The momentum grid is linear (cell 0.4) near the origin and exponential out to p = 400
-- so the particles accelerated during screening (p ~ 2 sqrt(xi)/3 = 330, xi = E0^3/S = 2.5e5)
-- stay on the grid. Reference: a Lagrangian k = 0 solution of the same equations (scratch
-- kzero.py). Figures of merit (tensor p1, 16 x 96 cells): screening at t = 8.42 (reference
-- 8.42, cold-beam estimate (sqrt(n0^2 + S E0) - n0)/S = 8.2), peak field 20.2, 12.4 and 9.4 at
-- t = 12.2, 19.5 and 24.8 (reference 20.1, 12.1 and 9.1), envelope decay rate 0.0613 over the
-- first three peaks against 0.0639 for the reference and 0.0624 for the PIC fit of Tolman et al.
-- (a = 0.43 xi^-0.47 per t0 = mc/(e E0)); the rms field tracks the reference to 2.3% of E0
-- (0.6% in L1); the injected mass is exact.
local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_pos = 1.0 -- Positron mass.
charge_pos = 1.0 -- Positron charge.

n0 = 1.0 -- Initial number density of each species.
T0 = 1.0 -- Initial temperature (units of mc^2).
E0 = 50.0 -- Initial electric field.
n_src = 0.5 -- Pair injection rate (number density per unit time, each species).
T_src = 0.1 -- Temperature of the injected pairs (units of mc^2).
u_src = 5.0 -- Drift (four-) velocity of the injected pairs (x-direction).

perturb = 1.0e-3 -- Relative amplitude of the seeded density modes.
num_modes = 4 -- Number of seeded density modes.

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Npx = 96 -- Cell count (momentum space: px-direction).
Lx = 1.0 -- Domain size (configuration space: x-direction).
px_max = 400.0 -- Domain boundary (momentum space: px-direction).
px_lin = 0.14 -- Momentum map: cell size of the linear part of the map at the origin.
poly_order = 1 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 30.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Seeded density modes: sum_j (perturb / j) cos(j k x + 2 pi j 0.3), j = 1 .. num_modes.
local function density_modes(x)
  local kx = 2.0 * pi / Lx
  local sum = 0.0
  for j = 1, num_modes do
    sum = sum + (perturb / j) * math.cos((j * kx * x) + (2.0 * pi * j * 0.3))
  end
  return sum
end

-- The electric field that balances the seeded charge density through Gauss's law.
local function field_modes(x)
  local kx = 2.0 * pi / Lx
  local sum = 0.0
  for j = 1, num_modes do
    sum = sum + (perturb / j) * math.sin((j * kx * x) + (2.0 * pi * j * 0.3)) / (j * kx)
  end
  return 2.0 * n0 * sum
end

-- Linear-to-exponential momentum map: cells of size px_lin at the origin, growing
-- exponentially to the domain boundary px_max.
local function mapc2p_px(vc)
  if vc < 0.0 then
    return (px_lin * vc * Npx) - math.exp(-vc * math.log(px_max)) + 1.0
  else
    return (px_lin * vc * Npx) + math.exp(vc * math.log(px_max)) - 1.0
  end
end

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
    modelID = G0.Model.SR,
    charge = charge_elc, mass = mass_elc,

    -- Momentum space grid.
    lower = { -1.0 },
    upper = { 1.0 },
    cells = { Npx },

    mapc2pVel = {
      -- px mapping
      {
        vmap = function (t, xn)
          return mapc2p_px(xn[1])
        end
      },
    },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]
          return n0 * (1.0 - density_modes(x)) -- Electron density, seeded modes removed.
        end,
        temperatureInit = function (t, xn)
          return T0 -- Initial temperature.
        end,
        driftVelocityInit = function (t, xn)
          return 0.0 -- The initial plasma is at rest.
        end,

        correctAllMoments = true,
        useLastConverged = true
      },
    },

    -- Pair injection.
    source = {
      sourceID = G0.Source.Proj,

      numSources = 1,
      projections = {
        {
          projectionID = G0.Projection.LTE,

          densityInit = function (t, xn)
            -- Pair injection rate. The LTE projection takes the rest-frame density, so divide
            -- by the Lorentz factor of the drift to inject n_src pairs per unit time (lab frame).
            return n_src / math.sqrt(1.0 + (u_src * u_src))
          end,
          temperatureInit = function (t, xn)
            return T_src -- Temperature of the injected pairs.
          end,
          driftVelocityInit = function (t, xn)
            return u_src -- Drift (four-) velocity of the injected pairs.
          end,

          correctAllMoments = true,
          useLastConverged = true
        },
      },
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.LTEMoments }
  },

  -- Positrons.
  pos = Vlasov.Species.new {
    modelID = G0.Model.SR,
    charge = charge_pos, mass = mass_pos,

    -- Momentum space grid.
    lower = { -1.0 },
    upper = { 1.0 },
    cells = { Npx },

    mapc2pVel = {
      -- px mapping
      {
        vmap = function (t, xn)
          return mapc2p_px(xn[1])
        end
      },
    },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]
          return n0 * (1.0 + density_modes(x)) -- Positron density, seeded modes added.
        end,
        temperatureInit = function (t, xn)
          return T0 -- Initial temperature.
        end,
        driftVelocityInit = function (t, xn)
          return 0.0 -- The initial plasma is at rest.
        end,

        correctAllMoments = true,
        useLastConverged = true
      },
    },

    -- Pair injection.
    source = {
      sourceID = G0.Source.Proj,

      numSources = 1,
      projections = {
        {
          projectionID = G0.Projection.LTE,

          densityInit = function (t, xn)
            -- Pair injection rate. The LTE projection takes the rest-frame density, so divide
            -- by the Lorentz factor of the drift to inject n_src pairs per unit time (lab frame).
            return n_src / math.sqrt(1.0 + (u_src * u_src))
          end,
          temperatureInit = function (t, xn)
            return T_src -- Temperature of the injected pairs.
          end,
          driftVelocityInit = function (t, xn)
            return u_src -- Drift (four-) velocity of the injected pairs.
          end,

          correctAllMoments = true,
          useLastConverged = true
        },
      },
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.LTEMoments }
  },

  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local Ex = E0 + field_modes(x) -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()
