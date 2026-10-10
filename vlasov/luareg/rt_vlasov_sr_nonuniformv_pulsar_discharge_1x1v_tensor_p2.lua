-- Pulsar polar-cap discharge with an adaptive pair source (special-relativistic Vlasov-Maxwell,
-- 1x1v), a scaled 1D version of the setup of Chernoglazov, Philippov and Timokhin (2024, ApJL
-- 974, L32) without photons: a positron atmosphere 100 exp(-x) held by gravity (acceleration -T
-- where the atmosphere density exceeds 0.5) over a magnetospheric pair plasma, electrons with the
-- Goldreich-Julian charge density -(1 + 0.8 x/L) and a flow that carries the applied
-- super-Goldreich-Julian current J0 = 2 on top of the positron flux, so that E = 0 satisfies
-- Gauss's law for the perturbation field and the plasma carries J0 at t = 0. Fixed-function
-- boundaries keep the star-side and magnetospheric states in the ghost cells. As the atmosphere
-- settles the electron flux can no longer carry J0 and a gap opens above the atmosphere,
-- dE/dt = -(j + J0); once particles pass p > 25 (a proxy for curvature-photon emission) each
-- species injects pairs at 0.5 (n_e + n_p above the threshold) per unit time, those bred by
-- particles moving away from the star drifting at p = 2 and those bred by particles moving toward
-- the star at p = -2, with T = 0.2 mc^2 and the rescaled densities smoothed by ten passes of the
-- Gaussian filter; the cascade screens the gap into a sustained discharge. The momentum grid is
-- linear (cell 0.3) at the origin and exponential out to p = 200.
-- Figures of merit (tensor p2, 256 x 192 cells): the gap field above the atmosphere reaches -0.73,
-- -1.05 and -1.38 at t = 10, 15 and 20 (x = 6 to 9; field energy 1.51, 5.90 and 17.1), the cascade
-- starts at t = 29.8, 370 pairs per unit area are created by t = 60, when the discharge injects
-- 24.0 pairs per unit time with the electron density above the threshold at 1.83 and the gap field
-- at -2.29 (x = 10).

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

T = 0.2 -- Temperature of the plasma and of the injected pairs (units of mc^2).
J0 = 2.0 -- Applied current (x-direction).
grav = 1.0 -- Strength of gravity (inverse scale height of the atmosphere).
n_atm = 100.0 -- Atmosphere density at the star.
n_mag = 5.0 -- Positron density of the magnetospheric plasma.
v_mag = 0.8 -- Flow velocity of the magnetospheric positrons (units of c).
rate_src = 0.5 -- Pair injection rate per unit density above the threshold (each species).
p_thresh = 25.0 -- Momentum above which particles breed pairs.
f_thresh = 1.0e-3 -- Distribution function below which the threshold moment is not accumulated.
u_src = 2.0 -- Drift (four-) velocity of the injected pairs (x-direction).
num_filters = 10 -- Number of Gaussian filter passes on the rescaled source density.

-- Simulation parameters.
Nx = 256 -- Cell count (configuration space: x-direction).
Npx = 192 -- Cell count (momentum space: px-direction).
Lx = 50.0 -- Domain size (configuration space: x-direction).
px_max = 200.0 -- Domain boundary (momentum space: px-direction).
px_lin = 0.14 -- Momentum map: cell size of the linear part of the map at the origin.
poly_order = 2 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 60.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Lab-frame positron density: the atmosphere over the magnetospheric plasma.
local function pos_density(x)
  return (n_atm * math.exp(-grav * x)) + n_mag
end

-- Lab-frame electron density: the positron density plus the Goldreich-Julian charge density
-- -(1 + 0.8 x/L), so that Gauss's law holds for the perturbation field E = 0 at t = 0.
local function elc_density(x)
  return (1.0 + ((0.8 * x) / Lx)) + pos_density(x)
end

-- Gravity acts where the atmosphere is denser than 0.5.
local function applied_accel(x)
  local n_atm_x = n_atm * math.exp(-grav * x)
  local accel = 0.0
  if n_atm_x > 0.5 then
    accel = -grav * T
  end
  return accel, 0.0, 0.0
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
  periodicDirs = { }, -- Periodic directions (none).

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
          -- The electron flux carries the applied current on top of the positron flux.
          local n = elc_density(x)
          local beta = (J0 + (n_mag * v_mag)) / n
          return n * math.sqrt(1.0 - (beta * beta)) -- Rest-frame electron density.
        end,
        temperatureInit = function (t, xn)
          return T -- Temperature.
        end,
        driftVelocityInit = function (t, xn)
          local x = xn[1]
          local beta = (J0 + (n_mag * v_mag)) / elc_density(x)
          return beta / math.sqrt(1.0 - (beta * beta)) -- Electron drift (four-) velocity.
        end,

        correctAllMoments = true,
        useLastConverged = true
      },
    },

    -- Adaptive pair injection, rescaled by the density above the threshold in both species.
    source = {
      sourceID = G0.Source.Adapt,
      -- Particles moving away from the star breed outgoing pairs, particles moving toward
      -- the star breed pairs that bombard the surface.
      numCrossSources = 2,
      sourceWith = { "pos", "pos" },
      sourceWithVThresh = { p_thresh, p_thresh },
      sourceWithFThresh = { f_thresh, f_thresh },
      sourceWithUpperHalf = { true, false },
      sourceWithProj = { 1, 2 },
      writeSource = true,
      filter = true,
      numFilters = num_filters,

      numSources = 2,
      projections = {
        {
          projectionID = G0.Projection.LTE,

          densityInit = function (t, xn)
            -- Pair injection rate per unit density above the threshold. The LTE projection
            -- takes the rest-frame density, so divide by the Lorentz factor of the drift to
            -- inject rate_src pairs per unit time in the lab frame.
            return rate_src / math.sqrt(1.0 + (u_src * u_src))
          end,
          temperatureInit = function (t, xn)
            return T -- Temperature of the injected pairs.
          end,
          driftVelocityInit = function (t, xn)
            return u_src -- Drift (four-) velocity of the pairs moving away from the star.
          end,

          correctAllMoments = true,
          useLastConverged = true
        },
        {
          projectionID = G0.Projection.LTE,

          densityInit = function (t, xn)
            return rate_src / math.sqrt(1.0 + (u_src * u_src)) -- As above.
          end,
          temperatureInit = function (t, xn)
            return T -- Temperature of the injected pairs.
          end,
          driftVelocityInit = function (t, xn)
            return -u_src -- Drift (four-) velocity of the pairs moving toward the star.
          end,

          correctAllMoments = true,
          useLastConverged = true
        },
      },
    },

    appliedAccel = function (t, xn)
      return applied_accel(xn[1])
    end,

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcFixedFunc
      },
      upper = {
        type = G0.SpeciesBc.bcFixedFunc
      },
    },

    evolve = true, -- Evolve species?
    -- The LTE moments are left out: the source diagnostics would divide by the zero source
    -- before any particle has passed the threshold.
    diagnostics = { G0.Moment.M0, G0.Moment.M1 }
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
          -- The magnetospheric positron flux is shared with the atmosphere.
          local n = pos_density(x)
          local beta = (n_mag * v_mag) / n
          return n * math.sqrt(1.0 - (beta * beta)) -- Rest-frame positron density.
        end,
        temperatureInit = function (t, xn)
          return T -- Temperature.
        end,
        driftVelocityInit = function (t, xn)
          local x = xn[1]
          local beta = (n_mag * v_mag) / pos_density(x)
          return beta / math.sqrt(1.0 - (beta * beta)) -- Positron drift (four-) velocity.
        end,

        correctAllMoments = true,
        useLastConverged = true
      },
    },

    -- Adaptive pair injection, rescaled by the density above the threshold in both species.
    source = {
      sourceID = G0.Source.Adapt,
      -- Particles moving away from the star breed outgoing pairs, particles moving toward
      -- the star breed pairs that bombard the surface.
      numCrossSources = 2,
      sourceWith = { "elc", "elc" },
      sourceWithVThresh = { p_thresh, p_thresh },
      sourceWithFThresh = { f_thresh, f_thresh },
      sourceWithUpperHalf = { true, false },
      sourceWithProj = { 1, 2 },
      writeSource = true,
      filter = true,
      numFilters = num_filters,

      numSources = 2,
      projections = {
        {
          projectionID = G0.Projection.LTE,

          densityInit = function (t, xn)
            -- Pair injection rate per unit density above the threshold. The LTE projection
            -- takes the rest-frame density, so divide by the Lorentz factor of the drift to
            -- inject rate_src pairs per unit time in the lab frame.
            return rate_src / math.sqrt(1.0 + (u_src * u_src))
          end,
          temperatureInit = function (t, xn)
            return T -- Temperature of the injected pairs.
          end,
          driftVelocityInit = function (t, xn)
            return u_src -- Drift (four-) velocity of the pairs moving away from the star.
          end,

          correctAllMoments = true,
          useLastConverged = true
        },
        {
          projectionID = G0.Projection.LTE,

          densityInit = function (t, xn)
            return rate_src / math.sqrt(1.0 + (u_src * u_src)) -- As above.
          end,
          temperatureInit = function (t, xn)
            return T -- Temperature of the injected pairs.
          end,
          driftVelocityInit = function (t, xn)
            return -u_src -- Drift (four-) velocity of the pairs moving toward the star.
          end,

          correctAllMoments = true,
          useLastConverged = true
        },
      },
    },

    appliedAccel = function (t, xn)
      return applied_accel(xn[1])
    end,

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcFixedFunc
      },
      upper = {
        type = G0.SpeciesBc.bcFixedFunc
      },
    },

    evolve = true, -- Evolve species?
    -- The LTE moments are left out: the source diagnostics would divide by the zero source
    -- before any particle has passed the threshold.
    diagnostics = { G0.Moment.M0, G0.Moment.M1 }
  },

  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    appliedCurrent = function (t, xn)
      return J0, 0.0, 0.0 -- Applied current (x-direction).
    end,

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()
