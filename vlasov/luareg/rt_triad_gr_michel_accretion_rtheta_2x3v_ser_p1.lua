-- Michel accretion of a Juttner gas onto a Schwarzschild black hole in ingoing Kerr-Schild
-- coordinates through the GR triad Vlasov solver with implicit BGK collisions (2x3v, r-theta preset
-- geometry with a single reflecting polar cell about the equator and cubic velocity maps). The
-- exact transonic Michel solution with sonic radius 8 M (temperature 0.037 at r = 15 M rising to
-- 0.15 inside the horizon, radial velocity -0.14 to -0.97) is tabulated and projected as the
-- initial state; the outer boundary holds the Michel state and the inner boundary at r = 1.8 M,
-- inside the horizon, absorbs the inflow. The driver runs the tensor p = 1 hybrid basis on a
-- coarser grid and shorter run (t = 5 M) than the radial driver and must reproduce it at equal
-- resolution. The accretion rate r^2 n U^r computed from the rest-frame LTE moments must equal the
-- exact -16 M^2 at every radius, and the profile must stay stationary over the run. Figures of
-- merit (serendipity p1, 10 x 1 x 6 x 6 x 6 cells on cubic velocity maps): accretion rate -16.001
-- at t = 0 (exact -16) and -16.40 at t = 5 M; the density, temperature and drift profiles change by
-- at most 3.3e-2, 5.0e-2 and 3.4e-2 inside r = 12 M.
local Vlasov = G0.Vlasov

mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.
mass_bh = 1.0 -- Black hole mass M.
nu = 15000.0 -- Collision frequency.
vt = 1.0 -- Momentum scale (units of m c).

Nr = 10 -- Cell count (configuration space: radial direction).
Ntheta = 1 -- Cell count (configuration space: polar direction).
Nvr = 6 -- Cell count (velocity space: radial momentum).
Nvtheta = 6 -- Cell count (velocity space: polar momentum).
Nvphi = 6 -- Cell count (velocity space: azimuthal momentum).
dtheta = 0.1 -- Polar extent of the single cell about the equator.
r_min = 1.8 -- Inner radius (inside the horizon; absorbing).
r_max = 15.0 -- Outer radius (fixed to the Michel state).
vr_max = 1.5 * vt -- Domain boundary (velocity space: radial momentum).
vtheta_max = 1.5 * vt -- Domain boundary (velocity space: polar momentum).
vphi_max = 1.5 * vt -- Domain boundary (velocity space: azimuthal momentum).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 5.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Michel accretion profile for a Juttner gas onto a Schwarzschild black hole (M = 1), sonic radius
-- 8 M: rest-frame density, temperature and the radial spatial four-velocity relative to the
-- Kerr-Schild normal observer, tabulated on 67 radii in [1.8, 15] (see the header).
michel_r = {
  1.800000, 2.000000, 2.200000, 2.400000, 2.600000, 2.800000,
  3.000000, 3.200000, 3.400000, 3.600000, 3.800000, 4.000000,
  4.200000, 4.400000, 4.600000, 4.800000, 5.000000, 5.200000,
  5.400000, 5.600000, 5.800000, 6.000000, 6.200000, 6.400000,
  6.600000, 6.800000, 7.000000, 7.200000, 7.400000, 7.600000,
  7.800000, 8.000000, 8.200000, 8.400000, 8.600000, 8.800000,
  9.000000, 9.200000, 9.400000, 9.600000, 9.800000, 10.000000,
  10.200000, 10.400000, 10.600000, 10.800000, 11.000000, 11.200000,
  11.400000, 11.600000, 11.800000, 12.000000, 12.200000, 12.400000,
  12.600000, 12.800000, 13.000000, 13.200000, 13.400000, 13.600000,
  13.800000, 14.000000, 14.200000, 14.400000, 14.600000, 14.800000,
  15.000000,
}
michel_n = {
  6.164246e+00, 5.373017e+00, 4.749860e+00, 4.247912e+00, 3.835976e+00, 3.492549e+00,
  3.202359e+00, 2.954288e+00, 2.740060e+00, 2.553398e+00, 2.389463e+00, 2.244466e+00,
  2.115403e+00, 1.999860e+00, 1.895884e+00, 1.801870e+00, 1.716494e+00, 1.638653e+00,
  1.567421e+00, 1.502014e+00, 1.441767e+00, 1.386110e+00, 1.334552e+00, 1.286670e+00,
  1.242095e+00, 1.200505e+00, 1.161620e+00, 1.125190e+00, 1.090997e+00, 1.058846e+00,
  1.028565e+00, 1.000000e+00, 9.730129e-01, 9.474801e-01, 9.232903e-01, 9.003430e-01,
  8.785475e-01, 8.578215e-01, 8.380903e-01, 8.192859e-01, 8.013460e-01, 7.842141e-01,
  7.678381e-01, 7.521704e-01, 7.371672e-01, 7.227883e-01, 7.089963e-01, 6.957572e-01,
  6.830390e-01, 6.708124e-01, 6.590503e-01, 6.477272e-01, 6.368197e-01, 6.263058e-01,
  6.161653e-01, 6.063789e-01, 5.969291e-01, 5.877990e-01, 5.789732e-01, 5.704371e-01,
  5.621770e-01, 5.541801e-01, 5.464342e-01, 5.389280e-01, 5.316509e-01, 5.245928e-01,
  5.177441e-01,
}
michel_T = {
  1.512292e-01, 1.409009e-01, 1.321331e-01, 1.245834e-01, 1.180048e-01, 1.122143e-01,
  1.070730e-01, 1.024735e-01, 9.833132e-02, 9.457910e-02, 9.116228e-02, 8.803625e-02,
  8.516412e-02, 8.251511e-02, 8.006333e-02, 7.778680e-02, 7.566677e-02, 7.368711e-02,
  7.183388e-02, 7.009497e-02, 6.845977e-02, 6.691900e-02, 6.546444e-02, 6.408886e-02,
  6.278578e-02, 6.154947e-02, 6.037475e-02, 5.925702e-02, 5.819210e-02, 5.717623e-02,
  5.620601e-02, 5.527834e-02, 5.439041e-02, 5.353964e-02, 5.272370e-02, 5.194042e-02,
  5.118784e-02, 5.046414e-02, 4.976765e-02, 4.909682e-02, 4.845022e-02, 4.782653e-02,
  4.722454e-02, 4.664308e-02, 4.608112e-02, 4.553766e-02, 4.501179e-02, 4.450264e-02,
  4.400940e-02, 4.353134e-02, 4.306774e-02, 4.261795e-02, 4.218134e-02, 4.175733e-02,
  4.134537e-02, 4.094494e-02, 4.055556e-02, 4.017677e-02, 3.980813e-02, 3.944923e-02,
  3.909968e-02, 3.875912e-02, 3.842720e-02, 3.810358e-02, 3.778796e-02, 3.748003e-02,
  3.717951e-02,
}
michel_u = {
  -5.141481e-02, -5.150117e-02, -5.133975e-02, -5.099059e-02, -5.049884e-02, -4.989886e-02,
  -4.921711e-02, -4.847410e-02, -4.768588e-02, -4.686505e-02, -4.602158e-02, -4.516340e-02,
  -4.429679e-02, -4.342677e-02, -4.255735e-02, -4.169170e-02, -4.083238e-02, -3.998138e-02,
  -3.914028e-02, -3.831033e-02, -3.749247e-02, -3.668743e-02, -3.589573e-02, -3.511776e-02,
  -3.435376e-02, -3.360388e-02, -3.286819e-02, -3.214667e-02, -3.143927e-02, -3.074587e-02,
  -3.006633e-02, -2.940047e-02, -2.874810e-02, -2.810899e-02, -2.748291e-02, -2.686962e-02,
  -2.626888e-02, -2.568042e-02, -2.510399e-02, -2.453932e-02, -2.398616e-02, -2.344425e-02,
  -2.291333e-02, -2.239314e-02, -2.188344e-02, -2.138398e-02, -2.089451e-02, -2.041480e-02,
  -1.994460e-02, -1.948371e-02, -1.903188e-02, -1.858890e-02, -1.815457e-02, -1.772866e-02,
  -1.731098e-02, -1.690134e-02, -1.649953e-02, -1.610537e-02, -1.571868e-02, -1.533928e-02,
  -1.496700e-02, -1.460167e-02, -1.424313e-02, -1.389122e-02, -1.354578e-02, -1.320666e-02,
  -1.287373e-02,
}

-- Linear interpolation of a tabulated profile in r.
function michel_interp(tab, r)
  local n = #michel_r
  if r <= michel_r[1] then
    return tab[1]
  end
  if r >= michel_r[n] then
    return tab[n]
  end
  local i = 1
  while michel_r[i + 1] < r do
    i = i + 1
  end
  local w = (r - michel_r[i]) / (michel_r[i + 1] - michel_r[i])
  return tab[i] * (1.0 - w) + tab[i + 1] * w
end

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { r_min, (math.pi / 2.0) - (dtheta / 2.0) },
  upper = { r_max, (math.pi / 2.0) + (dtheta / 2.0) },
  cells = { Nr, Ntheta },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = {  }, -- Periodic directions.

  -- Preset geometry.
  geom = Vlasov.Geom.new {
    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.GR_KS_rtheta,
    massBH = mass_bh,
    spinBH = 0.0,
  },

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.TriadGR,
    charge = charge, mass = mass,

    -- Velocity space grid.
    lower = { -vr_max, -vtheta_max, -vphi_max },
    upper = { vr_max, vtheta_max, vphi_max },
    cells = { Nvr, Nvtheta, Nvphi },
    -- Cubic velocity maps: linear core at the origin, stretched tails.
    mapc2pVel = {
      { vmap = function (t, xn)
        local vc = xn[1]
        return vc * (1.0 + 1.5 * (vc / vr_max) * (vc / vr_max)) / (1.0 + 1.5)
      end },
      { vmap = function (t, xn)
        local vc = xn[1]
        return vc * (1.0 + 1.5 * (vc / vtheta_max) * (vc / vtheta_max)) / (1.0 + 1.5)
      end },
      { vmap = function (t, xn)
        local vc = xn[1]
        return vc * (1.0 + 1.5 * (vc / vphi_max) * (vc / vphi_max)) / (1.0 + 1.5)
      end },
    },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,
        densityInit = function (t, xn)
          local r = xn[1]
          local theta = xn[2]

          -- Rest-frame density of the Michel profile times the metric determinant sqrt(h).
          return (r * math.sqrt((2.0 * mass_bh * r) + (r * r)) * math.sin(theta)) *
            michel_interp(michel_n, r)
        end,
        temperatureInit = function (t, xn)
          local r = xn[1]

          -- Rest-frame temperature of the Michel profile.
          return michel_interp(michel_T, r)
        end,
        driftVelocityInit = function (t, xn)
          local r = xn[1]

          -- Radial spatial four-velocity relative to the Kerr-Schild normal observer; no angular
          -- flow.
          return michel_interp(michel_u, r), 0.0, 0.0
        end,
        correctAllMoments = true,
        iterationEpsilon = 0.0,
        maxIterations = 0,
        useLastConverged = false
      }
    },

    collisions = {
      collisionID = G0.Collisions.BGK,
      selfNu = function (t, xn)
        return nu -- Collision frequency.
      end,
      useImplicitCollisionScheme = true
    },
    correct = {
      correctAllMoments = true,
      iterationEpsilon = 1e-12,
      maxIterations = 100,
      useLastConverged = false
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcAbsorb
      },
      upper = {
        type = G0.SpeciesBc.bcFixedFunc
      }
    },

    bcy = {
      lower = {
        type = G0.SpeciesBc.bcReflect
      },
      upper = {
        type = G0.SpeciesBc.bcReflect
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.LTEMoments, G0.Moment.EnergyMoment }
  },

  skipField = true,
}

vlasovApp:run()
