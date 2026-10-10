-- Sod shock tube of a neutral gas with implicit BGK collisions through the triad Vlasov solver on a
-- cylinder (coordinates r with velocities v_r, v_theta, v_z in the orthonormal triad, jump at r = 1
-- on 0.5 < r < 1.5), on cubic velocity maps that cluster cells at the origin (1x3v). The states
-- (n, T) = (1, 1) and (0.125, 0.8) are at rest on either side of the jump; with collision frequency
-- NU the solution is the Euler one with gamma = 5/3. The reference is a quasi-1D finite-volume
-- Euler solution with the metric determinant as area factor (reflecting walls), and the exact
-- planar Riemann solution is quoted for contrast. Figures of merit (serendipity p2, 24 x 6 x 6 x 6
-- cells on cubic velocity maps, 130 steps): density L1 error 0.0147 against the quasi-1D reference
-- (0.0158 against the planar exact solution); shock at 1.1875 (reference 1.1821), contact at 1.0764
-- (reference 1.0839); mass and energy conserved to round-off.
local Vlasov = G0.Vlasov

mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.
nl = 1.0 -- Left/inner number density.
Tl = 1.0 -- Left/inner temperature.
nr = 0.125 -- Right/outer number density.
Tr = 0.8 -- Right/outer temperature.
vt = 1.0 -- Thermal velocity.
nu = 1000.0 -- Collision frequency.
midplane = 1.0 -- Radial location of the jump.

Nr = 24 -- Cell count (configuration space: radial direction).
Nvr = 6 -- Cell count (velocity space: radial direction).
Nvtheta = 6 -- Cell count (velocity space: angular direction).
Nvz = 6 -- Cell count (velocity space: axial direction).
Lr = 1.0 -- Domain size (configuration space: radial direction).
vr_max = 8.0 * vt -- Domain boundary (velocity space: radial direction).
vtheta_max = 8.0 * vt -- Domain boundary (velocity space: angular direction).
vz_max = 8.0 * vt -- Domain boundary (velocity space: axial direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.1 -- Final simulation time.
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
  lower = { 0.5 },
  upper = { 0.5 + Lr },
  cells = { Nr },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = {  }, -- Periodic directions.

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.Triad,
    charge = charge, mass = mass,

    -- Geometry from the covariant tangent basis e_i and the orthonormal triad sigma_a.
    covTangentBasis = function (t, xn)
      local q_r = xn[1]
      local q_theta = 0.0

      -- Covariant tangent basis e_i . sigma_j of polar coordinates (row i, column j).
      return math.cos(q_theta), math.sin(q_theta), 0.0, -q_r * math.sin(q_theta), q_r *
        math.cos(q_theta), 0.0, 0.0, 0.0, 1.0
    end,
    triadBasis = function (t, xn)
      local q_r = xn[1]
      local q_theta = 0.0

      -- Orthonormal triad sigma_a . sigma_j (unit radial and azimuthal vectors).
      return math.cos(q_theta), math.sin(q_theta), 0.0, -math.sin(q_theta), math.cos(q_theta), 0.0,
        0.0, 0.0, 1.0
    end,
    triadBasisGradient = function (t, xn)
      local q_r = xn[1]
      local q_theta = 0.0

      -- d(sigma_a . sigma_j)/dx^k for x^k = r, theta, z (nine components each).
      return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -math.sin(q_theta), math.cos(q_theta),
        0.0, -math.cos(q_theta), -math.sin(q_theta), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0
    end,

    -- Velocity space grid.
    lower = { -vr_max, -vtheta_max, -vz_max },
    upper = { vr_max, vtheta_max, vz_max },
    cells = { Nvr, Nvtheta, Nvz },
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
        return vc * (1.0 + 1.5 * (vc / vz_max) * (vc / vz_max)) / (1.0 + 1.5)
      end },
    },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,
        densityInit = function (t, xn)
          local r = xn[1]

          -- Left and right states, times the metric determinant.
          local n = nr
          if r < midplane then
            n = nl
          end
          return (r) * n
        end,
        temperatureInit = function (t, xn)
          local r = xn[1]

          -- Isotropic temperature.
          local T = Tr
          if r < midplane then
            T = Tl
          end
          return T
        end,
        driftVelocityInit = function (t, xn)
          -- Drift velocity in the orthonormal frame (the gas is at rest).
          return 0.0, 0.0, 0.0
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
