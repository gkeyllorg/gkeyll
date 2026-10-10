-- Sod shock tube of a neutral gas with implicit BGK collisions through the triad Vlasov solver on
-- an annulus (polar coordinates r, theta with the orthonormal triad of unit radial and azimuthal
-- vectors, jump at r = 1 on 0.5 < r < 1.5) (1x2v). The states (n, T) = (1, 1) and (0.125, 0.8) are
-- at rest on either side of the jump; with collision frequency NU the solution is the Euler one
-- with gamma = 2. The Lua variants give the same geometry through the vierbein, the preset annulus
-- and a non-uniform radial mesh and must reproduce this driver. The reference is a quasi-1D
-- finite-volume Euler solution with the metric determinant as area factor (reflecting walls), and
-- the exact planar Riemann solution is quoted for contrast. Figures of merit (serendipity p2,
-- vierbein geometry input): identical to the tangent-basis input to round-off (density L1 error
-- 0.0077 against the quasi-1D reference).
local Vlasov = G0.Vlasov

mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.
nl = 1.0 -- Left/inner number density.
Tl = 1.0 -- Left/inner temperature.
nr = 0.125 -- Right/outer number density.
Tr = 0.8 -- Right/outer temperature.
vt = 1.0 -- Thermal velocity.
nu = 2000.0 -- Collision frequency.
midplane = 1.0 -- Radial location of the jump.

Nr = 64 -- Cell count (configuration space: radial direction).
Nvr = 12 -- Cell count (velocity space: radial direction).
Nvtheta = 12 -- Cell count (velocity space: angular direction).
Lr = 1.0 -- Domain size (configuration space: radial direction).
vr_max = 8.0 * vt -- Domain boundary (velocity space: radial direction).
vtheta_max = 8.0 * vt -- Domain boundary (velocity space: angular direction).
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

    -- Geometry from the vierbein e_i^a and its gradient.
    useVierbein = true,
    vierbein = function (t, xn)
      local q_r = xn[1]

      -- Vierbein e_i^a of polar coordinates: e_r = sigma_r, e_theta = r sigma_theta.
      return 1.0, 0.0, 0.0, q_r
    end,
    vierbeinGradient = function (t, xn)
      -- d(e_i^a)/dr (four components) and d(e_i^a)/dtheta (four components).
      return 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0
    end,

    -- Velocity space grid.
    lower = { -vr_max, -vtheta_max },
    upper = { vr_max, vtheta_max },
    cells = { Nvr, Nvtheta },

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
          return 0.0, 0.0
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
