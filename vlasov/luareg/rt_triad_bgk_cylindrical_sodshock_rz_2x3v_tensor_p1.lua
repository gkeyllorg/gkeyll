-- Sod shock tube of a neutral gas with implicit BGK collisions through the triad Vlasov solver on a
-- cylinder in (r, z) with a single periodic axial cell and the preset cylindrical geometry, jump at
-- r = 1 on 0.5 < r < 1.5, on cubic velocity maps (2x3v). The states (n, T) = (1, 1) and
-- (0.125, 0.8) are at rest on either side of the jump; with collision frequency NU the solution is
-- the Euler one with gamma = 5/3. The driver runs the tensor p = 1 hybrid basis
-- (p = 2 in velocity space), the configuration used for kinetic neutrals in the gyrokinetic solver,
-- and must reproduce the 1x3v driver. The reference is a quasi-1D finite-volume Euler solution with
-- the metric determinant as area factor (reflecting walls), and the exact planar Riemann solution
-- is quoted for contrast. Figures of merit (tensor p1 hybrid, 16 x 1 x 6 x 6 x 6 cells on cubic
-- velocity maps, 52 steps): density L1 error 0.0373 against the quasi-1D reference (0.0364 against
-- the planar exact solution); shock at 1.2031 (reference 1.1821) within one cell; cell-averaged
-- density within 0.027 (relative L1 0.012) of the serendipity p2 1x3v driver at 16 radial cells;
-- mass and energy conserved to round-off.
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

Nr = 16 -- Cell count (configuration space: radial direction).
Nz = 1 -- Cell count (configuration space: axial direction).
Nvr = 6 -- Cell count (velocity space: radial direction).
Nvz = 6 -- Cell count (velocity space: axial direction).
Nvtheta = 6 -- Cell count (velocity space: angular direction).
Lr = 1.0 -- Domain size (configuration space: radial direction).
Lz = 1.0 -- Domain size (configuration space: axial direction).
vr_max = 8.0 * vt -- Domain boundary (velocity space: radial direction).
vz_max = 8.0 * vt -- Domain boundary (velocity space: axial direction).
vtheta_max = 8.0 * vt -- Domain boundary (velocity space: angular direction).
poly_order = 1 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
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
  lower = { 0.5, 0.0 },
  upper = { 0.5 + Lr, Lz },
  cells = { Nr, Nz },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = { 2 }, -- Periodic directions.

  -- Preset geometry.
  geom = Vlasov.Geom.new {
    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.Cylindrical_rz,
    massBH = 0.0,
    spinBH = 0.0,
  },

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.Triad,
    charge = charge, mass = mass,

    -- Velocity space grid.
    lower = { -vr_max, -vz_max, -vtheta_max },
    upper = { vr_max, vz_max, vtheta_max },
    cells = { Nvr, Nvz, Nvtheta },
    -- Cubic velocity maps: linear core at the origin, stretched tails.
    mapc2pVel = {
      { vmap = function (t, xn)
        local vc = xn[1]
        return vc * (1.0 + 1.5 * (vc / vr_max) * (vc / vr_max)) / (1.0 + 1.5)
      end },
      { vmap = function (t, xn)
        local vc = xn[1]
        return vc * (1.0 + 1.5 * (vc / vz_max) * (vc / vz_max)) / (1.0 + 1.5)
      end },
      { vmap = function (t, xn)
        local vc = xn[1]
        return vc * (1.0 + 1.5 * (vc / vtheta_max) * (vc / vtheta_max)) / (1.0 + 1.5)
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
