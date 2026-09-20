local Vlasov = G0.Vlasov
ffi = require "ffi"

-- Axisymmetric spherical pulse test for GR Maxwell's equations.
-- Flat-space limit in spherical (r,theta), with discontinuous radial shells.
-- One pulse is outward-moving, one is inward-moving.

-- Mathematical constants (dimensionless).
pi = math.pi

-- Simulation parameters.
Nr = 96        -- Cell count (r-direction).
Ntheta = 48    -- Cell count (theta-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

-- Domain parameters.
r0, r1 = 2.0, 5.0 -- inner and outer radii

-- Final time/output.
t_end = 1.5
num_frames = 1 -- 30
field_energy_calcs = GKYL_MAX_INT
integrated_mom_calcs = GKYL_MAX_INT
integrated_L2_f_calcs = GKYL_MAX_INT
dt_failure_tol = 1.0e-4
num_failures_max = 20

massBH = 0.0
spinBH = 0.0

-- Pulse parameters.
A0 = 1.0
domain_len = r1 - r0

-- Pulse centers at 1/3 and 2/3 of the domain.
r_pulse1_center = r0 + domain_len / 3.0
r_pulse2_center = r0 + 2.0 * domain_len / 3.0

-- Half-width of each discontinuous shell.
pulse_halfwidth = 0.2

-- Pulse 1: outward-moving shell
r_pulse1_left  = r_pulse1_center - pulse_halfwidth
r_pulse1_right = r_pulse1_center + pulse_halfwidth

-- Pulse 2: inward-moving shell
r_pulse2_left  = r_pulse2_center - pulse_halfwidth
r_pulse2_right = r_pulse2_center + pulse_halfwidth

vlasovApp = Vlasov.App.new {

  -- Geometry.
  geom = Vlasov.Geom.new {
    massBH = massBH,
    spinBH = spinBH,

    -- Flat spherical r-theta geometry.
    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.Spherical_rtheta
  },

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,

  lower = { r0, 0.0 },
  upper = { r1, math.pi },
  cells = { Nr, Ntheta },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 },

  -- Boundary conditions for configuration space.
  periodicDirs = { },

  -- Field.
  field = Vlasov.Field.new {
    fieldID = G0.FieldModel.GR,

    -- Speed factors
    epsilon0 = 1.0,
    mu0 = 1.0,
    elcErrorSpeedFactor = 1.0, -- chi = c*elcErrorSpeedFactor = 1.
    mgnErrorSpeedFactor = 1.0, -- gamma = c*mgnErrorSpeedFactor = 1.
    K_phi = 1.0, -- Damping Constant (electric field).
    K_psi = 1.0, -- Damping Constant (magnetic field).

    -- Initial conditions.
    init = function (t, xn)
      local r, theta = xn[1], xn[2]

      -- Discontinuous shell indicators.
      local S1 = 0.0
      if r >= r_pulse1_left and r <= r_pulse1_right then
        S1 = 1.0
      end

      local S2 = 0.0
      if r >= r_pulse2_left and r <= r_pulse2_right then
        S2 = 1.0
      end

      -- Pulse 1: outward
      local Dtheta1 = A0 * math.sin(theta) * S1
      local Bphi1   = A0 * S1

      -- Pulse 2: inward
      local Dtheta2 = A0 * math.sin(theta) * S2
      local Bphi2   = - A0 * S2

      -- Total contravariant components.
      local Dr     = 0.0
      local Dtheta = Dtheta1 + Dtheta2
      local Dphi   = 0.0

      local Br     = 0.0
      local Btheta = 0.0
      local Bphi   = Bphi1 + Bphi2

      -- Hand off the primative variables (U^\xi)
      return Dr, Dtheta, Dphi, Br, Btheta, Bphi, 0.0, 0.0
    end,

    -- Radial boundaries.
    bcx = { G0.FieldBc.bcCopy, G0.FieldBc.bcCopy },

    -- Theta-pole boundaries.
    bcy = { G0.FieldBc.bcThetaPole, G0.FieldBc.bcThetaPole }
  }
}

-- Run application.
vlasovApp:run()