-- Axisymmetric spherical pulses with the general-relativistic Maxwell solver (2x). In the flat
-- spherical (r, theta) geometry between r = 2 and r = 5, two discontinuous shells of width 0.3
-- carry D_theta = A0 sin(theta) with B_phi = +A0 (outgoing shell at r = 3) and -A0
-- (incoming shell at r = 4). The shells propagate at the speed of light through open (copy) radial
-- boundaries and pole conditions at theta = 0 and pi. At t = 1.5 the leading edges of the two
-- shells have moved by 1.5 to within a cell
-- (the outgoing one to r = 4.65, the incoming one to r = 2.35), the incoming shell having steepened
-- as it converges; the energy changes by 0.5%.
local Vlasov = G0.Vlasov

pi = math.pi -- Pi.
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
A0 = 1.0 -- Pulse amplitude.
r0 = 2.0 -- Inner radius.
r1 = 5.0 -- Outer radius.
pulse_halfwidth = 0.15 -- Half width of each shell.
r_pulse1 = r0 + (r1 - r0) / 3.0 -- Center of the outgoing shell.
r_pulse2 = r0 + 2.0 * (r1 - r0) / 3.0 -- Center of the incoming shell.

Nr = 96 -- Cell count (r-direction).
Ntheta = 48 -- Cell count (theta-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 1.5 -- Final simulation time.
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
  lower = { r0, 0.0 },
  upper = { r1, math.pi },
  cells = { Nr, Ntheta },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = {  }, -- Periodic directions.

  -- Geometry.
  geom = Vlasov.Geom.new {
    massBH = 0.0, -- Black-hole mass.
    spinBH = 0.0, -- Black-hole spin parameter a = J / M (Kerr-Schild coordinates).
    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.Spherical_rtheta
  },

  -- Field (GR Maxwell, D and B).
  field = Vlasov.Field.new {
    fieldID = G0.FieldModel.GR,
    epsilon0 = epsilon0, mu0 = mu0,
    elcErrorSpeedFactor = 1.0, -- chi = c*elcErrorSpeedFactor = 1.
    mgnErrorSpeedFactor = 1.0, -- gamma = c*mgnErrorSpeedFactor = 1.
    K_phi = 1.0, -- Damping constant (electric field).
    K_psi = 1.0, -- Damping constant (magnetic field).

    -- Initial conditions function (D^i, B^i, phi, psi).
    init = function (t, xn)
      local r, theta = xn[1], xn[2]
      local hw = pulse_halfwidth

      -- Two discontinuous shells: S1 outgoing (D_theta and B_phi in phase), S2 incoming (opposite B_phi).
      local S1 = 0.0
      if r >= r_pulse1 - hw and r <= r_pulse1 + hw then
        S1 = 1.0
      end
      local S2 = 0.0
      if r >= r_pulse2 - hw and r <= r_pulse2 + hw then
        S2 = 1.0
      end
      local Dtheta = A0 * math.sin(theta) * (S1 + S2)
      local Bphi = A0 * (S1 - S2)

      return 0.0, Dtheta, 0.0, 0.0, 0.0, Bphi, 0.0, 0.0
    end,

    bcx = { G0.FieldBc.bcCopy, G0.FieldBc.bcCopy },
    bcy = { G0.FieldBc.bcThetaPole, G0.FieldBc.bcThetaPole },
  }
}

vlasovApp:run()
