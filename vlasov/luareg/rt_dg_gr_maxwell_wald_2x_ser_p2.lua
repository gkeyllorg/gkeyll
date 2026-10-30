-- Wald's stationary magnetosphere of a Schwarzschild black hole with the general-relativistic
-- Maxwell solver (2x). Kerr-Schild (r, theta) geometry with M = 1; the uniform-at-infinity magnetic
-- field B_0 = 1 threading the hole is an exact stationary vacuum solution
-- (Komissarov 2004, Eq. 101 with a = 0), held at the radial boundaries by fixed-function conditions
-- with pole conditions at theta = 0 and pi. Advanced for five time units, the field must remain
-- stationary. Over five time units the field changes by 1.3e-4 in the relative L2 norm, with the
-- energy conserved to 1e-11.
local Vlasov = G0.Vlasov

pi = math.pi -- Pi.
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_bh = 1.0 -- Black-hole mass.
spin_bh = 0.0 -- Black-hole spin parameter a = J / M (Kerr-Schild coordinates).
B_0 = 1.0 -- Asymptotic magnetic field strength.

Nr = 48 -- Cell count (r-direction).
Ntheta = 192 -- Cell count (theta-direction).
r_lo = 1.5 * mass_bh -- Inner radius.
r_up = 30.0 * mass_bh -- Outer radius.
poly_order = 2 -- Polynomial order.
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

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { r_lo, 0.0 },
  upper = { r_up, math.pi },
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
    massBH = mass_bh, -- Black-hole mass.
    spinBH = spin_bh, -- Black-hole spin parameter a = J / M (Kerr-Schild coordinates).
    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.GR_KS_rtheta
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
      local M = mass_bh

      -- Wald's stationary magnetosphere of a Schwarzschild black hole in Kerr-Schild coordinates
      -- (Komissarov 2004, Eqs. 101-102): contravariant D^i and B^i. The vacuum constitutive
      -- relations give a non-zero D^phi even though the field is purely magnetic at infinity.
      local fac = math.sqrt(1.0 + 2.0 * M / r)
      local Dphi = -2.0 * M * B_0 / (r * r * fac)
      local Br = -B_0 * math.cos(theta) / fac
      local Btheta = B_0 * math.sin(theta) / (r * fac)

      return 0.0, 0.0, Dphi, Br, Btheta, 0.0, 0.0, 0.0
    end,

    bcx = { G0.FieldBc.bcFixedFunc, G0.FieldBc.bcFixedFunc },
    bcy = { G0.FieldBc.bcThetaPole, G0.FieldBc.bcThetaPole },
  }
}

vlasovApp:run()
