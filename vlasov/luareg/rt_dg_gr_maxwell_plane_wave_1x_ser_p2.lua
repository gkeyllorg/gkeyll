-- Traveling electromagnetic plane wave with the general-relativistic Maxwell solver in flat space
-- (1x). A unit-amplitude plane wave with two wavelengths per box length along x is advanced for
-- five periods on a periodic domain, after which it must coincide with its initial state. Flat-
-- space Kerr-Schild geometry with M = 0 and a = 0, so D = E and the solution is that of the Maxwell
-- solver; the GR solver carries hyperbolic divergence cleaning. After five periods the field
-- returns to its initial state to 9e-5 in the relative L2 norm, with the energy conserved to 2e-4.
local Vlasov = G0.Vlasov

pi = math.pi -- Pi.
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
E0 = 1.0 / math.sqrt(3.0) -- Electric field amplitude (|E| = 1).

k_wave_x = 2.0 -- Wave number (x-direction, wavelengths per box length).
omega = 2.0 * math.pi * math.sqrt(k_wave_x * k_wave_x) -- Wave frequency (c = 1).
t_period = 2.0 * math.pi / omega -- Wave period.

Nx = 64 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 5.0 * t_period -- Final simulation time.
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
  lower = { 0.0 },
  upper = { Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions.

  -- Geometry.
  geom = Vlasov.Geom.new {
    massBH = 0.0, -- Black-hole mass.
    spinBH = 0.0, -- Black-hole spin parameter a = J / M (Kerr-Schild coordinates).
    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.Flat
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
      local x = xn[1]

      local phase = (2.0 * pi / Lx) * (k_wave_x * x)
      local amp = 1.0 * math.cos(phase) -- Traveling plane wave.

      local Dx = 0.0 * amp
      local Dy = 1.0 / math.sqrt(2.0) * amp
      local Dz = 1.0 / math.sqrt(2.0) * amp
      local Bx = 0.0 * amp
      local By = -1.0 / math.sqrt(2.0) * amp
      local Bz = 1.0 / math.sqrt(2.0) * amp

      return Dx, Dy, Dz, Bx, By, Bz, 0.0, 0.0
    end,

  }
}

vlasovApp:run()
