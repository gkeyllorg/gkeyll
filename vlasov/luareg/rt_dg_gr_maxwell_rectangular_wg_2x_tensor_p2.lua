-- Standing modes in a rectangular waveguide with the general-relativistic Maxwell solver in flat
-- space (2x). The TM mode E_z = sin(pi x) sin(pi y) and the TE mode B_z = cos(pi x) cos(pi y) of a
-- 7 by 5 box with perfectly conducting walls, both at frequency sqrt(2) pi, are advanced for ten
-- periods, after which the field must coincide with its initial state. After ten periods the field
-- returns to its initial state to 3e-4 in the relative L2 norm, with the energy conserved to 5e-4.
local Vlasov = G0.Vlasov

pi = math.pi -- Pi.
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
omega = math.sqrt(2.0) * math.pi -- Frequency of the mode.
t_period = 2.0 * math.pi / omega -- Period of the mode.

Nx = 70 -- Cell count (x-direction).
Ny = 50 -- Cell count (y-direction).
Lx = 7.0 -- Domain size (x-direction).
Ly = 5.0 -- Domain size (y-direction).
poly_order = 2 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 10.0 * t_period -- Final simulation time.
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
  lower = { 0.0, 0.0 },
  upper = { Lx, Ly },
  cells = { Nx, Ny },
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
      local x, y = xn[1], xn[2]

      -- Standing modes of the box with k = (pi, pi): a TM mode in D_z and a TE mode in B_z.
      local Dz = math.sin(pi * x) * math.sin(pi * y)
      local Bz = math.cos(pi * x) * math.cos(pi * y)

      return 0.0, 0.0, Dz, 0.0, 0.0, Bz, 0.0, 0.0
    end,

    bcx = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall },
    bcy = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall },
  }
}

vlasovApp:run()
