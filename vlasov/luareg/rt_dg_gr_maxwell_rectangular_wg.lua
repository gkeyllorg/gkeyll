local Vlasov = G0.Vlasov
ffi = require "ffi"
-- Rectangular waveguide test for GR Maxwell's equations. Tests flat space limit in 
-- Rectangular coordinate with M = 0, a = 0.
-- See: 

-- Mathematical constants (dimensionless).
pi = math.pi

-- Derived physical quantities (using normalized code units).
omega = math.sqrt(2.0) * pi -- Frequency of mode.
t_period = 2.0 * pi / omega -- Period of mode.

-- Simulation parameters.
Nx = 70 -- Cell count (x-direction).
Ny = 50 -- Cell count (y-direction).
Lx = 7.0 -- Domain size (x-direction).
Ly = 5.0 -- Domain size (y-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 10.0 * t_period -- Final simulation time.
num_frames = 1 -- 100 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

massBH = 0.0 -- Mass of the black hole
spinBH = 0.0 -- Spin parameter, a = J/M, of the black hole ( Kerr-Schild, coordinates , |a| <= M).

vlasovApp = Vlasov.App.new {

  -- Vlasov Geometry
  geom = Vlasov.Geom.new {
    
    -- Black hole parameters
    massBH = massBH,
    spinBH = spinBH,

    -- Specify the geometry ry choice
    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.Flat
    
  },

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
  decompCuts = { 1 }, -- Cuts in each coodinate direction (r-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (x- and y-directions only).

  -- Field.
  field = Vlasov.Field.new {

    -- Use GR field ID
    fieldID = G0.FieldModel.GR,

    -- Speed factors
    epsilon0 = 1.0,
    mu0 = 1.0,
    elcErrorSpeedFactor = 1.0, -- chi = c*elcErrorSpeedFactor = 1.
    mgnErrorSpeedFactor = 1.0, -- gamma = c*mgnErrorSpeedFactor = 1.
    K_phi = 1.0, -- Damping Constant (electric field).
    K_psi = 1.0, -- Damping Constant (magnetic field).

    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      -- D^i, and B^i (contravaraint components)
      local Dr = 0.0 -- Total electric field (r-direction).
      local Dy = 0.0 -- Total electric field (y-direction).
      local Dz = math.sin(pi * x) * math.sin(pi * y) -- Total electric field (z-direction). (Flat limit: E = D)

      local Br = 0.0 -- Total magnetic field (r-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = math.cos(pi * x) * math.cos(pi * y) -- Total magnetic field (z-direction).

      -- Hand off the primative variables (U^\xi)
      return Dr, Dy, Dz, Br, By, Bz, 0.0, 0.0
    end,

    bcx = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall }, -- PEC wall boundary conditions (x-direction).
    bcy = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall } -- PEC wall boundary conditions (y-direction).  }
  }
}

-- Run application.
vlasovApp:run()
