local Vlasov = G0.Vlasov
-- Rectangular Riemann discontinuity test for GR Maxwell's equations. Tests flat space limit in 
-- Rectangular coordinate with M = 0, a = 0.
-- See: C3.2 From: S. S. (2004). Electrodynamics of black hole magnetospheres. Monthly Notices of the Royal Astronomical Society, 350(1), 407-426.
-- Panel A.


-- Simulation parameters.
Nx = 200 -- Cell count (x-direction).
Lx = 1.5 -- Domain size (x-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 1.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

midplane = 0.0 -- Midplane inital discontinuity in moments

massBH = 0.0 -- Mass of the black hole
spinBH = 0.0 -- Spin parameter, a = J/M, of the black hole ( Kerr-Schild, coordinates , |a| <= M).
B0 = 0.5

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
  lower = { -Lx },
  upper = { Lx },
  cells = { Nx },
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

    -- Speed factors
    epsilon0 = 1.0,
    mu0 = 1.0,
    elcErrorSpeedFactor = 1.0, -- chi = c*elcErrorSpeedFactor = 1.
    mgnErrorSpeedFactor = 1.0, -- gamma = c*mgnErrorSpeedFactor = 1.
    K_phi = 1.0, -- Damping Constant (electric field).
    K_psi = 1.0, -- Damping Constant (magnetic field).

    -- Use GR field ID
    fieldID = G0.FieldModel.GR,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      -- D^i, and B^i (contravaraint components)
      local Dx = 0.0 -- Total electric field (r-direction).
      local Dy = 0.0 -- Total electric field (y-direction).
      local Dz = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (r-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      if x < midplane then
        Dx = 0.0
        Dy = 0.0
        Dz = 0.0
        Bx = 1.0
        By = B0
        Bz = 0.0
      else
        Dx = 0.0
        Dy = 0.0
        Dz = 0.0
        Bx = 1.0
        By = -B0
        Bz = 0.0
      end

      -- Hand off the primative variables (U^\xi)
      return Dx, Dy, Dz, Bx, By, Bz, 0.0, 0.0
    end,

    bcx = { G0.FieldBc.bcCopy, G0.FieldBc.bcCopy } -- Copy boundary conditions (x-direction).
  }
}

-- Run application.
vlasovApp:run()
