local Vlasov = G0.Vlasov
-- Plane-wave test for GR Maxwell's equations in the flat-space limit.
-- Uses rectangular coordinates with M = 0 and a = 0.

-- Derived physical quantities (using normalized code units).
J0 = 1.0 -- Applied current strength.

-- Simulation parameters.
Nx = 127 -- Cell count (x-direction).
Ny = 127 -- Cell count (y-direction).
Lx = 1.0 -- Domain size (x-direction).
Ly = 1.0 -- Domain size (y-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.
dx = Lx / Nx -- Cell size (x-direction).
dy = Ly / Ny -- Cell size (y-direction).

t_end = 2.0 -- Final simulation time.
num_frames = 30 -- 100 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

massBH = 0.0 -- Black-hole mass.
spinBH = 0.0 -- Black-hole spin parameter, a = J/M, in Kerr-Schild coordinates (|a| <= M).

vlasovApp = Vlasov.App.new {

  -- Geometry.
  geom = Vlasov.Geom.new {
    
    -- Black-hole parameters.
    massBH = massBH,
    spinBH = spinBH,

    -- Use the preset flat-space triad geometry.
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
  decompCuts = { 1, 1 }, -- Cuts in each coordinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- No periodic directions.

  -- Field.
  field = Vlasov.Field.new {

    -- Use GR Maxwell field model.
    fieldID = G0.FieldModel.GR,

    -- Field constants, divergence-cleaning speeds, and damping.
    epsilon0 = 1.0,
    mu0 = 1.0,
    elcErrorSpeedFactor = 1.0, -- chi = c*elcErrorSpeedFactor = 1.
    mgnErrorSpeedFactor = 1.0, -- gamma = c*mgnErrorSpeedFactor = 1.
    K_phi = 1.0, -- Damping constant (electric field).
    K_psi = 1.0, -- Damping constant (magnetic field).

    -- Initial conditions function.
    init = function (t, xn)
      local Ex = 0.0
      local Ey = 0.0
      local Ez = 0.0
    
      local Bx = 0.0
      local By = 0.0
      local Bz = 0.0

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    -- Applied current function.
    appliedCurrent = function (t, xn)
      local x, y = xn[1], xn[2]

      local app_x = 0.0 -- Applied current (x-direction).
      local app_y = 0.0 -- Applied current (y-direction).
      local app_z = 0.0 -- Applied current (z-direction).

      if (math.abs(x - Lx/2) < dx/2)  and (math.abs(y - Ly/2) < dy/2) then
        app_z = J0 -- Localized applied current (z-direction).
      else
        app_z = 0.0 -- No applied current outside the source cell.
      end

      return app_x, app_y, app_z
    end,
    evolveAppliedCurrent = false, -- Evolve applied current.
  }

  -- BCs
  bcx = { G0.FieldBc.bcCopy, G0.FieldBc.bcCopy }, -- Copy boundary conditions (x-direction).
  bcy = { G0.FieldBc.bcCopy, G0.FieldBc.bcCopy }, -- Copy boundary conditions (y-direction).
}

-- Run application.
vlasovApp:run()
