local Vlasov = G0.Vlasov
ffi = require "ffi"
-- Cylindrical waveguide test for GR Maxwell's equations. Tests flat non-curved space limit in 
-- cylindrical coordinate with M = 0, a = 0.
-- See: https://ammar-hakim.org/sj/je/je36/je36-maxwell-cylindrical-waveguide.html

-- Mathematical constants (dimensionless).
pi = math.pi

-- Simulation parameters.
Nr = 48 -- Cell count (r-direction).
Nphi = 192 -- Cell count (phi-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

-- Waveguide parameters --
r0, r1 = 2.0, 5.0 -- inner and outer radii
m = 4 -- aziumuthal mode number
kn = 0.0 -- axial mode
w = 2.430327042902498 -- frequency of mode
tperiod = 2*math.pi/w -- period of mode

t_end = 2 * tperiod -- Final simulation time.
num_frames = 1 -- 30 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

massBH = 0.0 -- Mass of the black hole
spinBH = 0.0 -- Spin parameter, a = J/M, of the black hole ( Kerr-Schild, coordinates , |a| <= M).


-- Load the C math library bessel functions
ffi.cdef [[
  double jn(int, double);
  double yn(int, double);
]]

vlasovApp = Vlasov.App.new {

  -- Vlasov Geometry
  geom = Vlasov.Geom.new {
    
    -- Black hole parameters
    massBH = massBH,
    spinBH = spinBH,

    -- Specify the geometry rphi choice
    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.Annulus
    
  },

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { r0, 0.0 },
  upper = { r1, 2.0 * math.pi },
  cells = { Nr, Nphi },
  cflFrac = cfl_frac,

  mapc2pPos = {
    -- r mapping: sinusoidal stretching (cell size varies by +/- 40%, endpoints fixed).
    {
      pmap = function (t, xn)
        local rc = xn[1]
        local A = 0.4
        local L = r1 - r0
        return rc - (A * L / (2.0 * pi)) * math.sin(2.0 * pi * (rc - r0) / L)
      end
    },
  },
    
  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (r-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 2 }, -- Periodic directions (x- and phi-directions only).

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
      local r, phi = xn[1], xn[2]

      local a = 1.0
      local wkn = math.sqrt(w*w-kn*kn)
      local b = -a*ffi.C.jn(m,r0*wkn)/ffi.C.yn(m,r0*wkn)
      local Ez_r = a*ffi.C.jn(m,r*wkn) + b*ffi.C.yn(m,r*wkn)
      local Ez = Ez_r*math.cos(m*phi)

      -- D^i, and B^i (contravaraint components)
      local Dr = 0.0 -- Total electric field (r-direction).
      local Dphi = 0.0 -- Total electric field (phi-direction).
      local Dz = Ez -- Total electric field (z-direction). (Flat limit: E = D)

      local Br = 0.0 -- Total magnetic field (r-direction).
      local Bphi = 0.0 -- Total magnetic field (phi-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      -- Hand off the primative variables (U^\xi)
      return Dr, Dphi, Dz, Br, Bphi, Bz, 0.0, 0.0
    end,

    bcx = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall }, -- boundary conditions (r-direction).
  }
}

-- Run application.
vlasovApp:run()
