local Vlasov = G0.Vlasov
ffi = require "ffi"
-- Spherical waveguide test for GR Maxwell's equations. Tests flat non-curved space limit in 
-- spherical coordinate with M = 0, a = 0.
-- Tests r-theta plane

-- Mathematical constants (dimensionless).
pi = math.pi

-- Simulation parameters.
Nr = 24 -- Cell count (r-direction).
Ntheta = 48 -- Cell count (theta-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

-- Waveguide parameters --
r0, r1 = 2.0, 5.0 -- inner and outer radii
w = 2.2281321786 -- frequency of mode
tperiod = 2*math.pi/w -- period of mode

t_end = 2 * tperiod -- Final simulation time.
num_frames = 30 -- 30 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

massBH = 0.0 -- Mass of the black hole
spinBH = 0.0 -- Spin parameter, a = J/M, of the black hole ( Kerr-Schild, coordinates , |a| <= M).


-- Local functions for spherical Bessel functions of order 2 (for TE, l=2, m=0 mode)
local function sph_j2(x)
  return (3.0/(x * x * x) - 1.0/x)*math.sin(x) - 3.0*math.cos(x)/(x * x)
end

local function sph_y2(x)
  return -(3.0/(x * x * x) - 1.0/x)*math.cos(x) - 3.0*math.sin(x)/(x * x)
end


vlasovApp = Vlasov.App.new {

  -- Vlasov Geometry
  geom = Vlasov.Geom.new {
    
    -- Black hole parameters
    massBH = massBH,
    spinBH = spinBH,

    -- Specify the geometry rtheta choice
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
  decompCuts = { 1 }, -- Cuts in each coodinate direction (r-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (x- and theta-directions only).

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
      local r, theta = xn[1], xn[2]

      -- TE, l=2, m=0 mode, spherical waveguide
      local E0 = 1.0
      local a = 1.0
      local b = -a * sph_j2(w*r0) / sph_y2(w*r0)
      local u = a*sph_j2(w*r) + b*sph_y2(w*r)

      -- D^i, and B^i (contravariant components)
      local Dr = 0.0 -- Total electric field (r-direction).
      local Dtheta = 0.0 -- Total electric field (theta-direction).
      local Dphi = 3.0 * E0 * u / r * math.cos(theta) -- Total electric field (phi-direction). (Flat limit: E = D)

      local Br = 0.0 -- Total magnetic field (r-direction).
      local Btheta = 0.0 -- Total magnetic field (theta-direction).
      local Bphi = 0.0 -- Total magnetic field (phi-direction).
      
      -- Hand off the primative varaibles (U^\xi)
      return Dr, Dtheta, Dphi, Br, Btheta, Bphi, 0.0, 0.0
    end,

    bcx = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall }, -- boundary conditions (r-direction).
    bcy = { G0.FieldBc.bcThetaPole, G0.FieldBc.bcThetaPole } -- boundary conditions (theta-direction).
  }
}

-- Run application.
vlasovApp:run()