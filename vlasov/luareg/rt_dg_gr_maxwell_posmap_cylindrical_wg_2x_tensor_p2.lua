-- Standing mode in a cylindrical (annular) waveguide with the general-relativistic Maxwell solver
-- (2x). The flat annulus geometry in (r, phi) with perfectly conducting walls at r = 2 and r = 5
-- carries the m = 4 TM mode D_z = R(r) cos(4 phi), R a combination of Bessel functions of order 4
-- vanishing at both walls, at frequency 2.4303. It is advanced for two periods, after which the
-- field must coincide with its initial state. A sinusoidal position map stretches the cells in r by
-- 40% with the endpoints fixed. After two periods the field returns to its initial state to 9e-6 in
-- the relative L2 norm, with the energy conserved to 2e-6.
local Vlasov = G0.Vlasov

pi = math.pi -- Pi.
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
omega = 2.430327042902498 -- Frequency of the mode (root of the cross product of Bessel functions).
t_period = 2.0 * math.pi / omega -- Period of the mode.
m_mode = 4 -- Azimuthal mode number.
kn = 0.0 -- Axial wave number.

Nr = 48 -- Cell count (r-direction).
Nphi = 192 -- Cell count (phi-direction).
r0 = 2.0 -- Inner radius.
r1 = 5.0 -- Outer radius.

-- Bessel functions of integer order from the C math library.
ffi = require "ffi"
ffi.cdef [[
  double jn(int, double);
  double yn(int, double);
]]

poly_order = 2 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 2.0 * t_period -- Final simulation time.
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
  upper = { r1, 2.0 * math.pi },
  cells = { Nr, Nphi },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Sinusoidal position map in r: cell size varies by +/- 40% with the endpoints fixed.
  mapc2pPos = {
    {
      pmap = function (t, xn)
        local rc = xn[1]
        local A = 0.4
        local L = r1 - r0

        return rc - (A * L / (2.0 * pi)) * math.sin(2.0 * pi * (rc - r0) / L)
      end
    },
  },

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = { 2 }, -- Periodic directions.

  -- Geometry.
  geom = Vlasov.Geom.new {
    massBH = 0.0, -- Black-hole mass.
    spinBH = 0.0, -- Black-hole spin parameter a = J / M (Kerr-Schild coordinates).
    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.Annulus
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
      local r, phi = xn[1], xn[2]

      -- TM mode of the annulus: D_z = R(r) cos(m phi) with R a combination of Bessel functions vanishing at r0 and r1.
      local a = 1.0
      local wkn = math.sqrt(omega * omega - kn * kn)
      local b = -a * ffi.C.jn(m_mode, r0 * wkn) / ffi.C.yn(m_mode, r0 * wkn)
      local Dz = (a * ffi.C.jn(m_mode, r * wkn) + b * ffi.C.yn(m_mode, r * wkn)) * math.cos(m_mode * phi)

      return 0.0, 0.0, Dz, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    bcx = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall },
  }
}

vlasovApp:run()
