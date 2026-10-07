-- Standing mode in a spherical shell with the general-relativistic Maxwell solver (2x). The flat
-- spherical (r, theta) geometry with perfectly conducting walls at r = 2 and r = 5 and pole
-- conditions at theta = 0 and pi carries the TE l = 2, m = 0 mode D_phi = 3 u(r) / r cos(theta), u
-- a combination of spherical Bessel functions of order 2 vanishing at both walls, at frequency
-- 2.2281. It is advanced for two periods, after which the field must coincide with its initial
-- state. The mode frequency is reproduced to 0.1% and after two periods the field returns to its
-- initial state to 8e-4 in the relative L2 norm, with the energy conserved to 3e-5.
local Vlasov = G0.Vlasov

pi = math.pi -- Pi.
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
omega = 2.2281321786 -- Frequency of the mode (root of the cross product of spherical Bessel functions).
t_period = 2.0 * math.pi / omega -- Period of the mode.

Nr = 24 -- Cell count (r-direction).
Ntheta = 48 -- Cell count (theta-direction).
r0 = 2.0 -- Inner radius.
r1 = 5.0 -- Outer radius.

-- Spherical Bessel functions of order 2.
local function sph_j2(x)
  return (3.0 / (x * x * x) - 1.0 / x) * math.sin(x) - 3.0 * math.cos(x) / (x * x)
end

local function sph_y2(x)
  return -(3.0 / (x * x * x) - 1.0 / x) * math.cos(x) - 3.0 * math.sin(x) / (x * x)
end
poly_order = 1 -- Polynomial order.
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

      -- TE l = 2, m = 0 mode of the spherical shell: D_phi = 3 u(r) / r cos(theta) with u a combination of
      -- spherical Bessel functions of order 2 vanishing at r0 and r1.
      local a = 1.0
      local b = -a * sph_j2(omega * r0) / sph_y2(omega * r0)
      local u = a * sph_j2(omega * r) + b * sph_y2(omega * r)
      local Dphi = 3.0 * u / r * math.cos(theta)

      return 0.0, 0.0, Dphi, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    bcx = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall },
    bcy = { G0.FieldBc.bcThetaPole, G0.FieldBc.bcThetaPole },
  }
}

vlasovApp:run()
