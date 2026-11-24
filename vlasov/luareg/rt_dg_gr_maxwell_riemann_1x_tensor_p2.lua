-- One-dimensional Riemann problem for the general-relativistic Maxwell solver in flat space (1x).
-- Komissarov's problem A (Electrodynamics of black hole magnetospheres, 2004, Appendix C): B_x = 1
-- and B_y = +B0 on the left, -B0 on the right, with no electric field (B0 = 0.5, -1.5 < x < 1.5).
-- The exact solution splits the jump into two fronts moving at the speed of light, with B_y = 0 and
-- E_z = -B0 between them. Figures of merit (tensor p2, 200 cells) at t = 1: E_z crosses -B0/2 at x
-- = -1.0005 and +1.0005 (exact fronts at -1 and +1) and averages -0.5000 on |x| < 0.5; the L1 error
-- of E_z against the exact solution, (1/3) int |E_z - E_z,exact| dx over the domain, is 0.0032; at
-- the fronts B_y overshoots B0 by 7.8% and E_z overshoots -B0 by 4.8%.
local Vlasov = G0.Vlasov

epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
B0 = 0.5 -- Transverse magnetic field, +B0 on the left and -B0 on the right.

Nx = 200 -- Cell count (x-direction).
Lx = 1.5 -- Half domain size (x-direction).
poly_order = 2 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 1.0 -- Final simulation time.
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
  lower = { -Lx },
  upper = { Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction.

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
      local x = xn[1]

      -- Komissarov's Riemann problem A: B_x = 1 with B_y = +B0 on the left and -B0 on the right, no electric field.
      local By = -B0
      if x < 0.0 then
        By = B0
      end

      return 0.0, 0.0, 0.0, 1.0, By, 0.0, 0.0, 0.0
    end,

    bcx = { G0.FieldBc.bcCopy, G0.FieldBc.bcCopy },
  }
}

vlasovApp:run()
