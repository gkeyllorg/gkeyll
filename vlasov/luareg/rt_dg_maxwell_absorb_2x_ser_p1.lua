-- Pulse leaving a vacuum box through absorbing boundaries of the Maxwell field on its four sides
-- (2x). A Gaussian pulse of width 0.08 at the centre of [0, 1]^2 in E_z and B_z (the two
-- polarizations of a cylindrical wave) and in E_x, whose divergence the correction potential
-- carries away, with c = 1 and both divergence-error speeds equal to c. The four boundaries are
-- GKYL_FIELD_ABSORB: the flux through a boundary face is that of the outgoing characteristic
-- variables normal to the face, which is exact for a wave at normal incidence and reflects part of
-- a wave that meets the face at an angle.
--
-- Figures of merit (48 x 48 cells, serendipity p1, 640 steps): the energy of E and B is 1.1e-2 of
-- the initial one at t = 1, 5.2e-3 at t = 1.5 and 3.8e-4 at t = 2; the eight components are finite
-- and the largest cell average is 6.8e-3 at t = 2, against 0.956 at t = 0. With 96 x 96 cells the
-- energy left is 1.1e-2, 5.3e-3 and 3.9e-4 at the same times: what is left is set by the reflection
-- at oblique incidence and not by the resolution. With copy boundaries in place of the absorbing
-- ones (GKYL_FIELD_COPY, central flux) the energy reaches 20 times the initial one and is 13 times
-- it at t = 2.

local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.

x0 = 0.5 -- Centre of the pulse (x-direction).
y0 = 0.5 -- Centre of the pulse (y-direction).
width = 0.08 -- Width of the pulse.

-- Simulation parameters.
Nx = 48 -- Cell count (x-direction).
Ny = 48 -- Cell count (y-direction).
Lx = 1.0 -- Domain size (x-direction).
Ly = 1.0 -- Domain size (y-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 2.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = 200 -- Number of times to calculate field energy.
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
  periodicDirs = { }, -- Periodic directions.

  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,
    elcErrorSpeedFactor = 1.0,
    mgnErrorSpeedFactor = 1.0,

    -- Initial conditions function (E, B, phi, psi).
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local c = 1.0 / math.sqrt(epsilon0 * mu0)
      local r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0)
      local pulse = math.exp(-r2 / (width * width))

      -- A pulse of E_z and one of B_z (the two polarizations of a cylindrical wave), and a pulse
      -- of E_x, whose divergence the correction potential phi carries out of the box.
      return pulse, 0.0, pulse, 0.0, 0.0, pulse / c, 0.0, 0.0
    end,

    bcx = { G0.FieldBc.bcAbsorb, G0.FieldBc.bcAbsorb },
    bcy = { G0.FieldBc.bcAbsorb, G0.FieldBc.bcAbsorb },
  },
}

vlasovApp:run()
