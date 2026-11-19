-- Pulse leaving a vacuum box through absorbing boundaries of the Maxwell field (1x). A Gaussian
-- pulse of width 0.08 at the centre of [0, 1] in E_x, E_y, E_z and B_x, with c = 1 and both
-- divergence-error speeds equal to c: E_y and E_z split into plane waves moving left and right, and
-- E_x and B_x into waves of the divergence correction, (E_x, phi) and (B_x, psi). Both boundaries
-- are GKYL_FIELD_ABSORB: the flux through the boundary face is that of the outgoing characteristic
-- variables of the cell next to it (the interior faces keep the central flux), for the light waves
-- and for the correction waves, and no ghost cell enters it.
--
-- Figures of merit (64 cells, serendipity p1, 427 steps): the energy of E and B is 2.4e-6 of the
-- initial one at t = 1, when the pulses have left, and 5.3e-9 at t = 2; the eight components are
-- finite and the largest cell average is 4.8e-5 at t = 2, against 0.987 at t = 0.

local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.

x0 = 0.5 -- Centre of the pulse.
width = 0.08 -- Width of the pulse.
right_frac = 0.5 -- Fraction of the pulse that travels to the right.

-- Simulation parameters.
Nx = 64 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
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
  lower = { 0.0 },
  upper = { Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions.

  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,
    elcErrorSpeedFactor = 1.0,
    mgnErrorSpeedFactor = 1.0,

    -- Initial conditions function (E, B, phi, psi).
    init = function (t, xn)
      local x = xn[1]

      local c = 1.0 / math.sqrt(epsilon0 * mu0)
      local pulse = math.exp(-(x - x0) * (x - x0) / (width * width))
      local fr, fl = right_frac, 1.0 - right_frac

      -- Right-going part: (E_y, B_z) = (1, 1/c) pulse; left-going part: (1, -1/c) pulse. The same
      -- pulses in (E_z, B_y), with the opposite sign of B, and in the cleaning pairs (E_x, phi) and
      -- (B_x, psi).
      return (fr + fl) * pulse, (fr + fl) * pulse, (fr + fl) * pulse,
        (fr + fl) * pulse / c, -(fr - fl) * pulse / c, (fr - fl) * pulse / c,
        (fr - fl) * pulse / c, (fr - fl) * pulse
    end,

    bcx = { G0.FieldBc.bcAbsorb, G0.FieldBc.bcAbsorb },
  },
}

vlasovApp:run()
