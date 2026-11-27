-- Propagation into a plasma wave beach with kinetic electrons and the Maxwell field (1x2v), the
-- kinetic counterpart of vlasov/creg/rt_dg_5m_mom_beach_1x.c, about a held background density.
-- The setup follows entry JE8 of Ammar's Simulation Journal
-- (https://ammar-hakim.org/sj/je/je8/je8-plasmabeach.html), adapted from Section III. A. of D. N.
-- Smithe (2007), "Finite-difference time-domain simulation of fusion plasmas at radiofrequency time
-- scales", Physics of Plasmas 14 (5): 056104: the plasma frequency falls as omega_p deltaT = 25 (1
-- - x)^5 with deltaT = 0.01 / c, here on a uniform background with omega_p deltaT = 0.025, and an
-- antenna current J_y = -J_0 sin(omega_drive t) on 1 - 1/400 < x < 1 drives a wave with
-- omega_drive deltaT = pi / 10 toward the cutoff omega_p = omega_drive at x = 0.583. The upper
-- boundary of the field, next to the antenna, is absorbing (G0.FieldBc.bcAbsorb). The electron
-- current is explicit in the SSP-RK3 stages, so the time step resolves the plasma frequency of
-- every cell (omega_p dt = 0.9 in the cell at x = 0), where the 5-moment beach, with its implicit
-- source step, takes a quarter of the steps. The Maxwellian electrons (1 eV, velocity grid +-6
-- v_te) are held about their initial state
-- (First_Fourier_Harmonic_Notes/background_density_gradient.tex): the app derives the body force
-- that holds them and subtracts what the discretization leaves. Without the hold the pressure
-- gradient pushes the electrons until their space charge stops them, and E_x rings about the
-- ambipolar field -(T_e/e) d ln n/dx at the local plasma frequency, up to 54 V/m. The response is
-- linear in the antenna current: J_0 = 1 A/m^3 (1e-12 in the 5-moment beach) puts the electron
-- drift at 2.4e-5 v_te, above the round-off of the current of a distribution function.
--
-- Figures of merit at t = 5e-9 (100 x 8 x 8 cells, serendipity p2, 4062 steps), from the cell
-- averages: on 0.3 < x < 0.95, E_y is that of the 5-moment beach, per unit antenna current, times
-- 1.1105 to 3.9e-3 in L2. The factor is the 10/9 of its current that the 3-point quadrature of
-- the applied current gives the antenna, which covers a quarter of the last cell. The maximum of
-- |E_y| is 0.992 V/m at x = 0.625, 2.11 times the amplitude mu_0 c J_0 / 800 the antenna radiates
-- to each side in free space (1.90 times that of the projected current, as in the 5-moment beach),
-- and |E_y| stays below 0.3% of the maximum for x < 0.485, beyond the cutoff. |E_x| stays below
-- 1e-7 V/m, the density at its initial value to 6e-11 and |u_x| below 3e-11 v_te. 12 x 12
-- velocity cells give the same E_y to four digits.

local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi -- Pi.

-- Physical constants (using non-normalized physical units).
epsilon0 = 8.854187817620389850536563031710750260608e-12 -- Permittivity of free space.
mu0 = 12.56637061435917295385057353311801153679e-7 -- Permeability of free space.
mass_elc = 9.10938215e-31 -- Electron mass.
charge_elc = -1.602176487e-19 -- Electron charge.

Te = 1.0 * (-charge_elc) -- Electron temperature (1 eV).
J0 = 1.0 -- Reference current density (Amps / m^3).
omega_p_floor = 0.025 -- Plasma frequency of the uniform background (times deltaT).

-- Derived physical quantities (using non-normalized physical units).
light_speed = 1.0 / math.sqrt(mu0 * epsilon0) -- Speed of light.
vte = math.sqrt(Te / mass_elc) -- Electron thermal velocity.

-- Simulation parameters.
Nx = 100 -- Cell count (configuration space: x-direction).
Nvx = 8 -- Cell count (velocity space: vx-direction).
Nvy = 8 -- Cell count (velocity space: vy-direction).
Lx = 1.0 -- Domain size (configuration space: x-direction).
vx_max = 6.0 * vte -- Domain boundary (velocity space: vx-direction).
vy_max = 6.0 * vte -- Domain boundary (velocity space: vy-direction).
Lx100 = Lx / 100.0 -- Domain size over 100 (x-direction).
x_last_edge = Lx - Lx / 400.0 -- Lower edge of the antenna, which reaches the upper boundary (low density side).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

deltaT = Lx100 / light_speed -- Arbitrary constant, with units of time.
factor = deltaT * deltaT * charge_elc * charge_elc / (mass_elc * epsilon0) -- Numerical factor for calculation of electron number density.
omega_drive = pi / 10.0 / deltaT -- Drive current angular frequency.

t_end = 5.0e-9 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Electron distribution function: a Maxwellian at rest with the plasma frequency profile (times
-- deltaT) on the uniform background. The initial condition and the held background.
function elcDist (t, xn)
  local x, vx, vy = xn[1], xn[2], xn[3]

  local omegaPdt = 25.0 * (1.0 - x) * (1.0 - x) * (1.0 - x) * (1.0 - x) * (1.0 - x)
  local ne = (omegaPdt * omegaPdt + omega_p_floor * omega_p_floor) / factor -- Electron number density.

  return ne / (2.0 * pi * vte * vte) * math.exp(-(vx * vx + vy * vy) / (2.0 * vte * vte))
end

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

  -- Electrons.
  elc = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_elc, mass = mass_elc,

    -- Velocity space grid.
    lower = { -vx_max, -vy_max },
    upper = { vx_max, vy_max },
    cells = { Nvx, Nvy },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = elcDist,
      },
    },

    -- The initial state is held as the background.
    background = {
      projectionID = G0.Projection.Func,

      init = elcDist,
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcCopy
      },
      upper = {
        type = G0.SpeciesBc.bcCopy
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    init = function (t, xn)
      return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    appliedCurrent = function (t, xn)
      local x = xn[1]

      -- Applied current (y-direction): the antenna next to the upper boundary.
      local app_y = 0.0
      if x > x_last_edge then
        app_y = -J0 * math.sin(omega_drive * t)
      end

      return 0.0, app_y, 0.0
    end,
    evolveAppliedCurrent = true,

    bcx = { G0.FieldBc.bcCopy, G0.FieldBc.bcAbsorb },
  },
}

vlasovApp:run()
