-- Propagation into a plasma wave beach with PKPM electrons in a magnetic field along the direction
-- of propagation (1x, normalized parallel velocity w = v_par/xi), the magnetized counterpart of
-- vlasov/creg/rt_dg_5m_mom_beach_1x.c. The plasma frequency falls as omega_p deltaT = 25 (1 - x)^5
-- with deltaT = 0.01 / c, on a uniform background with omega_p deltaT = 0.025, and an antenna
-- current J_y = -J_0 sin(omega_drive t) on 1 - 1/400 < x < 1 drives a wave with omega_drive deltaT
-- = pi / 10 into the density gradient. The external field B_x gives Omega_ce = 0.5 omega_drive: the
-- wave splits into its two circular polarizations, with the cutoff of the right-handed one at x =
-- 0.611 (omega_p^2 = omega (omega - Omega_ce)) and of the left-handed one at x = 0.566 (omega_p^2 =
-- omega (omega + Omega_ce)), and its plane of polarization rotates. J_0 is such that the wave the
-- antenna radiates to each side in free space has a magnetic field of 0.02 B_x, so the magnetic
-- unit vector changes in time up to the edge. The upper boundary of the field, next to the antenna,
-- is absorbing (G0.FieldBc.bcAbsorb). The fluid-field coupling is the staggered operator split with
-- the implicit fluid-Maxwell source step, which takes the applied current.
--
-- Figures of merit at t = 5e-9 (200 x 8 cells, serendipity p1, 1000 steps), from the cell averages:
-- the maximum of |E_perp| = sqrt(E_y^2 + E_z^2) is 2.79e5 V/m at x = 0.6275, 1.74 times the
-- amplitude mu_0 c J_0 / 800 the antenna radiates to each side in free space; the maxima of |E_y|
-- and |E_z| are 1.33 and 1.17 times that amplitude. |E_perp| stays below 0.7% of its maximum for x
-- < 0.466, beyond both cutoffs. The norm of the first-harmonic blocks of the distribution is at
-- most 0.42 of that of F (at x = 0.9075), and 0.24 in the cell at the absorbing edge.

local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using non-normalized physical units).
epsilon0 = 8.854187817620389850536563031710750260608e-12 -- Permittivity of free space.
mu0 = 12.56637061435917295385057353311801153679e-7 -- Permeability of free space.
mass_elc = 9.10938215e-31 -- Electron mass.
charge_elc = -1.602176487e-19 -- Electron charge.

Te = 1.602176487e-19 -- Electron temperature (1 eV).
omega_c_frac = 0.5 -- Electron cyclotron frequency over the drive frequency.
dB_B0 = 0.02 -- Amplitude of the wave the antenna radiates in free space, over B0.
omega_p_floor = 0.025 -- Plasma frequency of the uniform background (times deltaT).

-- Derived physical quantities (using non-normalized physical units).
light_speed = 1.0 / math.sqrt(mu0 * epsilon0) -- Speed of light.

-- Simulation parameters.
Nx = 200 -- Cell count (configuration space: x-direction).
Nvx = 8 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 1.0 -- Domain size (configuration space: x-direction).
Lx100 = Lx / 100.0 -- Domain size over 100 (x-direction).
x_last_edge = Lx - Lx / 400.0 -- Lower edge of the antenna, which reaches the upper boundary.
vx_max = 8.0 -- Domain boundary (velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

deltaT = Lx100 / light_speed -- Arbitrary constant, with units of time.
factor = deltaT * deltaT * charge_elc * charge_elc / (mass_elc * epsilon0) -- Numerical factor for calculation of electron number density.
omega_drive = pi / 10.0 / deltaT -- Drive current angular frequency.

B0 = omega_c_frac * omega_drive * mass_elc / math.abs(charge_elc) -- External magnetic field.
-- A current sheet of width Lx - x_last_edge radiates B = mu0 J0 (Lx - x_last_edge) / 2 to each
-- side in free space.
J0 = 2.0 * dB_B0 * B0 / (mu0 * (Lx - x_last_edge)) -- Antenna current density.

t_end = 5.0e-9 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate integrated L2 norm of distribution function.
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
  decompCuts = { 1 }, -- Cuts in each coordinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Electrons.
  elc = Vlasov.PKPMSpecies.new {
    charge = charge_elc, mass = mass_elc,

    kinetic = {
      -- Velocity space grid.
      lower = { -vx_max },
      upper = { vx_max },
      cells = { Nvx },

      -- Initial conditions.
      numInit = 1,
      projections = {
        {
          projectionID = G0.Projection.Func,

          init = function (t, xn)
            local w = xn[2]

            local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * pi) -- Normalized Maxwellian in w = v_par/xi.

            return F, F, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
          end
        }
      },

      bcx = {
        lower = {
          type = G0.SpeciesBc.bcCopy
        },
        upper = {
          type = G0.SpeciesBc.bcCopy
        }
      }
    },

    fluid = {
      -- Initial conditions.
      init = function (t, xn)
        local x = xn[1]

        -- Plasma frequency profile (times deltaT), on the uniform background.
        local omegaPdt = 25.0 * (1.0 - x) * (1.0 - x) * (1.0 - x) * (1.0 - x) * (1.0 - x)
        local ne = (omegaPdt * omegaPdt + omega_p_floor * omega_p_floor) / factor -- Electron number density.

        return mass_elc * ne, 0.0, 0.0, 0.0, 1.5 * ne * Te
      end,

      bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }
    },

    evolve = true -- Evolve species?
  },

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    externalFieldInit = function (t, xn)
      return 0.0, 0.0, 0.0, B0, 0.0, 0.0
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

    bcx = { G0.FieldBc.bcCopy, G0.FieldBc.bcAbsorb }, -- Open boundary and absorbing boundary (x-direction).

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()
