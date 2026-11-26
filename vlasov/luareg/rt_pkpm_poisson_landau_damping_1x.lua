-- Linear Landau damping of a Langmuir wave with the PKPM system of equations and the Poisson field
-- (1x, normalized parallel velocity w = v_par/xi), the Poisson counterpart of
-- rt_pkpm_landau_damping_1x.c and of rt_vp_landau_damping_1x1v.c. A Maxwellian electron plasma
-- (vt = 1, lambda_D = 1) on a neutralizing ion background, which the periodic solve removes as the
-- mean charge density, is perturbed in density at k lambda_D = 0.5 with amplitude 1e-4. A uniform
-- external magnetic field along x gives PKPM its field direction and does not act on the parallel
-- electrostatic wave. The fluid takes E = -grad(phi) in the kicks of the kick-drift-kick split
-- between the SSP-RK3 transport steps, with phi solved from the density at the start of each kick.
-- Collisionless linear theory at k lambda_D = 0.5: omega = 1.4157 omega_pe, damping rate gamma =
-- 0.1533 omega_pe.
-- Figures of merit (32 x 32 cells, 931 steps at the CFL step, omega_pe dt = 0.021): the
-- field-energy peaks over t = 0.5-12 decay at gamma = 0.1556 and recur at omega = 1.413; the field
-- energy at t = 20 is 0.00086 of its initial value (the Maxwell twin: 0.1549, 1.414, 0.00090).

local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

vt = 1.0 -- Thermal velocity.
B0 = 1.0 -- Magnetic field along x (sets the field direction b).
alpha = 1.0e-4 -- Applied perturbation amplitude.
k0 = 0.5 -- Perturbed wave number.

-- Simulation parameters.
Nx = 32 -- Cell count (configuration space: x-direction).
Nvx = 32 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 4.0 * pi -- Domain size (configuration space: x-direction).
vx_max = 6.0 -- Domain boundary (velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 20.0 -- Final simulation time.
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
  lower = { -0.5 * Lx },
  upper = { 0.5 * Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coordinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

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
      }
    },

    fluid = {
      -- Initial conditions.
      init = function (t, xn)
        local x = xn[1]

        local n = 1.0 + alpha * math.cos(k0 * x) -- Electron number density.
        local rho = mass_elc * n -- Electron mass density.
        local E = 1.5 * n * mass_elc * vt * vt -- Electron total energy density (at rest).

        return rho, 0.0, 0.0, 0.0, E
      end
    },

    evolve = true -- Evolve species?
  },

  -- Field.
  -- Field: the electrostatic potential (the periodic solve removes the mean charge density, the
  -- ion background) and the static external magnetic field that sets b.
  field = Vlasov.Field.new {
    fieldID = G0.FieldModel.Phi,
    epsilon0 = epsilon0,

    poissonBcs = {
      lowerType = {
        G0.PoissonBc.bcPeriodic
      },
      upperType = {
        G0.PoissonBc.bcPeriodic
      }
    },

    externalFieldInit = function (t, xn)
      return 0.0, 0.0, 0.0, B0, 0.0, 0.0
    end
  }
}

vlasovApp:run()
