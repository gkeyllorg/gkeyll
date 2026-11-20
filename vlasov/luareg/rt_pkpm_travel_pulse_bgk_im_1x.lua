-- Advection of an entropy (contact) wave by a collisional neutral gas with BGK collisions with the
-- PKPM system of equations (1x, normalized parallel velocity w = v_par/xi). A density pulse n = 1 +
-- 0.2 sin(pi x) in a uniform flow u = 1 at uniform pressure p = 1 crosses a periodic box of length
-- 2 once (t = 2), along the magnetic field.
-- In the Euler limit the pulse returns unchanged to its starting position. With collisions at nu =
-- 100 the heat flux relaxes at nu, so the gas conducts heat with thermal diffusivity chi = T/(m nu)
-- and the pulse, a temperature perturbation at constant pressure, decays as exp(-chi k^2 t), k^2 =
-- pi^2. The BGK relaxation takes one backward-Euler step in each SSP-RK3 stage, over the time weight
-- of the stage. Every stage has the relaxed heat flux as its fixed point, so the conduction does not
-- depend on the time step.
-- Figures of merit (16 x 16 cells): at t = 2 the pulse is back at its starting position (phase
-- error below 2e-3 rad) with amplitude 0.836, against exp(-chi k^2 t) = 0.821. At nu dt = 0.26 the
-- decay rate is 0.9996 times that of rt_pkpm_travel_pulse_bgk_1x (amplitude 0.836). Particle number
-- is conserved to 1e-13 and the constraints to round-off.

local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

n0 = 1.0 -- Reference number density.
alpha = 0.2 -- Perturbation amplitude.
u0 = 1.0 -- Reference fluid velocity (each direction).
pr0 = 1.0 -- Reference pressure.
nu = 100.0 -- Collision frequency.
B0 = 1.0 -- Magnetic field along x (sets the field direction b).

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Nvx = 16 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 2.0 -- Domain size (configuration space: x-direction).
vx_max = 6.0 -- Domain boundary (velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 2.0 -- Final simulation time.
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
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Neutrals.
  neut = Vlasov.PKPMSpecies.new {
    charge = charge, mass = mass,

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

      collisions = {
        collisionID = G0.Collisions.BGK,
        useImplicitCollisionScheme = true,

        selfNu = function (t, xn)
          return nu
        end
      }
    },

    fluid = {
      -- Initial conditions.
      init = function (t, xn)
        local x = xn[1]

        local n = n0 + alpha * math.sin(pi * x) -- Number density.
        local rho = mass * n -- Mass density.
        local E = 0.5 * rho * u0 * u0 + 1.5 * pr0 -- Total energy density (uniform pressure).

        return rho, rho * u0, 0.0, 0.0, E
      end
    },

    evolve = true -- Evolve species?
  },

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      return 0.0, 0.0, 0.0, B0, 0.0, 0.0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()
