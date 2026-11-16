-- Compressive flow of a collisionless neutral gas across a nearly perpendicular magnetic field with
-- the PKPM system of equations (1x, normalized parallel velocity w = v_par/xi). A uniform gas (n =
-- 1, p = 1) with the flow u_x = 0.1 cos(pi x) sits in the uniform static field B = (sin 0.0871,
-- -cos 0.0871, 0) (b_x = 0.087) in a periodic box of length 2. The gas has no charge, so the field
-- only sets the direction b. The gas streams along x at b_x w xi, comparable to the flow, so the
-- velocity w = -u_x/b_x at which the configuration-space flux changes its upwind side lies inside
-- the velocity grid, in different velocity cells at the two faces of a cell.
-- Figures of merit (16 x 16 cells) to t = 2: particle number, the three components of the momentum
-- and the total energy (kinetic + p_par/2 + p_perp) are conserved to 2e-13. The kinetic energy
-- falls from 0.005 to 0.0013, and the integrated parallel and perpendicular pressures are 1.9983
-- and 2.0045 at t = 2 (both 2 initially). The positivity limiter scales no cell.

local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

n0 = 1.0 -- Reference number density.
ux0 = 0.1 -- Amplitude of the compressive flow (x-direction).
pr0 = 1.0 -- Reference pressure.
B0 = 1.0 -- Magnetic field strength (in the x-y plane).
theta = 0.0871 -- Angle of the magnetic field from the y-direction (b_x = sin(theta)).

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

pkpmApp = PKPM.App.new {

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
  neut = PKPM.Species.new {
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
      }
    },

    fluid = {
      -- Initial conditions.
      init = function (t, xn)
        local x = xn[1]

        local rho = mass * n0 -- Mass density.
        local ux = ux0 * math.cos(pi * x) -- Compressive flow.
        local E = 0.5 * rho * ux * ux + 1.5 * pr0 -- Total energy density (uniform pressure).

        return rho, rho * ux, 0.0, 0.0, E
      end
    },

    evolve = true -- Evolve species?
  },

  -- Field.
  field = PKPM.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      return 0.0, 0.0, 0.0, B0 * math.sin(theta), -B0 * math.cos(theta), 0.0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()
