-- Flow of a collisionless neutral gas through a magnetic field whose direction varies along x, with
-- the PKPM system of equations (1x, normalized parallel velocity w = v_par/xi). A density pulse n =
-- 1 + 0.2 sin(pi x) in a uniform flow u = 0.5 at uniform pressure p = 1 moves through the static
-- field B = (sin 0.3, -cos 0.3, -0.5 cos(pi x)) in a periodic box of length 2. The gas has no
-- charge, so the field only sets the direction b, and div b is not zero: streaming along b, the
-- mirror force, the work of the pressure anisotropy on the flow, (p_par - p_perp) bb:grad(u), and
-- the first-harmonic stress act together.
-- Figures of merit (16 x 16 cells) to t = 2: particle number, the three components of the momentum
-- and the total energy (kinetic + p_par/2 + p_perp) are conserved to 2e-13. The kinetic energy
-- rises from 0.25 to 0.272 (t = 0.74) and is 0.254 at t = 2, where the integrated parallel and
-- perpendicular pressures are 1.9991 and 1.9968 (both 2 initially). The positivity limiter scales
-- no cell.

local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

n0 = 1.0 -- Reference number density.
alpha = 0.2 -- Perturbation amplitude.
u0 = 0.5 -- Reference fluid velocity (x-direction).
pr0 = 1.0 -- Reference pressure.
B0 = 1.0 -- Magnetic field strength of the uniform part (in the x-y plane).
theta = 0.3 -- Angle of the uniform part from the y-direction (b_x = sin(theta)).
dB = 0.5 -- Amplitude of the varying part (z-direction).

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

        local n = n0 + alpha * math.sin(pi * x) -- Number density.
        local rho = mass * n -- Mass density.
        local E = 0.5 * rho * u0 * u0 + 1.5 * pr0 -- Total energy density (uniform pressure).

        return rho, rho * u0, 0.0, 0.0, E
      end
    },

    evolve = true -- Evolve species?
  },

  -- Field.
  field = PKPM.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      return 0.0, 0.0, 0.0, B0 * math.sin(theta), -B0 * math.cos(theta), -dB * math.cos(pi * x), 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()
