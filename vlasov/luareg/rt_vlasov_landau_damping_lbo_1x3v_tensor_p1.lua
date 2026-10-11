-- Linear Landau damping of a Langmuir wave with the Vlasov-Maxwell system (1x3v).
-- A Maxwellian electron plasma (vt = 1, lambda_D = 1) on a neutralizing ion background is perturbed
-- in density at k lambda_D = 0.5 with amplitude 1e-4; the electric field follows from Gauss's law.
-- The 3V Maxwellian keeps vy and vz on 4 cells each. LBO self-collisions at nu = 0.05 omega_pe act
-- on the electrons.
-- Collisionless linear theory at k lambda_D = 0.5: omega = 1.4157 omega_pe, damping rate gamma =
-- 0.1533 omega_pe.
-- Figures of merit (tensor p1 (hybrid), 16 x 16 x 4 x 4 cells): the field-energy peaks over t =
-- 0.5-12 decay at gamma = 0.1436 and recur at omega = 1.387; the field energy at t = 20 is 0.00059
-- of its initial value. The rate differs from the 1x1v value 0.123 because the 3V collisions also
-- scatter into the perpendicular directions and relax the isotropic temperature.

local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

vt = 1.0 -- Thermal velocity.
nu = 0.05 -- Collision frequency.

alpha = 1.0e-4 -- Applied perturbation amplitude.
k0 = 0.5 -- Perturbed wave number.

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Nvx = 16 -- Cell count (velocity space: vx-direction).
Nvy = 4 -- Cell count (velocity space: vy-direction).
Nvz = 4 -- Cell count (velocity space: vz-direction).
Lx = 4.0 * pi -- Domain size (configuration space: x-direction).
vx_max = 6.0 * vt -- Domain boundary (velocity space: vx-direction).
vy_max = 3.5 * vt -- Domain boundary (velocity space: vy-direction).
vz_max = 3.5 * vt -- Domain boundary (velocity space: vz-direction).
poly_order = 1 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 20.0 -- Final simulation time.
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
  lower = { -0.5 * Lx },
  upper = { 0.5 * Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Electrons.
  elc = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_elc, mass = mass_elc,
    
    -- Velocity space grid.
    lower = { -vx_max, -vy_max, -vz_max },
    upper = { vx_max, vy_max, vz_max },
    cells = { Nvx, Nvy, Nvz },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local x, vx, vy, vz = xn[1], xn[2], xn[3], xn[4]

          local v_sq = (vx * vx) + (vy * vy) + (vz * vz)
          local n = (1.0 + alpha * math.cos(k0 * x)) *
            (1.0 / math.pow(math.sqrt(2.0 * pi * vt * vt), 3.0)) * (math.exp(-v_sq / (2.0 * vt * vt))) -- Distribution function.

          return n
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu -- Collision frequency.
      end,
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local Ex = -alpha * math.sin(k0 * x) / k0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()