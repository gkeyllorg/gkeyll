-- Jeans instability of a self-gravitating collisionless medium with the Vlasov-Poisson system
-- (1x1v). A Maxwellian of neutral dark matter (vt = 1, m = 1, n0 = 1) with gravitational coupling
-- alpha_g = 4 pi G = 1, so that omega_J = sqrt(alpha_g m n0) = 1 and the Jeans wave number is k_J =
-- omega_J/vt = 1, is perturbed in density at k/k_J = 0.5 with amplitude 1e-4. The species is
-- neutral, so only the self-gravity Poisson equation nabla^2 phi_g = alpha_g m n is solved (the
-- periodic solve removes the mean density, i.e. the Jeans swindle).
-- Kinetic linear theory at k/k_J = 0.5 (Landau's dispersion relation with omega_pe^2 ->
-- -omega_J^2): a purely growing mode with gamma = 0.6872 omega_J, so the field energy grows at 2
-- gamma = 1.374 omega_J.
-- Figures of merit (tensor p2, 32 x 32 cells): the gravitational field energy grows at 1.3754 (fit
-- of ln E over t = 2-8) from 2.51e-07 at t = 0 to 0.1065 at t = 10; mass conserved to 5e-14 and
-- total energy (kinetic minus half the gravitational energy file, int |grad phi_g|^2/alpha_g) to
-- 8e-10.

local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
alpha_g = 1.0 -- Gravitational coupling 4 pi G.
mass_dm = 1.0 -- Dark matter mass.
charge_dm = 0.0 -- Dark matter charge (neutral).

n0 = 1.0 -- Reference number density.
vt = 1.0 -- Thermal velocity.

alpha = 1.0e-4 -- Applied perturbation amplitude.

-- Derived physical quantities (using normalized code units).
omega_J = math.sqrt(alpha_g * mass_dm * n0) -- Jeans frequency.
k_J = omega_J / vt -- Jeans wave number.

k0 = 0.5 * k_J -- Perturbed wave number.

-- Simulation parameters.
Nx = 32 -- Cell count (configuration space: x-direction).
Nvx = 32 -- Cell count (velocity space: vx-direction).
Lx = 2.0 * pi / k0 -- Domain size (configuration space: x-direction).
vx_max = 6.0 * vt -- Domain boundary (velocity space: vx-direction).
poly_order = 2 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 10.0 / omega_J -- Final simulation time.
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

  -- Dark matter.
  dm = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_dm, mass = mass_dm,

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
          local x, vx = xn[1], xn[2]

          local n = n0 * (1.0 + alpha * math.cos(k0 * x)) *
            (1.0 / math.sqrt(2.0 * pi * vt * vt)) * (math.exp(-(vx * vx) / (2.0 * vt * vt))) -- Distribution function.

          return n
        end
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  isElectrostatic = true,

  -- Field: self-gravity only (the species is neutral, so no electrostatic solve).
  field = Vlasov.Field.new {
    alphaG = alpha_g,

    poissonBcs = {
      lowerType = {
        G0.PoissonBc.bcPeriodic
      },
      upperType = {
        G0.PoissonBc.bcPeriodic
      }
    }
  }
}

vlasovApp:run()
