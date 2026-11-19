-- Electrons in an external electrostatic potential well with the Vlasov-Poisson system (1x1v). A
-- uniform Maxwellian (n = 1, v_t = 1) starts in the periodic external potential phi_ext = -cos(x),
-- i.e. the field E = -sin(x) of the original Vlasov-Maxwell test, and the self-consistent potential
-- of the plasma response (omega_p = 1, Debye length 1) is solved from Poisson's equation on top of
-- it. Electrons bunch at the bottom of the well and partially screen it. Particle number is
-- conserved to round-off and the total energy (kinetic, external potential q int n phi_ext from
-- frames 0 and 1, and self-consistent field energy epsilon0/2 int |grad phi|^2) to 1e-8 relative.
-- By t = 3 the density ranges from 0.65 to 1.37 and the peak-to-peak self-consistent potential
-- reaches 0.37 of the external one.
local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

-- Simulation parameters.
Nx = 32 -- Cell count (configuration space: x-direction).
Nvx = 24 -- Cell count (velocity space: vx-direction).
Lx = 2.0 * pi -- Domain size (configuration space: x-direction).
vx_max = 6.0 -- Domain boundary (velocity space: vx-direction).
poly_order = 2 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 3.0 -- Final simulation time.
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
  lower = { 0.0 },
  upper = { Lx },
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
    lower = { -vx_max },
    upper = { vx_max },
    cells = { Nvx },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local vx = xn[2]

          local n = (1.0 / math.sqrt(2.0 * pi)) * (math.exp(-(vx * vx) / 2.0)) -- Distribution function.
          return n
        end
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Field.
  -- Field: Poisson solve with the external potential well.
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

    -- External potentials (phi, A_x, A_y, A_z): the well phi_ext = -cos(x).
    externalPotentialInit = function (t, xn)
      local x = xn[1]

      local phi = -math.cos(x) -- External electrostatic potential (the well).

      return phi, 0.0, 0.0, 0.0
    end,
    evolveExternalPotential = false
  }
}

vlasovApp:run()