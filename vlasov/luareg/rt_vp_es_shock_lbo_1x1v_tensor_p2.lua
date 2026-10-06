-- Electrostatic shock with the Vlasov-Poisson system of equations.
-- Half domain: a plasma flowing toward a reflecting wall at x = 0 piles up against its own
-- reflection and launches a shock back upstream. The inflow boundary at x = Lx holds the upstream
-- Maxwellians fixed (reservoir). Ions slower than the shock potential are reflected ahead of it.
-- Weak LBO self-collisions act on both species.
-- Figures of merit at t_end (serendipity p2): shock front at x = 58 lambda_D (speed 0.97 c_s),
-- downstream density 2.4 n0, potential jump across the front e*dphi/Te = 1.3-1.4, reflected ions ahead
-- of the front 0.25 n0 at 3.0 c_s.

local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 100.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.
Te_over_Ti = 4.0 -- Electron temperature / ion temperature.
mach_num = 2.0 -- Mach number of the inflow (relative to the ion sound speed).

-- Derived physical quantities (using normalized code units).
Te = 1.0 * charge_ion -- Electron temperature.
Ti = Te / Te_over_Ti -- Ion temperature.

vte = math.sqrt(Te / mass_elc) -- Electron thermal velocity.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.
cs = math.sqrt(Te / mass_ion) -- Ion sound speed.
Vx_drift = -mach_num * cs -- Drift velocity of the inflow (x-direction).

lambda_D = math.sqrt(epsilon0 * Te / (n0 * charge_ion * charge_ion)) -- Electron Debye length.
omega_pe = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_elc)) -- Electron plasma frequency.

nu_ee = 0.002 * omega_pe -- Electron-electron collision frequency.
nu_ii = nu_ee * math.sqrt(mass_elc / mass_ion) * math.pow(Ti / Te, 1.5) -- Ion-ion collision frequency.

-- Simulation parameters.
Nx = 96 -- Cell count (configuration space: x-direction).
Nvx_elc = 16 -- Cell count (electron velocity space: vx-direction).
Nvx_ion = 32 -- Cell count (ion velocity space: vx-direction).
Lx = 192.0 * lambda_D -- Domain size (configuration space: x-direction).
vx_max_elc = 6.0 * vte -- Domain boundary (electron velocity space: vx-direction).
vx_max_ion = 6.0 * cs -- Domain boundary (ion velocity space: vx-direction).
poly_order = 2 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 600.0 / omega_pe -- Final simulation time.
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
  periodicDirs = { }, -- Periodic directions (none).

  -- Electrons.
  elc = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_elc, mass = mass_elc,

    -- Velocity space grid.
    lower = { -vx_max_elc },
    upper = { vx_max_elc },
    cells = { Nvx_elc },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          return n0 -- Electron total number density.
        end,
        temperatureInit = function (t, xn)
          return Te -- Electron isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift -- Electron drift velocity.
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_ee -- Electron collision frequency.
      end
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcReflect
      },
      upper = {
        type = G0.SpeciesBc.bcFixedFunc
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Ions.
  ion = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_ion, mass = mass_ion,

    -- Velocity space grid.
    lower = { -vx_max_ion },
    upper = { vx_max_ion },
    cells = { Nvx_ion },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          return n0 -- Ion total number density.
        end,
        temperatureInit = function (t, xn)
          return Ti -- Ion isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift -- Ion drift velocity.
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_ii -- Ion collision frequency.
      end
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcReflect
      },
      upper = {
        type = G0.SpeciesBc.bcFixedFunc
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  isElectrostatic = true,

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0,

    poissonBcs = {
      lowerType = {
        G0.PoissonBc.bcNeumann
      },
      upperType = {
        G0.PoissonBc.bcDirichlet
      },
      lowerValue = {
        0.0
      },
      upperValue = {
        0.0
      }
    }
  }
}

vlasovApp:run()
