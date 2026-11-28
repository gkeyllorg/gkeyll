-- Langmuir oscillation about a held, nonuniform electron background with kinetic electrons and the
-- Poisson field (1x1v). Electrons (vt = 0.01, lambda_D = 0.01, no ions) with the density n0 (1 +
-- 0.5 cos(k x)), k = 1, are held about that state at rest: the potential is solved from the
-- departure of the electrons from their background, so the background, which is not neutral,
-- exerts no force (its own field would be 5e3 times that of the oscillation), and the app derives
-- the body force that balances its pressure gradient. The initial velocity u0 sin(k x), u0 = 1e-4,
-- starts an oscillation of every fluid element at its local plasma frequency: for a cold plasma
-- u_x = u0 sin(k x) cos(omega_p(x) t) and E_x = (m / |q|) u0 sin(k x) omega_p(x) sin(omega_p(x)
-- t), with omega_p(x)^2 = omega_pe^2 (1 + 0.5 cos(k x)). The oscillation phase-mixes along x:
-- after two periods of omega_pe, the end of the run, the phase differs by 6.5 radians across the
-- domain. Thermal corrections are of the order of (k lambda_D)^2. LBO self-collisions at nu = 0.01
-- omega_pe act on the electrons, as in the PKPM test (they conserve the momentum of the
-- oscillation). The potential is solved in every SSP-RK3 stage (cfl_frac = 0.5, so that the
-- third-order stepper resolves the oscillation: omega_p dt up to 0.16).
-- rt_dg_5m_poisson_langmuir_background_1x_ser_p1.lua and
-- rt_pkpm_poisson_langmuir_background_1x.lua are the same problem with fluid and PKPM electrons.
--
-- Figures of merit at t = 4 pi / omega_pe (64 x 16 cells, serendipity p2, velocity grid +-6 vt, 100
-- steps): the cell averages of E_x differ from the cold local solution by 8.7e-3 in L2 and those
-- of u_x by 1.1e-2. The density changes by at most 3.1e-4 of itself, the compression of the
-- oscillation.

local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi -- Pi.

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

vt = 0.01 -- Thermal velocity.
n0 = 1.0 -- Number density.
a_bg = 0.5 -- Relative amplitude of the density variation of the background.
u0 = 1.0e-4 -- Amplitude of the initial velocity.
k = 1.0 -- Wavenumber.

-- Derived physical quantities (using normalized code units).
omega_pe = math.sqrt(n0 * charge_elc * charge_elc / (epsilon0 * mass_elc)) -- Plasma frequency.
nu_ee = 0.01 * omega_pe -- Electron-electron collision frequency.

-- Simulation parameters.
Nx = 64 -- Cell count (configuration space: x-direction).
Nvx = 16 -- Cell count (velocity space: vx-direction).
Lx = 2.0 * pi / k -- Domain size (configuration space: x-direction).
vx_max = 6.0 * vt -- Domain boundary (velocity space: vx-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.5 -- CFL coefficient.

t_end = 4.0 * pi / omega_pe -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Maxwellian of the electrons with the density of the background and the flow velocity ux.
local function elcMaxwellian (x, vx, ux)
  local n = n0 * (1.0 + a_bg * math.cos(k * x)) -- Number density.

  return n / math.sqrt(2.0 * pi * vt * vt) * math.exp(-(vx - ux) * (vx - ux) / (2.0 * vt * vt))
end

-- Background of the electrons: a density that varies along x, at rest, with a uniform temperature.
function elcBackground (t, xn)
  return elcMaxwellian(xn[1], xn[2], 0.0)
end

-- Initial condition: the background with the velocity u0 sin(k x).
function elcInit (t, xn)
  return elcMaxwellian(xn[1], xn[2], u0 * math.sin(k * xn[1]))
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
  periodicDirs = { 1 }, -- Periodic directions.

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

        init = elcInit
      }
    },

    -- The electrons are held about their background.
    background = {
      projectionID = G0.Projection.Func,

      init = elcBackground
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_ee -- Electron-electron collision frequency.
      end
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Field: the electrostatic potential of the departure of the electrons from their background.
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
    }
  }
}

vlasovApp:run()
