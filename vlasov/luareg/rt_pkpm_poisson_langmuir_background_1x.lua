-- Langmuir oscillation about a held, nonuniform electron background with the PKPM system of
-- equations and the Poisson field (1x, normalized parallel velocity w = v_par/xi). Electrons (vt =
-- 0.01, lambda_D = 0.01, no ions) with the density n0 (1 + 0.5 cos(k x)), k = 1, are held about
-- that state at rest: the potential is solved from the departure of the electrons from their
-- background, so the background, which is not neutral, exerts no force (its own field would be 5e3
-- times that of the oscillation), and the app derives the body force that balances its pressure
-- gradient. The initial velocity u0 sin(k x), u0 = 1e-4, starts an oscillation of every fluid
-- element at its local plasma frequency: for a cold plasma u_x = u0 sin(k x) cos(omega_p(x) t) and
-- E_x = (m / |q|) u0 sin(k x) omega_p(x) sin(omega_p(x) t), with omega_p(x)^2 = omega_pe^2 (1 +
-- 0.5 cos(k x)). The oscillation phase-mixes along x: after two periods of omega_pe, the end of
-- the run, the phase differs by 6.5 radians across the domain. Thermal corrections are of the
-- order of (k lambda_D)^2. A uniform external magnetic field along x gives PKPM its field
-- direction and does not act on the oscillation; LBO self-collisions at nu = 0.01 omega_pe act on
-- the electrons (they conserve the momentum of the oscillation). The coupling is the
-- kick-drift-kick split, whose kicks take E = -grad(phi) (cfl_frac = 0.2).
-- rt_dg_5m_poisson_langmuir_background_1x_ser_p1.lua and rt_vp_langmuir_background_1x1v_ser_p2.lua
-- are the same problem with fluid and kinetic electrons.
--
-- Figures of merit at t = 4 pi / omega_pe (64 x 16 cells, 248 steps, omega_pe dt = 0.051): the
-- cell averages of E_x differ from the cold local solution by 8.4e-3 in L2 and those of u_x by
-- 1.1e-2. The density changes by at most 3.1e-4 of itself, the compression of the oscillation.

local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi -- Pi.

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

vt = 0.01 -- Thermal velocity.
B0 = 1.0 -- Magnetic field along x (sets the field direction b).
n0 = 1.0 -- Number density.
a_bg = 0.5 -- Relative amplitude of the density variation of the background.
u0 = 1.0e-4 -- Amplitude of the initial velocity.
k = 1.0 -- Wavenumber.

-- Derived physical quantities (using normalized code units).
omega_pe = math.sqrt(n0 * charge_elc * charge_elc / (epsilon0 * mass_elc)) -- Plasma frequency.
nu_ee = 0.01 * omega_pe -- Electron-electron collision frequency.

-- Simulation parameters.
Nx = 64 -- Cell count (configuration space: x-direction).
Nvx = 16 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 2.0 * pi / k -- Domain size (configuration space: x-direction).
vx_max = 6.0 -- Domain boundary (velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.2 -- CFL coefficient.

t_end = 4.0 * pi / omega_pe -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Normalized electron distributions [F, G, M0_perp (x, y, z), N_perp (x, y, z)]: a Maxwellian in
-- w = v_par/xi and no first harmonic. The initial condition and the held background.
function elcDist (t, xn)
  local w = xn[2]

  local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * pi) -- Normalized Maxwellian in w = v_par/xi.

  return F, F, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
end

-- Fluid background of the electrons: a density that varies along x, at rest, with a uniform
-- temperature.
function elcFluidBackground (t, xn)
  local x = xn[1]

  local n = n0 * (1.0 + a_bg * math.cos(k * x)) -- Electron number density.

  return mass_elc * n, 0.0, 0.0, 0.0, 1.5 * n * mass_elc * vt * vt
end

-- Initial fluid state: the background with the velocity u0 sin(k x).
function elcFluidInit (t, xn)
  local x = xn[1]

  local rho, _, _, _, E = elcFluidBackground(t, xn)
  local ux = u0 * math.sin(k * x) -- Flow velocity.

  return rho, rho * ux, 0.0, 0.0, E + 0.5 * rho * ux * ux
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
          init = elcDist
        }
      },
      -- The electrons are held about their background: the normalized Maxwellian here, the fluid
      -- background in the fluid block.
      background = {
        projectionID = G0.Projection.Func,
        init = elcDist
      },

      collisions = {
        collisionID = G0.Collisions.LBO,

        selfNu = function (t, xn)
          return nu_ee -- Electron-electron collision frequency.
        end
      }
    },

    fluid = {
      init = elcFluidInit,
      background = elcFluidBackground,
      evalOnNodes = true
    },

    evolve = true -- Evolve species?
  },

  -- Field: the electrostatic potential of the departure of the electrons from their background,
  -- and the static external magnetic field that sets b.
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
