-- Buneman instability with kinetic electrons and a DG Euler fluid for the ions (1x1v), the hybrid
-- twin of vlasov/luareg/rt_vlasov_buneman_1x1v_tensor_p1.lua. Parameters follow entry JE33 of
-- Ammar's Simulation Journal (https://ammar-hakim.org/sj/je/je33/je33-buneman.html): m_i/m_e = 25,
-- electron drift 0.159 so that k v_d = omega_pe for the box mode k = 2 pi, vte = 0.02, density
-- perturbation 1e-5 of the electrons with the electric field from Gauss's law. The ions are a fluid
-- (gamma = 5/3, HLLC flux) with the temperature of the kinetic test (vti = 0.001), coupled to the
-- field by the implicit source step of the staggered operator split; the electrons deposit their
-- current in that step as a known current. The ghost current cancels the net current of the
-- electrons and of the ion fluid, so that Ampere's law reduces to the periodic Poisson problem (see
-- JE33). LBO self-collisions act on the electrons at nu_ee = 1e-4 omega_pe, the frequency of the
-- PKPM tests of this problem (vlasov/luareg/rt_pkpm_hybrid_buneman_1x.lua and
-- rt_pkpm_hybrid_buneman_explicit_1x.lua). Collisions between a kinetic species and a fluid species
-- are not supported, so there is no electron-ion friction here, while the PKPM electrons also
-- exchange momentum and energy with the ion fluid (nu_ei = nu_ee), which heats them by 16% up to t
-- = 40. In one velocity dimension the self-collisions conserve the drift and the temperature of the
-- electrons: the figures below are close to those without collisions and are not the values the
-- PKPM tests should reproduce.
--
-- Collisionless linear theory: 1 + chi_e + chi_i = 0 with chi_e = (1 + zeta Z(zeta)) / (k vte /
-- omega_pe)^2, zeta = (omega - k v_d) / (sqrt(2) k vte), of the drifting Maxwellian and chi_i =
-- -omega_pi^2 / (omega^2 - gamma k^2 vti^2) of the adiabatic fluid gives the growth rate 0.19764
-- omega_pe and omega_r = 0.14251 omega_pe at k = 2 pi (0.19758 and 0.14254 with kinetic ions).
--
-- Figures of merit (tensor p1, hybrid basis of the electrons, 16 x 64 cells, 4114 steps, omega_pe
-- dt = 0.021, the CFL limit of the light wave): the field energy epsilon_0 / 2 int E^2 dx grows at
-- 0.1975 omega_pe over t = 20-50 (-0.1%; 0.1972 over t = 25-55), and at 0.1946 over t = 20-65, the
-- window of the kinetic test (0.1900 in the tensor p1 twin of that test on the same grid, which has
-- the electron-ion friction; 0.1945 without collisions). It peaks at 3.1e-3 at t = 72.1 and is
-- 6.8e-4 at t = 80 (that twin: 3.1e-3 at t = 73.0 and 5.1e-4). Without collisions the three rates
-- are the same to four digits and the peak is 0.4% higher, at t = 72.2. The total energy (field and
-- electrons from the histories, ion fluid from the frames) changes by 3e-10 up to t = 40 and by
-- 4.9e-6 up to t = 80, the electron number by 3e-13, and the mean of E_x stays below 1e-13 of its
-- maximum. The total momentum changes by 7e-5 up to t = 60 and by +14.3% through saturation (+0.6%
-- in the kinetic twin).

local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
mass_ion = 25.0 -- Ion mass.
charge_elc = -1.0 -- Electron charge.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.
Vx_drift_elc = 0.159 -- Electron drift velocity.
vte = 0.02 -- Electron thermal velocity.
vti = 0.001 -- Ion thermal velocity.
gas_gamma = 5.0 / 3.0 -- Adiabatic index of the ion fluid.

alpha = 1.0e-5 -- Applied perturbation amplitude.

-- Derived physical quantities (using normalized code units).
Te = (vte * vte) * mass_elc -- Electron temperature.
Ti = (vti * vti) * mass_ion -- Ion temperature.
omega_pe = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_elc)) -- Electron plasma frequency.
nu_ee = 1.0e-4 * omega_pe -- Electron-electron collision frequency.

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Nvx = 64 -- Cell count (velocity space: vx-direction).
Lx = 1.0 -- Domain size (configuration space: x-direction).
kx = 2.0 * pi / Lx -- Perturbed wave number (x-direction); kx * Vx_drift_elc is close to the electron plasma frequency.
vx_max_elc = 6.0 * Vx_drift_elc -- Domain boundary (velocity space: vx-direction).
poly_order = 1 -- Polynomial order.
basis_type = "tensor" -- Basis function set (hybrid phase-space basis of the kinetic species).
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 80.0 -- Final simulation time.
num_frames = 2 -- Number of output frames.
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

  -- Electrons (kinetic).
  elc = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_elc, mass = mass_elc,

    -- Velocity space grid.
    lower = { -vx_max_elc },
    upper = { vx_max_elc },
    cells = { Nvx },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]
          return n0 * (1.0 + alpha * math.cos(kx * x)) -- Electron total number density.
        end,
        temperatureInit = function (t, xn)
          return Te -- Electron isotropic temperature.
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift_elc -- Electron drift velocity.
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_ee -- Self-collision frequency.
      end
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Ions (adiabatic fluid, HLLC fluxes by default).
  ion = Vlasov.FluidSpecies.new {
    charge = charge_ion, mass = mass_ion,
    equation = Euler.new { gasGamma = gas_gamma },

    init = function (t, xn)
      -- Uniform ions at rest: mass density, momentum density and total energy density.
      return mass_ion * n0, 0.0, 0.0, 0.0, n0 * Ti / (gas_gamma - 1.0)
    end,
  },

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,
    useGhostCurrent = true,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local Ex = -alpha * math.sin(kx * x) / kx -- Total electric field (x-direction).
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
