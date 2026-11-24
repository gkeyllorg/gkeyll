-- Buneman instability with PKPM electrons and a DG Euler fluid for the ions (1x), the mixed
-- counterpart of vlasov/creg/rt_pkpm_buneman_explicit_1x.c (both species PKPM) and
-- vlasov/creg/rt_vlasov_hybrid_buneman_1x1v.c (kinetic electrons with the same ion fluid).
-- Parameters follow entry JE33 of Ammar's Simulation Journal
-- (https://ammar-hakim.org/sj/je/je33/je33-buneman.html): m_i/m_e = 25, electron drift 0.159 so
-- that k v_d = omega_pe for the box mode k = 2 pi, vte = 0.02, vti = 0.001, density perturbation
-- 1e-5 of the electrons with the electric field from Gauss's law. The electrons are a PKPM species:
-- the drift is their fluid velocity and F is a Maxwellian in w = v_par/xi; a uniform magnetic field
-- along x sets the field direction. The ions are a fluid (gamma = 5/3, HLLC flux) with the
-- temperature of the kinetic test. The ghost current of each species cancels its mean current, so
-- that Ampere's law equals the periodic Poisson solve.
-- The fluid-field coupling is advanced explicitly within the SSP-RK3 step (use_explicit_sources):
-- the Lorentz force is in the right-hand side of both species and their currents are on the Ampere
-- right-hand side of every stage; rt_pkpm_hybrid_buneman_1x.lua is the same problem with the
-- staggered operator split.
-- Linear theory with kinetic electrons and the adiabatic ion fluid at k = 2 pi: gamma = 0.19764
-- omega_pe, omega_r = 0.14251 omega_pe (0.19758 and 0.14254 with kinetic ions).
-- Figures of merit (16 cells, 64 in w): the field-energy history int E_x^2 dx grows at gamma =
-- 0.1896 over t = 20-65 (both species PKPM: 0.1895; kinetic electrons with the same ion fluid:
-- 0.1946). It peaks at 3.4e-3 at t = 72.3 (both PKPM: 3.2e-3 at t = 72.3) and is 1.8e-3 at t = 80.
-- Electron and ion number are conserved to 1e-12. The total energy (field and the fluid energies of
-- both species from the histories) changes by 7e-14 of itself up to t = 40 and by 5e-13 up to t =
-- 80. The mean of E_x is at round-off (below 1e-15) at t = 80, where the largest cell average of
-- E_x is 0.083.

local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 25.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.
Vx_drift_elc = 0.159 -- Electron drift velocity (x-direction).
Vx_drift_ion = 0.0 -- Ion drift velocity (x-direction).
vte = 0.02 -- Electron thermal velocity.
vti = 0.001 -- Ion thermal velocity.
gas_gamma = 5.0 / 3.0 -- Adiabatic index of the ion fluid.
alpha = 1.0e-5 -- Applied perturbation amplitude.
B0 = 1.0 -- Magnetic field along x (sets the field direction b).

Lx = 1.0 -- Domain size (configuration space: x-direction).
kx = 2.0 * pi / Lx -- Perturbed wave number (x-direction).
omega_pe = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_elc)) -- Electron plasma frequency.

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Nvx = 64 -- Cell count (electron velocity space: normalized parallel velocity w).
vx_max_elc = 8.0 -- Domain boundary (electron velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 80.0 / omega_pe -- Final simulation time.
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

  useExplicitSources = true,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coordinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Electrons.
  elc = Vlasov.PKPMSpecies.new {
    charge = charge_elc, mass = mass_elc,

    kinetic = {
      -- Velocity space grid.
      lower = { -vx_max_elc },
      upper = { vx_max_elc },
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

        local rho = mass_elc * n0 * (1.0 + alpha * math.cos(kx * x)) -- Electron mass density.
        local ux = Vx_drift_elc -- Electron drift.
        local E = 0.5 * rho * ux * ux + 1.5 * rho * vte * vte -- Electron total energy density.

        return rho, rho * ux, 0.0, 0.0, E
      end
    },

    evolve = true -- Evolve species?
  },

  -- Ions (adiabatic fluid, HLLC fluxes by default).
  ion = Vlasov.FluidSpecies.new {
    charge = charge_ion, mass = mass_ion,
    equation = Euler.new { gasGamma = gas_gamma },

    -- Initial conditions.
    init = function (t, xn)
      local rho = mass_ion * n0 -- Ion mass density.
      local ux = Vx_drift_ion -- Ion drift.
      local E = 0.5 * rho * ux * ux + rho * vti * vti / (gas_gamma - 1.0) -- Ion total energy density.

      return rho, rho * ux, 0.0, 0.0, E
    end,
  },

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local Ex = -alpha * math.sin(kx * x) / kx -- Electric field from Gauss's law.

      return Ex, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    externalFieldInit = function (t, xn)
      return 0.0, 0.0, 0.0, B0, 0.0, 0.0
    end,

    useGhostCurrent = true,

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()
