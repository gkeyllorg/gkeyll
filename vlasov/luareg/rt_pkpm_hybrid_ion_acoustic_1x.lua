-- Ion acoustic wave with PKPM ions and fluid electrons (1x), the counterpart of
-- vlasov/creg/rt_vlasov_hybrid_ion_acoustic_1x1v.c (kinetic ions in one velocity dimension with the
-- same electron fluid and the same ion self-collisions; collisions between a kinetic and a fluid
-- species are not supported, so that test has no ion-electron exchange).
-- Maxwellian ions (mi/me = 100, Te/Ti = 2) and an adiabatic electron fluid (gamma = 5/3, HLLC flux)
-- carry a travelling wave at k lambda_De = 0.05, started from the linear densities and velocities
-- of the wave with an ion density amplitude of 1e-4. The ions are a PKPM species: the density and
-- the velocity of the wave are in their fluid state and F is a Maxwellian in w = v_par/xi, which is
-- not the kinetic eigenfunction: the rest of the perturbation phase-mixes during the first period,
-- so the figures of merit are taken after it. A uniform magnetic field along x sets the field
-- direction; along it the ions are kinetic, so the collisionless linear theory is that of the
-- kinetic test. The ions have LBO collisions with constant frequencies: nu_ii = 2e-4 omega_pi among
-- themselves and nu_ie = (me/mi) nu_ei with the electron fluid, nu_ei = 2e-3 omega_pi. The pair
-- exchanges momentum and energy through the fluid states, which relaxes Te - Ti by 1.6% during the
-- run (Te/Ti = 1.98 at the end).
-- The fluid-field coupling is the staggered operator split with the implicit fluid-Maxwell source
-- step, which solves the momenta of both species with E;
-- rt_pkpm_hybrid_ion_acoustic_explicit_1x.lua is the same problem with the coupling advanced
-- explicitly within the SSP-RK3 step.
-- Collisionless linear theory at Te/Ti = 2, 1 + chi_i + chi_e = 0 with chi_i = (1 + zeta Z(zeta)) /
-- (k lambda_Di)^2, zeta = omega / (sqrt(2) k v_ti), and chi_e = -omega_pe^2 / (omega^2 - gamma k^2
-- Te / me): omega = 0.09470 omega_pi, damping rate 0.01333 omega_pi.
-- Figures of merit (32 cells, 48 in w, 1962 steps, omega_pe dt = 1.014 rising to 1.016 with the
-- falling sound speed of the electron fluid): between t = 1 and 3 wave periods the ion density wave
-- decays at 0.01324 omega_pi (-0.6%) with omega = 0.09439 omega_pi (-0.3%). The kinetic test, with
-- its one velocity dimension, gives 0.01307 and 0.09448, so its figures are not the values to
-- reproduce here. The total energy changes by -5e-13 of itself (7e-5 of twice the initial kinetic
-- energy of the ion flow).

local Vlasov = G0.Vlasov
local Euler = G0.Vlasov.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 100.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.
Te = 1.0 -- Electron temperature.
Ti = 0.5 -- Ion temperature.
gas_gamma = 5.0 / 3.0 -- Adiabatic index of the electron fluid.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.
B0 = 1.0 -- Magnetic field along x (sets the field direction b).
omega_pi = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_ion)) -- Ion plasma frequency.
nu_ii = 2.0e-4 * omega_pi -- Ion-ion collision frequency.
nu_ei = 2.0e-3 * omega_pi -- Electron-ion collision frequency.
nu_ie = (mass_elc / mass_ion) * nu_ei -- Ion-electron collision frequency.

alpha = 1.0e-4 -- Amplitude of the ion density of the wave, relative to n0.
k0 = 0.05 -- Wave number (k lambda_De = 0.05).
-- Root of the dispersion relation of kinetic ions and adiabatic fluid electrons (see the header),
-- in units of the electron plasma frequency.
omega_r = 0.009469560873 -- Real frequency of the wave.
omega_i = -0.001332945847 -- Growth rate of the wave.

-- Densities of the wave per unit electric field (real and imaginary parts), all perturbations
-- proportional to exp(i (k x - omega t)). Electrons (continuity and momentum of an adiabatic fluid):
-- n_e = -i q_e n0 k E / (gamma Te k^2 - m_e omega^2). Ions (Gauss's law):
-- i k eps0 E = q_i n_i + q_e n_e.
local den_r = gas_gamma * Te * k0 * k0 - mass_elc * (omega_r * omega_r - omega_i * omega_i)
local den_i = -2.0 * mass_elc * omega_r * omega_i
local den_sq = den_r * den_r + den_i * den_i
elc_resp_r = -charge_elc * n0 * k0 * den_i / den_sq
elc_resp_i = -charge_elc * n0 * k0 * den_r / den_sq
ion_resp_r = -charge_elc * elc_resp_r / charge_ion
ion_resp_i = (k0 * epsilon0 - charge_elc * elc_resp_i) / charge_ion
E0 = alpha * n0 / math.sqrt(ion_resp_r * ion_resp_r + ion_resp_i * ion_resp_i) -- Field amplitude.

-- Simulation parameters.
Nx = 32 -- Cell count (configuration space: x-direction).
Nvx = 48 -- Cell count (ion velocity space: normalized parallel velocity w).
Lx = 2.0 * pi / k0 -- Domain size (configuration space: x-direction).
vx_max = 6.0 -- Domain boundary (ion velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 3.0 * (2.0 * pi / omega_r) -- Final simulation time (three wave periods).
num_frames = 3 -- Number of output frames (one per wave period).
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

  -- Ions (PKPM).
  ion = Vlasov.PKPMSpecies.new {
    charge = charge_ion, mass = mass_ion,

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
      },

      collisions = {
        collisionID = G0.Collisions.LBO,

        selfNu = function (t, xn)
          return nu_ii
        end,

        numCrossCollisions = 1,
        collideWith = { "elc" },
        collideWithCrossNu = {
          {
            crossNu = function (t, xn)
              return nu_ie
            end
          }
        }
      }
    },

    fluid = {
      -- Initial conditions.
      init = function (t, xn)
        local x = xn[1]

        -- Density of the wave and the velocity that moves it at the phase speed (continuity:
        -- u = omega_r n / (k n0)).
        local n1 = E0 * (ion_resp_r * math.cos(k0 * x) - ion_resp_i * math.sin(k0 * x))
        local n = n0 + n1
        local ux = omega_r * n1 / (k0 * n0)

        local rho = mass_ion * n -- Ion mass density.
        local E = 0.5 * rho * ux * ux + 1.5 * rho * vti * vti -- Ion total energy density.

        return rho, rho * ux, 0.0, 0.0, E
      end
    },

    evolve = true -- Evolve species?
  },

  -- Electrons (adiabatic fluid): they collide with the ions through the exchange.
  elc = Vlasov.FluidSpecies.new {
    charge = charge_elc, mass = mass_elc,
    equation = Euler.new { gasGamma = gas_gamma },

    init = function (t, xn)
      local x = xn[1]

      -- Density, velocity (continuity: u = omega n / (k n0)) and pressure (adiabatic) of the wave.
      local cos_kx, sin_kx = math.cos(k0 * x), math.sin(k0 * x)
      local n1 = E0 * (elc_resp_r * cos_kx - elc_resp_i * sin_kx)
      local ux = E0 / (k0 * n0) * ((omega_r * elc_resp_r - omega_i * elc_resp_i) * cos_kx -
        (omega_r * elc_resp_i + omega_i * elc_resp_r) * sin_kx)
      local rho = mass_elc * (n0 + n1)
      local p = n0 * Te + gas_gamma * Te * n1

      return rho, rho * ux, 0.0, 0.0, p / (gas_gamma - 1.0) + 0.5 * rho * ux * ux
    end,

    collisions = {
      numCrossCollisions = 1,
      collideWith = { "ion" },
      collideWithCrossNu = {
        {
          crossNu = function (t, xn)
            return nu_ei
          end
        }
      }
    }
  },

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      return E0 * math.cos(k0 * x), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    externalFieldInit = function (t, xn)
      return 0.0, 0.0, 0.0, B0, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
  }
}

vlasovApp:run()
