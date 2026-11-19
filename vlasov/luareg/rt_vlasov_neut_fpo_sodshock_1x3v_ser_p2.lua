-- Sod shock tube for a neutral gas with self-collisions with the Rosenbluth/Fokker-Planck
-- operator (FPO, Vlasov, 1x3v), on a uniform velocity grid.
-- Density 1 and pressure 1 on the left, density 1/8 and pressure 1/10 on the right, at rest. The
-- prefactor of the operator, Gamma, is uniform and set from the collision frequency nu of the left
-- state through nu = n Gamma/(3 sqrt(pi) v_t^3) (Braginskii's ion collision time), so the
-- collision frequency follows the local density and temperature: it is n/T^(3/2) = 0.175 times
-- smaller in the right state.
-- Exact solution of the Euler equations (adiabatic index 5/3) at t = 0.1: shock at x = 0.684,
-- contact at 0.584, post-shock density 0.230, velocity 0.841.
-- The collision frequency of the left state is 10, a mean free path of 3.2 cells there (the
-- explicit time step of the diffusion term scales as 1/nu, which makes more collisional cases
-- expensive). The LBO and BGK tests rt_vlasov_neut_lbo_sodshock_1x3v and
-- rt_vlasov_neut_bgk_sodshock_1x3v use the same grid, initial condition and nu.
-- Figures of merit (serendipity p2, 32 x 6 x 6 x 6 cells): at t = 0.1 the density n of frame 1
-- follows the exact solution with the discontinuities smoothed over several mean free paths:
-- n = 0.403 at x = 0.547 (exact plateau behind the rarefaction 0.480) and 0.196 at x = 0.641
-- (exact post-shock density 0.230); the density L1 error int |n - n_exact| dx is 0.022. The
-- integrated density and energy change by less than 2e-6 of their initial values, and the
-- integrated x-momentum grows to 0.0900 = (p_l - p_r) t.
local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
mass_neut = 1.0 -- Neutral mass.
charge_neut = 0.0 -- Neutral charge.

nl = 1.0 -- Left number density.
Tl = 1.0 -- Left temperature.
nr = 0.125 -- Right number density.
Tr = 0.8 -- Right temperature.
Vx_drift = 0.0 -- Drift velocity (x-direction).
nu = 10.0 -- Collision frequency of the left state.

-- Derived physical quantities (using normalized code units).
vt = math.sqrt(Tl / mass_neut) -- Thermal velocity (left).
-- Prefactor of the FPO, Gamma = 3 sqrt(pi) v_t^3 nu/n, from the left state.
gamma = 3.0 * math.sqrt(math.pi) * math.pow(vt, 3.0) * nu / nl

-- Simulation parameters.
Nx = 32 -- Cell count (configuration space: x-direction).
Nvx = 6 -- Cell count (velocity space: vx-direction).
Nvy = 6 -- Cell count (velocity space: vy-direction).
Nvz = 6 -- Cell count (velocity space: vz-direction).
Lx = 1.0 -- Domain size (configuration space: x-direction).
vx_max = 6.0 * vt -- Domain boundary (velocity space: vx-direction).
vy_max = 6.0 * vt -- Domain boundary (velocity space: vy-direction).
vz_max = 6.0 * vt -- Domain boundary (velocity space: vz-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.1 -- Final simulation time.
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

  -- Neutrals.
  neut = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_neut, mass = mass_neut,

    -- Velocity space grid.
    lower = { -vx_max, -vy_max, -vz_max },
    upper = { vx_max, vy_max, vz_max },
    cells = { Nvx, Nvy, Nvz },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]

          local n = 0.0
          if x < 0.5 then
            n = nl -- Total number density (left).
          else
            n = nr -- Total number density (right).
          end

          return n
        end,
        temperatureInit = function (t, xn)
          local x = xn[1]

          local T = 0.0
          if x < 0.5 then
            T = Tl -- Isotropic temperature (left).
          else
            T = Tr -- Isotropic temperature (right).
          end

          return T
        end,
        driftVelocityInit = function (t, xn)
          return Vx_drift, 0.0, 0.0 -- Total drift velocity.
        end,

        correctAllMoments = true,
        useLastConverged = true
      }
    },

    collisions = {
      collisionID = G0.Collisions.FPO,

      selfNu = function (t, xn)
        return gamma -- Prefactor of the FPO.
      end
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.LTEMoments }
  },

  skipField = true,
}

vlasovApp:run()
