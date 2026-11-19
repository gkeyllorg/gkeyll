-- Isotropization of the temperature of a bi-Maxwellian under self-collisions with the
-- Rosenbluth/Fokker-Planck operator (FPO, Vlasov, 1x3v), with the Rosenbluth potentials of a
-- Maxwellian (J. Rodman et al., J. Comput. Phys. 2026). A homogeneous species of unit density
-- with thermal speed 1 along v_x and v_z and 1.3 along v_y relaxes for one isotropization time,
-- 1/nu_T, with nu_T = Gamma n/(2 sqrt(pi) v_t^3) A^-2 (-3 + (A + 3) a_F), A = T_perp/T_par - 1 and
-- a_F = atanh(sqrt(-A))/sqrt(-A) for A < 0 (NRL Plasma Formulary), where par is the y direction.
-- Figures of merit (serendipity p2, 2 x 10 x 10 x 10 cells): density, momentum and energy are
-- conserved to 6e-14 (largest change of the integrated M0, M1 and M2 over the run, relative to M0,
-- M0 v_rms and M2). From the M2ij moments of frames 0 and 1, T_x = T_z rises from 1.000 to 1.2165
-- and T_y falls from 1.690 to 1.2569, so the anisotropy T_y - T_x falls to 0.0585 of its initial
-- value. Integrating dT_perp/dt = -(1/2) dT_par/dt = -nu_T (T_perp - T_par) with the evolving
-- nu_T, which assumes the distribution stays bi-Maxwellian, gives 1.2195, 1.2509 and 0.0455.
local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
mass = 1.0 -- Species mass.
charge = 0.0 -- Species charge.

n0 = 1.0 -- Reference number density.
vt_perp = 1.0 -- Thermal velocity along v_x and v_z.
vt_par = 1.3 -- Thermal velocity along v_y.
gamma = 1.0 -- Prefactor of the FPO.

-- Temperature isotropization rate (NRL Plasma Formulary), T_perp < T_par.
A = vt_perp * vt_perp / (vt_par * vt_par) - 1.0
aF = 0.5 * math.log((1.0 + math.sqrt(-A)) / (1.0 - math.sqrt(-A))) / math.sqrt(-A)
nu_T = gamma * n0 / (2.0 * math.sqrt(math.pi) * math.pow(vt_par, 3.0)) * (-3.0 + (A + 3.0) * aF) / (A * A)

-- Simulation parameters.
Nx = 2 -- Cell count (configuration space: x-direction).
Nvx = 10 -- Cell count (velocity space: vx-direction).
Nvy = 10 -- Cell count (velocity space: vy-direction).
Nvz = 10 -- Cell count (velocity space: vz-direction).
Lx = 5.0 -- Domain size (configuration space: x-direction).
vx_max = 6.0 -- Domain boundary (velocity space: vx-direction).
vy_max = 6.0 -- Domain boundary (velocity space: vy-direction).
vz_max = 6.0 -- Domain boundary (velocity space: vz-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 1.0 / nu_T -- Final simulation time.
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

  -- Bi-Maxwellian species.
  bimax = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge, mass = mass,

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
          local vx = xn[2]
          local vy = xn[3]
          local vz = xn[4]

          -- Bi-Maxwellian, hotter along v_y.
          local f = n0 / (math.pow(2.0 * math.pi, 3.0 / 2.0) * vt_perp * vt_perp * vt_par) *
            math.exp(-(vx * vx + vz * vz) / (2.0 * vt_perp * vt_perp) - vy * vy / (2.0 * vt_par * vt_par))

          return f
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.FPO,

      selfNu = function (t, xn)
        return gamma -- Prefactor of the FPO.
      end
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2ij }
  },

  skipField = true,
}

vlasovApp:run()
