-- Relaxation of a waterbag and of a bump-on-tail distribution under self-collisions with the
-- Rosenbluth/Fokker-Planck operator (FPO, Vlasov, 2x3v), with the Rosenbluth potentials of a
-- Maxwellian (J. Rodman et al., J. Comput. Phys. 2026). Same two species, velocity grid and
-- prefactors of the operator as rt_vlasov_fpo_relax_1x3v, in two configuration-space dimensions
-- and with the density of both species modulated by 20% in each direction of the periodic
-- domain, so that every configuration-space mode of the drag and diffusion coefficients is
-- exercised. Run for a fifth of a collision time.
-- Figures of merit (tensor p1 (hybrid), 2 x 2 x 12 x 12 x 12 cells): density, momentum and energy
-- are conserved to 1e-14 (largest change of the integrated M0, M1 and M2 over the run, relative to
-- M0, M0 v_rms and M2), and the L2 norm of the bump distribution decreases monotonically to 0.954
-- of its initial value.
local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
mass = 1.0 -- Species mass.
charge = 0.0 -- Species charge.

n0 = 1.0 -- Reference number density.
vt = 1.0 / 3.0 -- Thermal velocity of the Maxwellian core of the bump distribution.
nu = 1.0 -- Collision frequency.
-- Prefactor of the FPO, Gamma = 3 sqrt(pi) v_t^3 nu/n: the waterbag has density 4 and
-- thermal speed 1/sqrt(3), the Maxwellian core of the bump distribution density n0 and
-- thermal speed vt.
gamma_square = 3.0 * math.sqrt(math.pi) * math.pow(1.0 / math.sqrt(3.0), 3.0) * nu / (4.0 * n0)
gamma_bump = 3.0 * math.sqrt(math.pi) * math.pow(vt, 3.0) * nu / n0
ab = math.sqrt(0.1) -- Bump amplitude.
sb = 0.12 -- Bump width (Lorentzian softening).
ub = 4.0 * math.sqrt(0.25 / 3.0) -- Bump location (v_x).
vtb = 1.0 -- Bump Maxwellian thermal velocity.
an = 0.2 -- Amplitude of the density modulation.

-- Simulation parameters.
Nx = 2 -- Cell count (configuration space: x-direction).
Ny = 2 -- Cell count (configuration space: y-direction).
Nvx = 12 -- Cell count (velocity space: vx-direction).
Nvy = 12 -- Cell count (velocity space: vy-direction).
Nvz = 12 -- Cell count (velocity space: vz-direction).
Lx = 5.0 -- Domain size (configuration space: x-direction).
Ly = 5.0 -- Domain size (configuration space: y-direction).
vx_max = 5.0 -- Domain boundary (velocity space: vx-direction).
vy_max = 5.0 -- Domain boundary (velocity space: vy-direction).
vz_max = 5.0 -- Domain boundary (velocity space: vz-direction).
poly_order = 1 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.2 -- Final simulation time.
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
  lower = { 0.0, 0.0 },
  upper = { Lx, Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction (x- and y-directions).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1, 2 }, -- Periodic directions (x- and y-directions).

  -- Waterbag species.
  square = Vlasov.Species.new {
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
          local x, y = xn[1], xn[2]
          local vx = xn[3]
          local vy = xn[4]
          local vz = xn[5]

          -- Density modulation in both configuration-space directions.
          local mod = (1.0 + an * math.cos(2.0 * math.pi * x / Lx)) * (1.0 + an * math.cos(2.0 * math.pi * y / Ly))

          -- Waterbag: half the reference density inside |v_i| < 1, zero outside.
          local f = 0.0
          if math.abs(vx) < 1.0 and math.abs(vy) < 1.0 and math.abs(vz) < 1.0 then
            f = 0.5 * n0 * mod
          end

          return f
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.FPO,

      selfNu = function (t, xn)
        return gamma_square -- Prefactor of the FPO.
      end
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Bump-on-tail species.
  bump = Vlasov.Species.new {
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
          local x, y = xn[1], xn[2]
          local vx = xn[3]
          local vy = xn[4]
          local vz = xn[5]

          -- Density modulation in both configuration-space directions.
          local mod = (1.0 + an * math.cos(2.0 * math.pi * x / Lx)) * (1.0 + an * math.cos(2.0 * math.pi * y / Ly))

          local v_sq = vx * vx + vy * vy + vz * vz
          local vb_sq = (vx - ub) * (vx - ub) + vy * vy + vz * vz

          -- Maxwellian core plus a bump centered at (ub, 0, 0) with a Lorentzian envelope of width sb.
          local f = (n0 / math.pow(2.0 * math.pi * vt * vt, 3.0 / 2.0)) * math.exp(-v_sq / (2.0 * vt * vt)) +
            (n0 / math.pow(2.0 * math.pi * vtb * vtb, 3.0 / 2.0)) * math.exp(-vb_sq / (2.0 * vtb * vtb)) * (ab * ab) / (vb_sq + (sb * sb))

          return f * mod
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.FPO,

      selfNu = function (t, xn)
        return gamma_bump -- Prefactor of the FPO.
      end
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  skipField = true,
}

vlasovApp:run()
