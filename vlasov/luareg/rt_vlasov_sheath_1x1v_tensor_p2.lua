-- Sheath formation with the Vlasov-Maxwell system of equations, run to a steady state.
-- Half domain: symmetry boundary at x = 0, absorbing conducting wall at x = Lx. A boundary-flux
-- source replaces the particles lost to the wall, and BGK collisions relaxing to a fixed
-- temperature act upstream (|x| < L_nu) to hold the pre-sheath plasma at its source temperature;
-- the region next to the wall is collisionless. The position map clusters the cells at the wall.
-- Figures of merit at steady state (serendipity p2): midplane-to-wall potential drop e*dphi/Te =
-- 1.93, ion wall flux / (n_mid c_s) = 0.51, ion speed at the wall / c_s = 1.64.
local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 / 25.0 -- Permeability of free space (speed of light = 5 vte, the electron velocity-grid edge).
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 100.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.
Vx_drift_elc = 0.0 -- Electron drift velocity (x-direction).
Vx_drift_ion = 0.0 -- Ion drift velocity (x-direction).

-- Derived physical quantities (using normalized code units).
Te = 1.0 * charge_ion -- Electron temperature.
Ti = 1.0 * charge_ion -- Ion temperature.

vte = math.sqrt(Te / mass_elc) -- Electron thermal velocity.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.

lambda_D = math.sqrt(epsilon0 * Te / (n0 * charge_ion * charge_ion)) -- Electron Debye length.
omega_pe = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_elc)) -- Electron plasma frequency.

nu_ee = vte / (5.0 * lambda_D) -- Electron-electron collision frequency (upstream).
nu_ii = vti / (5.0 * lambda_D) -- Ion-ion collision frequency (upstream).

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Nvx = 24 -- Cell count (velocity space: vx-direction).
Lx = 64.0 * lambda_D -- Domain size (configuration space: x-direction).
Ls = 40.0 * lambda_D -- Domain size (source).
L_nu = 32.0 * lambda_D -- Extent of the collisional (fixed-temperature) region.
vx_max_elc = 5.0 * vte -- Domain boundary (electron velocity space: vx-direction).
vx_max_ion = 8.0 * vti -- Domain boundary (ion velocity space: vx-direction).
poly_order = 2 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 1500.0 / omega_pe -- Final simulation time.
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

  -- Position map clustering the cells at the wall (x = Lx): the cell size goes from
  -- (1 + A) * Lx / Nx at the symmetry boundary to (1 - A) * Lx / Nx at the wall.
  mapc2pPos = {
    {
      pmap = function (t, xn)
        local xc = xn[1]
        local A = 0.8
        return (1.0 + A) * xc - A * (xc * xc) / Lx
      end
    },
  },
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
    cells = { Nvx },

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
          return Vx_drift_elc -- Electron drift velocity.
        end
      }
    },

    source = {
      sourceID = G0.Source.BoundaryFlux,
      sourceLength = Ls,
      sourceSpecies = "ion",

      numSources = 1,
      projections = {
        {
          projectionID = G0.Projection.LTE,

          densityInit = function (t, xn)
            local x = xn[1]

            local n = 0.0

            if x < Ls then
              n = 2.0 * (Ls - x) / Ls -- Electron source total number density (inside the source region).
            else
              n = 0.0 -- Electron source total number density (outside the source region).
            end

            return n
          end,
          temperatureInit = function (t, xn)
            return Te -- Electron source isotropic temperature.
          end,
          driftVelocityInit = function (t, xn)
            return Vx_drift_elc -- Electron source drift velocity.
          end
        }
      }
    },

    collisions = {
      collisionID = G0.Collisions.BGK,

      selfNu = function (t, xn)
        local x = xn[1]

        local nu = nu_ee / (1.0 + math.exp((x - L_nu) / (6.0 * lambda_D))) -- Electron collision frequency.

        return nu
      end,

      fixedTempRelax = true
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcReflect
      },
      upper = {
        type = G0.SpeciesBc.bcAbsorb
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
    cells = { Nvx },

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
          return Vx_drift_ion -- Ion drift velocity.
        end
      }
    },

    source = {
      sourceID = G0.Source.BoundaryFlux,
      sourceLength = Ls,
      sourceSpecies = "ion",

      numSources = 1,
      projections = {
        {
          projectionID = G0.Projection.LTE,

          densityInit = function (t, xn)
            local x = xn[1]

            local n = 0.0

            if x < Ls then
              n = 2.0 * (Ls - x) / Ls -- Ion source total number density (inside the source region).
            else
              n = 0.0 -- Ion source total number density (outside the source region).
            end

            return n
          end,
          temperatureInit = function (t, xn)
            return Ti -- Ion source isotropic temperature.
          end,
          driftVelocityInit = function (t, xn)
            return Vx_drift_ion -- Ion source drift velocity.
          end
        }
      }
    },

    collisions = {
      collisionID = G0.Collisions.BGK,

      selfNu = function (t, xn)
        local x = xn[1]

        local nu = nu_ii / (1.0 + math.exp((x - L_nu) / (6.0 * lambda_D))) -- Ion collision frequency.

        return nu
      end,

      fixedTempRelax = true
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcReflect
      },
      upper = {
        type = G0.SpeciesBc.bcAbsorb
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = 0.0 -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    bcx = {
      lower = {
        type = G0.FieldBc.bcSymWall
      },
      upper = {
        type = G0.FieldBc.bcWall
      }
    },

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()