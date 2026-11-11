-- Sheath formation with the PKPM system of equations (1x, normalized parallel velocity
-- w = v_par/xi), the counterpart of vlasov/creg/rt_vlasov_sheath_1x1v.c.
-- Half domain: symmetry boundary at x = 0 (reflecting species, symmetry-wall field),
-- absorbing conducting wall at x = Lx. A boundary-flux source replaces the particles lost
-- to the wall, and a BGK-form temperature source at the source temperature acts upstream
-- (x < L_nu) to hold the pre-sheath plasma at that temperature; the region next to
-- the wall is collisionless. The position map clusters the cells at the wall.
-- Run to t = 1000/omega_pe, when the wall loss is 0.6% below its steady value.
-- Figures of merit at t = 1000 (Vlasov 1x1v in parentheses): midplane-to-wall potential drop
-- e*dphi/Te = 1.929 (1.925), ion wall flux / (n_mid c_s) = 0.506 (0.509), ion speed at the
-- wall / c_s = 1.635 (1.629), with c_s = sqrt((Te + Ti)/mi). Over the presheath (the 13 cells
-- upstream of the three sheath cells) the cell averages of the density and the flow velocity
-- of both species are within 1.1% of Vlasov (rms 0.4%), the largest difference being in the
-- cell at the sheath edge; the wall cell has the electron density 2.2% low and the electron
-- velocity 2.0% high.

local PKPM = G0.PKPM

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 / 25.0 -- Permeability of free space (speed of light = 5 vte).
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 100.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.
Te_over_Ti = 1.0 -- Electron temperature / ion temperature.
B0 = 1.0 -- Reference magnetic field strength.

-- Derived physical quantities (using normalized code units).
vte = 1.0 -- Electron thermal velocity.
vti = vte / math.sqrt(Te_over_Ti * mass_ion) -- Ion thermal velocity.

lambda_D = math.sqrt(epsilon0 * vte * vte * mass_elc / (n0 * charge_ion * charge_ion)) -- Electron Debye length.
omega_pe = math.sqrt(n0 * charge_ion * charge_ion / (epsilon0 * mass_elc)) -- Electron plasma frequency.
nu_Q_elc = vte / (5.0 * lambda_D) -- Rate of the electron temperature source (upstream).
nu_Q_ion = vti / (5.0 * lambda_D) -- Rate of the ion temperature source (upstream).

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Nvx = 16 -- Cell count (velocity space: normalized parallel velocity w).
Lx = 64.0 * lambda_D -- Domain size (configuration space: x-direction).
Ls = 40.0 * lambda_D -- Extent of the source region.
L_nu = 32.0 * lambda_D -- Extent of the temperature-source region.
map_A = 0.8 -- Amplitude of the position map clustering the cells at the wall.
vx_max_elc = 5.0 -- Domain boundary (electron velocity space: w = v_par/xi).
vx_max_ion = 5.0 -- Domain boundary (ion velocity space: w = v_par/xi).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 1000.0 / omega_pe -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate integrated L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

pkpmApp = PKPM.App.new {

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
  -- (1 + A) Lx/Nx at the symmetry boundary to (1 - A) Lx/Nx at the wall.
  mapc2pPos = {
    {
      pmap = function (t, xn)
        local xc = xn[1]

        return (1.0 + map_A) * xc - map_A * (xc * xc) / Lx
      end
    },
  },

  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coordinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Electrons.
  elc = PKPM.Species.new {
    charge = charge_elc, mass = mass_elc,

    -- Velocity space grid.
    lower = { -vx_max_elc },
    upper = { vx_max_elc },
    cells = { Nvx },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[2] -- Normalized parallel velocity w = v_par/xi.

      -- Normalized isotropic Maxwellian: F = exp(-w^2/2)/sqrt(2 pi) (M0 = 1, M2 = 1) and
      -- G = (T_perp/T) F = F; the density and temperature enter through the fluid state.
      local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * math.pi)

      return F, F
    end,

    initFluid = function (t, xn)
      local rho = mass_elc * n0 -- Electron mass density.
      local E = 1.5 * rho * vte * vte -- Electron total energy density (no flow).

      return rho, 0.0, 0.0, 0.0, E
    end,

    -- Temperature source: rate nu_Q_elc upstream (x < L_nu), vanishing over 6 lambda_D towards
    -- the collisionless region next to the wall, at the source temperature.
    sourceBGK = {
      rateProfile = function (t, xn)
        local x = xn[1]

        return nu_Q_elc / (1.0 + math.exp((x - L_nu) / (6.0 * lambda_D)))
      end,

      tempShape = function (t, xn)
        return mass_elc * vte * vte -- Electron source temperature.
      end
    },

    -- Source profile: a density 2 (Ls - x)/Ls inside the source region (it integrates to Ls,
    -- so the boundary-flux source replaces the wall losses exactly), at the plasma
    -- temperature and without drift.
    source = {
      sourceID = G0.Source.BoundaryFlux,
      sourceLength = Ls,
      sourceSpecies = "ion",
      writeSource = true,

      numSources = 1,
      projections = {
        {
          densityInit = function (t, xn)
            local x = xn[1]
            local n = 0.0

            if x < Ls then
              n = 2.0 * (Ls - x) / Ls -- Source number density (inside the source region).
            end

            return n
          end,

          temperatureInit = function (t, xn)
            return mass_elc * vte * vte -- Electron source temperature.
          end,

          driftVelocityInit = function (t, xn)
            return 0.0, 0.0, 0.0 -- Source drift velocity.
          end
        }
      }
    },

    bcx = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcAbsorb }, -- Symmetry boundary and absorbing wall (x-direction).

    evolve = true -- Evolve species?
  },

  -- Ions.
  ion = PKPM.Species.new {
    charge = charge_ion, mass = mass_ion,

    -- Velocity space grid.
    lower = { -vx_max_ion },
    upper = { vx_max_ion },
    cells = { Nvx },

    -- Initial conditions.
    initDist = function (t, xn)
      local w = xn[2] -- Normalized parallel velocity w = v_par/xi.

      -- Normalized isotropic Maxwellian: F = exp(-w^2/2)/sqrt(2 pi) (M0 = 1, M2 = 1) and
      -- G = (T_perp/T) F = F; the density and temperature enter through the fluid state.
      local F = math.exp(-0.5 * w * w) / math.sqrt(2.0 * math.pi)

      return F, F
    end,

    initFluid = function (t, xn)
      local rho = mass_ion * n0 -- Ion mass density.
      local E = 1.5 * rho * vti * vti -- Ion total energy density (no flow).

      return rho, 0.0, 0.0, 0.0, E
    end,

    -- Temperature source: rate nu_Q_ion upstream (x < L_nu), vanishing over 6 lambda_D towards
    -- the collisionless region next to the wall, at the source temperature.
    sourceBGK = {
      rateProfile = function (t, xn)
        local x = xn[1]

        return nu_Q_ion / (1.0 + math.exp((x - L_nu) / (6.0 * lambda_D)))
      end,

      tempShape = function (t, xn)
        return mass_ion * vti * vti -- Ion source temperature.
      end
    },

    -- Source profile: the same density as the electrons' (the rate follows the ion wall
    -- losses), at the plasma temperature and without drift.
    source = {
      sourceID = G0.Source.BoundaryFlux,
      sourceLength = Ls,
      sourceSpecies = "ion",
      writeSource = true,

      numSources = 1,
      projections = {
        {
          densityInit = function (t, xn)
            local x = xn[1]
            local n = 0.0

            if x < Ls then
              n = 2.0 * (Ls - x) / Ls -- Source number density (inside the source region).
            end

            return n
          end,

          temperatureInit = function (t, xn)
            return mass_ion * vti * vti -- Ion source temperature.
          end,

          driftVelocityInit = function (t, xn)
            return 0.0, 0.0, 0.0 -- Source drift velocity.
          end
        }
      }
    },

    bcx = { G0.SpeciesBc.bcReflect, G0.SpeciesBc.bcAbsorb }, -- Symmetry boundary and absorbing wall (x-direction).

    evolve = true -- Evolve species?
  },

  -- Field.
  field = PKPM.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function (the guide field is external).
    init = function (t, xn)
      return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    -- Static external guide field along x (b = x_hat): not subject to the field boundary
    -- conditions, so the wall's image field does not reverse it.
    externalFieldInit = function (t, xn)
      return 0.0, 0.0, 0.0, B0, 0.0, 0.0
    end,
    evolveExternalField = false,

    bcx = { G0.FieldBc.bcSymWall, G0.FieldBc.bcPECWall }, -- Symmetry boundary and conducting wall (x-direction).

    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

pkpmApp:run()
