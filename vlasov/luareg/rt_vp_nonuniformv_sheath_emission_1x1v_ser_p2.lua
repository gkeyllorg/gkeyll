-- Sheath formation at an electron-emitting wall with the Vlasov-Poisson system of equations, run
-- to a steady state. The setup of the sheath test (symmetry boundary at x = 0, conducting wall
-- at x = Lx, a boundary-flux source replacing the particles lost to the wall, BGK collisions
-- holding the pre-sheath at its source temperature, the position map clustering the cells at
-- the wall) with the real ion mass and secondary electron emission from the wall: the energy-
-- dependent emission model of Bradshaw and Srinivasan (2024, PSST 33, 035008) with the copper
-- fits of Furman and Pivi (true secondaries with a Gaussian spectrum, elastic backscattering)
-- scaled to Te = 10 eV. A quadratic velocity map clusters the electron cells at v = 0.
-- Figures of merit at steady state (serendipity p2, 16 x 24 cells, t = 4000/omega_pe): the flux-
-- weighted true-secondary yield is 0.190 and the total yield (emitted over incoming electron flux,
-- from the net electron wall flux) 0.62, below the Hobbs-Wesson critical yield 1 - 8.3 sqrt(me/mi)
-- = 0.81, so the sheath stays monotonic; the midplane-to-wall potential drop e*dphi/Te = 2.27
-- against 3.23 without emission, a reduction of 0.96 where the Hobbs-Wesson floating-potential
-- shift -ln(1 - delta) gives 0.95; ion wall flux / (n_mid c_s) = 0.72, ion speed at the wall / c_s
-- = 2.33; the net electron and ion wall fluxes balance to 1e-4. The emitted distribution and the
-- yield diagnostic (<species>_bc_up_yield.gkyl: yield and incoming flux) are written at every
-- frame.

local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 1836.153 -- Ion mass.
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

-- Secondary electron emission from the wall: the copper fits of Furman and Pivi (2002) with
-- every energy expressed in units of the electron temperature (the fits are in eV; here
-- Te = 10 eV). The emitted spectrum is the Gaussian in log-energy, the true-secondary yield
-- is the Furman-Pivi fit and the elastic backscattering is the Furman-Pivi low-energy fit.
Te_eV = 10.0 -- Electron temperature the emission fits are scaled with (eV).
E_0 = 1.97 / Te_eV -- Spectrum: peak energy (units of Te).
tau = 0.88 -- Spectrum: width in log-energy.
deltahat_ts = 1.885 -- Yield: peak true-secondary yield.
Ehat_ts = 276.8 / Te_eV -- Yield: energy of the peak yield (units of Te).
t1, t2, t3, t4 = 0.66, 0.8, 0.7, 1.0 -- Yield: angular-dependence parameters (unused in 1V).
s = 1.54 -- Yield: shape parameter.
P1_inf = 0.02 -- Elastic: backscattering probability at high energy.
P1_hat = 0.496 -- Elastic: backscattering probability at zero energy.
E_hat = 1.0e-6 / Te_eV -- Elastic: energy of the peak backscattering (units of Te).
W = 60.86 / Te_eV -- Elastic: decay energy of the backscattering (units of Te).
p = 1.0 -- Elastic: shape parameter.

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Nvx = 24 -- Cell count (velocity space: vx-direction).
Lx = 64.0 * lambda_D -- Domain size (configuration space: x-direction).
Ls = 40.0 * lambda_D -- Domain size (source).
L_nu = 32.0 * lambda_D -- Extent of the collisional (fixed-temperature) region.
vx_max_elc = 5.0 * vte -- Domain boundary (electron velocity space: vx-direction).
vx_max_ion = 8.0 * vti -- Domain boundary (ion velocity space: vx-direction).
vx_lin_elc = 0.3 -- Electron velocity map: cell size at the origin relative to a uniform grid.
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 4000.0 / omega_pe -- Final simulation time.
t_bound = 500.0 / omega_pe -- Time over which the emission is ramped up from zero.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Secondary electron emission at the wall, driven by the electrons hitting it. The models
-- take the unit charge that converts m v^2 / 2 to the energy unit of the fits; with the fits
-- scaled to Te and Te = 1 in code units, that charge is 1.
emission = G0.Vlasov.Emission.new {
  tBound = t_bound,
  inSpecies = { "elc" },
  spectrum = { G0.Vlasov.EmissionSpectrum.Gaussian.new { charge = 1.0, E0 = E_0, tau = tau } },
  yield = {
    G0.Vlasov.EmissionYield.FurmanPivi.new {
      charge = 1.0, deltaHat = deltahat_ts, EHat = Ehat_ts, t1 = t1, t2 = t2, t3 = t3, t4 = t4, s = s
    }
  },
  elastic = G0.Vlasov.EmissionElastic.FurmanPivi.new {
    charge = 1.0, P1Inf = P1_inf, P1Hat = P1_hat, EHat = E_hat, W = W, p = p
  },
}

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
    lower = { -1.0 },
    upper = { 1.0 },
    cells = { Nvx },

    -- Quadratic electron velocity map: finest cells at the origin (where the emitted electrons
    -- sit), stretching to the domain boundary.
    mapc2pVel = {
      {
        vmap = function (t, xn)
          local vc = xn[1]
          local vx_lin = vx_lin_elc * vx_max_elc
          return vx_lin * vc + (vx_max_elc - vx_lin) * vc * math.abs(vc)
        end
      },
    },

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
        type = G0.SpeciesBc.bcEmission,
        emission = emission
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

  isElectrostatic = true,

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = epsilon0,

    poissonBcs = {
      lowerType = {
        G0.PoissonBc.bcNeumann
      },
      upperType = {
        G0.PoissonBc.bcDirichlet
      },
      lowerValue = {
        0.0
      },
      upperValue = {
        0.0
      }
    }
  }
}

vlasovApp:run()