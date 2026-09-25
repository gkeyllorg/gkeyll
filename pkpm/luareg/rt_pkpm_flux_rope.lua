local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants in SI units.
epsilon0_phys = 8.8541878128e-12 -- Permittivity of free space.
epsilon0 = 1000 * epsilon0_phys -- Artificially increased permittivity of free space.
mu0 = 4.0 * pi * 1.0e-7 -- Permeability of free space.
eV = 1.602176634e-19
charge_elc = -1.602176634e-19 -- Electron charge.
charge_ion = 1.602176634e-19 -- Ion charge.
mass_p = 1.67262192369e-27 -- Proton mass.
mass_ion = 4.0 * mass_p -- Ion mass.
mass_elc = mass_ion / 400.0 -- Electron mass.

-- Plasma parameters.
nb0 = 2.0e18 -- Peak background density.
nr0 = 1.0e18 -- Peak rope density.

Te_b = 5.0 * eV -- Electron background temperature.
Te_r = 17.3 * eV -- Electron rope temperature.
Ti_b = 1.0 * eV -- Ion background temperature.
Ti_r = 3.4 * eV -- Ion rope temperature.

lambda_elc = 10
lambda_ion = 10
nu_coeff = 0.1 -- Artificial collision scaling.
nu_elc = nu_coeff * (4.0 * math.sqrt(2.0 * pi) * nb0 * charge_elc^4 * lambda_elc) / (3.0 * (4.0 * pi * epsilon0_phys)^2 * math.sqrt(mass_elc) * Te_b^1.5) -- Electron collision frequency.
nu_ion = nu_coeff * (4.0 * math.sqrt(2.0 * pi) * nb0 * charge_ion^4 * lambda_ion) / (3.0 * (4.0 * pi * epsilon0_phys)^2 * math.sqrt(mass_ion) * Ti_b^1.5) -- Ion collison frequency.

B0 = 330.0e-4 -- Background magnetic field strength (Tesla).
I0 = -750.0 -- Total current. 

rs = 4.45e-2 -- Characteristic flux rope radius.
Lr = 30.0e-2 -- Radius of background plasma column.
rope_sep = 22.1e-2 -- Flux rope separation. 
a = 0.1
alpha = 0.1

-- Rope centers.
x1 = -0.5 * rope_sep -- x-coordinate of flux rope 1 center.
y1 = 0.0 -- y-coordinate of flux rope 1 center.

x2 = 0.5 * rope_sep -- x-coordinate of flux rope 2 center.
y2 = 0.0 -- y-coordinate of flux rope 2 center.

-- Derived quantities.
Omega_ci = charge_ion * B0 / mass_ion -- Ion cyclotron frequency.
rho_s = math.sqrt(Te_b / mass_ion) / Omega_ci -- Ion sound gyroradius.
V0e = -I0 / (2.0 * pi * charge_ion * nr0 * rs * rs) --  lab-frame fluid velocity of the rope population.
nu_hyp = 1.0e-3 -- Diffusion dimensionless scale quantity.
omega_pe = math.sqrt(nb0 * charge_elc * charge_elc / (epsilon0 * mass_elc)) -- Electron plasma frequency.
c_sim = 1.0 / math.sqrt(mu0 * epsilon0) -- Artificial speed of light.
de = c_sim / omega_pe -- Electron skin depth
D_hyp = nu_hyp * (omega_pe / (2.0 * pi)) * de ^ 4 -- Hyperdiffusion coefficient.

-- Simulation parameters.
Nx = 8 -- Cell count (configuration space: x-direction).
Ny = 8 -- Cell count (configuration space: y-direction).
Nz = 16 -- Cell count (configuration space: z-direction).
Nvpar = 8 -- Cell count (velocity space: v_parallel-direction).
Lx = 0.8 -- Domain size (configuration space: x-direction).
Ly = 0.8 -- Domain size (configuration space: y-direction).
Lz = 10.0 -- Domain size (configuration space: z-direction).
vpar_max_elec = 6.0 * math.sqrt(Te_b / mass_elc) -- Domain boundary (velocity space: vpar-direction).
vpar_max_ion = 6.0 * math.sqrt(Ti_b / mass_ion) -- Domain boundary (velocity space: vpar-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.1 / Omega_ci -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Helper functions.
-- Density of the background.
function background_density(x, y)
    local r = math.sqrt(x*x + y*y)

    return 0.5 * nb0 * (1.0 - math.tanh((r - Lr) / rho_s)) + a * nb0
end

-- Density of the two flux ropes.
function rope_density(x, y)
    local dx1 = x - x1
    local dy1 = y - y1

    local dx2 = x - x2
    local dy2 = y - y2

    local r1_sq = dx1*dx1 + dy1*dy1
    local r2_sq = dx2*dx2 + dy2*dy2

    return nr0 * (math.exp(-r1_sq / (rs*rs)) + math.exp(-r2_sq / (rs*rs)))
end

-- Magnetic field from one rope.
function single_rope_Bxy(x, y, xc, yc, I_rope)
    local dx = x - xc
    local dy = y - yc

    local r_sq = dx*dx + dy*dy
    if r_sq < 1.0e-30 then  -- Safety to avoid division by zero.
        return 0.0, 0.0
    end

    local r = math.sqrt(dx*dx + dy*dy)

    local Bphi = (1.0 - math.exp(-r_sq / (rs*rs))) * mu0 * I_rope / (2.0 * pi * r)

    local sin_phi = dy / r
    local cos_phi = dx / r

    local Bx = -Bphi * sin_phi
    local By =  Bphi * cos_phi

    return Bx, By
end

-- Magnetic field from both ropes.
function rope_Bxy(x, y)
    local I_rope = 0.5 * I0 -- Discharge current of one rope

    local Bx1, By1 = single_rope_Bxy(x, y, x1, y1, I_rope)
    local Bx2, By2 = single_rope_Bxy(x, y, x2, y2, I_rope)

    return Bx1 + Bx2, By1 + By2
end

-- Inital magnetic field perturbation in z-direction.
function delta_Bz(x, y)
    local r = math.sqrt(x*x + y*y)

    return (mu0 * nb0 / (2.0 * B0)) * (Te_b + Ti_b) * (1.0 + math.tanh((r - Lr) / rho_s))
end

-- z-component of b-hat.
function bhat_z(x, y)
    local Bx, By = rope_Bxy(x, y)
    local Bz = B0 + delta_Bz(x, y)

    local Bmag = math.sqrt(Bx*Bx + By*By + Bz*Bz)

    return Bz / Bmag
end

-- 1D Maxwellian.
function maxwellian_1v(n, mass, T, v_shift)
    return n * math.sqrt(mass / (2.0 * pi * T)) * math.exp(-mass * v_shift * v_shift / (2.0 * T))
end



pkpmApp = PKPM.App.new{

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { -0.5 * Lx, -0.5 * Ly, 0.0 },
  upper = { 0.5 * Lx, 0.5 * Ly, Lz },
  cells = { Nx, Ny, Nz },
  cflFrac = cfl_frac,

  useExplicitSource = true,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 , 1, 1}, 

  -- Boundary condtions.
  periodicDirs = {},

  -- Electrons.
  elc = PKPM.Species.new{
    modelID = G0.Model.Default,
    charge = charge_elc, mass = mass_elc,

    -- Velocity space grid
    lower = {-vpar_max_elec},
    upper = {vpar_max_elec},
    cells = {Nvpar},
    
    -- Initial conditions (distribution function).
    initDist = function (t, xn)
      local x, y, z, vpar = xn[1], xn[2], xn[3], xn[4]
      
      local nb = background_density(x, y) -- Density of the background.
      local nr = rope_density(x, y) -- Density of the flux ropes.
      
      local bz = bhat_z(x, y) -- z-component of b-hat.
      
      local ue = nr * V0e / (nb + nr) -- Electron bulk flow.
      
      -- Parallel projections.
      local ue_par = ue * bz
      local V0e_par = V0e * bz
      
      local Fe_b = maxwellian_1v(nb, mass_elc, Te_b, vpar + ue_par) -- Background electron population written in the PKPM moving frame.
      local Fe_r = maxwellian_1v(nr, mass_elc, Te_r, vpar + ue_par - V0e_par) -- Rope electron population written in the PKPM moving frame.
      
      local F0 = Fe_b + Fe_r -- Distribution function (F0).
      local G = (Te_b * Fe_b + Te_r * Fe_r) / mass_elc -- Distribution function (G).
      
      return F0, G
    end,
    
    -- Initial conditions (fluid).
    initFluid = function (t, xn)
        local x, y, z = xn[1], xn[2], xn[3]

        local nb = background_density(x, y) -- Density of the background.
        local nr = rope_density(x,y) -- Density of the flux ropes.
        local ue = nr * V0e / (nb + nr) -- Electron bulk flow.

        local mom_x = 0.0 -- Electron total momentum density (x-direction).
        local mom_y = 0.0 -- Electron total momentum density (y-direction).
        local mom_z = (nb + nr) * ue * mass_elc-- Electron total momentum density (z-direction).

        return mom_x, mom_y, mom_z
    end,
    
    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_elc -- Electron collision frequency.
      end,
      correctAllMoments = true
    },

    diffusion = {
      D = D_hyp,
      order = 4
    },

    evolve = true, -- Evolve species? 

    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy },
    bcy = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy },
    bcz = { G0.SpeciesBc.bcLineTied, G0.SpeciesBc.bcCopy },

    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 } 
  },

  -- Ions.
  ion = PKPM.Species.new {
    modelID = G0.Model.Default,
    charge = charge_ion, mass = mass_ion,
    
    -- Velocity space grid.
    lower = { -vpar_max_ion },
    upper = { vpar_max_ion },
    cells = { Nvpar },

    -- Initial conditions (distribution function).
    initDist = function (t, xn)
      local x, y, z, vpar = xn[1], xn[2], xn[3], xn[4]
      
      local nb = background_density(x, y) -- Density of the background.
      local nr = rope_density(x, y) -- Density of the flux ropes.

      local Fi_b = maxwellian_1v(nb, mass_ion, Ti_b, vpar) -- Background ion population written in the PKPM moving frame.
      local Fi_r = maxwellian_1v(nr, mass_ion, Ti_r, vpar) -- Rope ion population written in the PKPM moving frame.

      local F0 = Fi_b + Fi_r -- Ion distribution function (F0).
      local G = (Ti_b * Fi_b + Ti_r * Fi_r) / mass_ion -- Ion distribution function (G).
      
      return F0, G
    end,

    -- Initial conditions (fluid).
    initFluid = function (t, xn)
      local mom_x = 0.0 -- Ion total momentum density (x-direction).
      local mom_y = 0.0 -- Ion total momentum density (y-direction).
      local mom_z = 0.0 -- Ion total momentum density (z-direction).

      return mom_x, mom_y, mom_z
    end,

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_ion -- Ion collision frequency.
      end,
      correctAllMoments = true
    },

    diffusion = {
      D = D_hyp,
      order = 4
    },

    evolve = true, -- Evolve species?

    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy },
    bcy = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy },
    bcz = { G0.SpeciesBc.bcLineTied, G0.SpeciesBc.bcCopy },

    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2 }
  },

  
  -- Field.
  field = PKPM.Field.new {
    epsilon0 = epsilon0,
    mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x, y, z = xn[1], xn[2], xn[3]

      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx, By = rope_Bxy(x, y) -- Total magnetic field (x and y directions).
      local Bz = delta_Bz(x, y) -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    externalFieldInit = function (t, xn)
      local Ex = 0.0 -- External electric field (x-direction).
      local Ey = 0.0 -- External electric field (y-direction).
      local Ez = 0.0 -- External electric field (z-direction).

      local Bx = 0.0 -- External magnetic field (x-direction).
      local By = 0.0 -- External magnetic field (y-direction).
      local Bz = B0 -- External magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz
    end,

    sigma = function (t,xn)
      local x, y, z = xn[1], xn [2], xn[3]
      local dx1 = x - x1
      local dy1 = y - y1
      local dx2 = x - x2
      local dy2 = y - y2
      local r = math.sqrt( x*x + y*y )
      local r1 = math.sqrt( dx1*dx1 + dy1*dy1 )
      local r2 = math.sqrt( dx2*dx2 + dy2*dy2 )

      -- hyperbolic tanh profile for damping EM waves
      local tau = 2.5e-10 -- decay timescale of fields in boundary region (must be greater than dt)

      -- Define sigma shaping functions noramlized to 1 (1 -> sigma is HIGH, 0 -> sigma is LOW)
      local column = 0.5*(1 + math.tanh((r - (Lr))/rho_s)) -- kill fields outside of top hat density column
      local boundary = 0.25*(1 + math.tanh((r1 - rs/math.sqrt(alpha))/rho_s))*(1 + math.tanh((r2 - rs/math.sqrt(alpha))/rho_s))*0.5*(1 - math.tanh((z - 30*de)/(10*de))) -- kill fields away from flux rope centers near z = 0

      return (1/tau)*(1 - (1 - column)*(1 - boundary)) -- continuous verison of logical OR operator
    end,

    evolve = true,

    bcx = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall },
    bcy = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall },
    bcz = { G0.FieldBc.bcPECWall, G0.FieldBc.bcCopy },

    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }


}
pkpmApp:run()