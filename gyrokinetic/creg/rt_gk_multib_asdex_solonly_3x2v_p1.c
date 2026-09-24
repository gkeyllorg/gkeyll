#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_efit.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_mpi_comm.h>
#include <gkyl_null_comm.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>

struct gk_app_ctx {
  int cdim, vdim;
  // Plasma parameters
  int num_species;
  double me, qe, mi, qi;
  // Initial conditions.
  double den_floor       ; 
  double Te_floor        ; 
  double Ti_floor        ; 
  double n0, Te0, Ti0, B0; // Reference parameters.
  // Collision parameters
  double nuFrac, nuElc, nuIon;
  // Source parameters
  int num_sources;
  bool adapt_energy_srcCORE, adapt_particle_srcCORE;
  double center_srcCORE[3], sigma_srcCORE[3];
  double energy_srcCORE, particle_srcCORE;
  double floor_srcCORE;
  bool adapt_energy_srcWALL, adapt_particle_srcWALL;
  double center_srcWALL[3], sigma_srcWALL[3];
  double energy_srcWALL, particle_srcWALL;
  double floor_srcWALL;
  double temp_recycle_srcWALL;
  // Grid parameters
  double psi_axis;
  double psi_sep;
  double Lx, Ly, Lz;
  double x_min, x_max, y_min, y_max, z_min, z_max, rho_min, rho_max;
  int num_cell_x, num_cell_y, num_cell_z, num_cell_vpar, num_cell_mu;
  int cells[GKYL_MAX_DIM], poly_order;
  double vpar_max_elc, mu_max_elc, vpar_max_ion, mu_max_ion;
  // Simulation control parameters
  double t_end, write_phase_freq;
  int num_frames, int_diag_calc_num, num_failures_max;
  double dt_failure_tol;

};

double rho_psi(double psi, double psi_axis, double psi_sep)
{
  // Normalized radial coordinate.
  return sqrt((psi-psi_axis) / (psi_sep - psi_axis));
}

double psi_rho(double rho, double psi_axis, double psi_sep)
{
  // Poloidal flux given the normalized radial coordinate.
  return pow(rho,2) * (psi_sep - psi_axis) + psi_axis;
}

// Density initial condition (like AUG exp profile)
void eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double psi_axis = app->psi_axis;
  double psi_sep = app->psi_sep;
  double den_floor = app->den_floor;

  double rho = rho_psi(x, psi_axis, psi_sep);

  double profile = (-2.17282515*pow(rho, 3.0) + 6.80801422*pow(rho, 2.0) - 7.12279432*rho + 2.48870613) * 1.0e22;

  fout[0] = fmax(profile, den_floor);
}

// Flow initial condition
void
eval_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

// Electron temperature initial conditions
void eval_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double psi_axis = app->psi_axis;
  double psi_sep = app->psi_sep;
  double Te_floor = app->Te_floor;
  double eV = GKYL_ELEMENTARY_CHARGE;

  double rho = rho_psi(x, psi_axis, psi_sep);

  double profile = (-140.97018851*rho + 156.25649838) * eV;

  // Multiply by a function that smoothly drops the temperature by 2X towards
  // the divertor plates.
  profile *= 0.5+0.5*exp(-0.25*pow(z/(M_PI/1.35),10));

  fout[0] = fmax(profile, Te_floor);
}

// Ion temperature initial conditions
void eval_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[2];

  struct gk_app_ctx *app = ctx;
  double psi_axis = app->psi_axis;
  double psi_sep = app->psi_sep;
  double Ti_floor = app->Ti_floor;
  double eV = GKYL_ELEMENTARY_CHARGE;

  double rho = rho_psi(x, psi_axis, psi_sep);

  double profile = (125.6039 * exp(-206.2917*(rho-1.0)) + 20.9689) * eV;

  // Multiply by a function that smoothly drops the temperature by 5X towards
  // the divertor plates.
  profile *= 0.5+0.5*exp(-0.25*pow(z/(M_PI/1.35),10));

  fout[0] = fmax(profile, Ti_floor);
}

// Taken from rt gk d3d 3x2c, is this the non uniform v grid mapping?
void mapc2p_vel_elc(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double vpar_max_elc = app->vpar_max_elc;
  double mu_max_elc = app->mu_max_elc;
  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_elc*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_elc*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_elc*2.0*pow(cvpar,2);
  // Quadratic map in mu.
  vp[1] = mu_max_elc*pow(cmu,2);
}

void mapc2p_vel_ion(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;
  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_ion*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_ion*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_ion*2.0*pow(cvpar,2);
  // Quadratic map in mu.
  vp[1] = mu_max_ion*pow(cmu,2);
}


struct gkyl_gk_block_geom*
create_gk_block_geom(void *ctx)
{
  struct gk_app_ctx *params = ctx;

  struct gkyl_gk_block_geom *bgeom = gkyl_gk_block_geom_new(3,3);

  /* Block layout and coordinates.

    z  
    ^  
    |

    |  +------------------+
    |  |                  | 
    |  |b2                |
    |  |lower inner sol   |
    |  |                  | 
    |  |                  |
    |  +------------------+
    |  |$                 |
    |  |$                 |
    |  |$ b1              |
    |  +$ upper sol       |
    |  |$                 |
    |  |$                 |
    |  +------------------+
    |  |                  |
    |  |                  | 
    |  |b0                |
    |  |lower outer sol   |
    |  |                  |
    0  +------------------+

       0 -------------------------------------------------------- -> x

    Edges that touch coincide are physically connected unless
    otherwise indicated by a special symbol. Edges with a special
    symbol such as o,x,%, or % are instead connected to the other
    edge with the same symbol. Edges that do not coincide with
    another edge are a physical boundary.
  */  

  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .filepath = "./gyrokinetic/data/eqdsk/asdex_33292_3.800.eqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };

  // Theta limits are actually set internally by the code.
  double theta_min = -1.0, theta_max = 1.0;

  double psi_sep  = params->psi_sep; // Psi at the separatrix.
  double psi_axis = params->psi_axis; // Psi at the magnetic axis.

  // Number of cells.
  int Npsi_sol        = 12;
  int Npsi_pf         = 12;
  int Ntheta_divertor = 2;
  int Ntheta_sol      = 8;

  int nalpha = 8;

  // Block 0: lower outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
      .lower = { params->x_min, params->y_min, params->z_min },
      .upper = { params->x_max, params->y_max, params->z_max },
      .cells = { Npsi_sol, nalpha, Ntheta_divertor},
      .cuts = {  1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO,
          .rclose = 2.5,           // Closest R to region of interest.
          .rright = 2.5,           // Closest R to outboard SOL.
          .rleft = 0.7,            // Closest R to inboard SOL.
          .rmax = 2.5,             // Largest R in machine.
          .rmin = 0.7,             // Smallest R in machine.
          .zmin = -1.3,            // Lower Z boundary.
          .zmax = 1.0,             // Upper Z boundary.
          .zmin_left = -1.0,       // Z of inboard divertor plate.
          .zmin_right = -1.0,      // Z of outboard divertor plate.
        }
      },
      
      .connections[0] = { // x-direction.
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL },
      },
      .connections[1] = { // y-direction.
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL },
      },
      .connections[2] = { // z-direction.
        { .bid = 0, .dir = 2, .edge = GKYL_PHYSICAL},
        { .bid = 1, .dir = 2, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // Block 1: mid SOL.
  gkyl_gk_block_geom_set_block(bgeom, 1, &(struct gkyl_gk_block_geom_info) {
      .lower = { params->x_min, params->y_min, params->z_min },
      .upper = { params->x_max, params->y_max, params->z_max },
      .cells = { Npsi_sol, nalpha, Ntheta_sol },
      .cuts = {  1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID,
          .rclose = 2.5,           // Closest R to region of interest.
          .rright = 2.5,           // Closest R to outboard SOL.
          .rleft = 0.7,            // Closest R to inboard SOL.
          .rmax = 2.5,             // Largest R in machine.
          .rmin = 0.7,             // Smallest R in machine.
          .zmin = -1.3,            // Lower Z boundary.
          .zmax = 1.0,             // Upper Z boundary.
          .zmin_left = -1.0,       // Z of inboard divertor plate.
          .zmin_right = -1.0,      // Z of outboard divertor plate.
        }
      },
      
      .connections[0] = { // x-direction.
        { .bid = 1, .dir = 0, .edge = GKYL_PHYSICAL},
        { .bid = 1, .dir = 0, .edge = GKYL_PHYSICAL },
      },
      .connections[1] = { // y-direction.
        { .bid = 1, .dir = 1, .edge = GKYL_PHYSICAL},
        { .bid = 1, .dir = 1, .edge = GKYL_PHYSICAL },
      },
      .connections[2] = { // z-direction.
        { .bid = 0, .dir = 2, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 2, .dir = 2, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // Block 2: lower inner SOL.
  gkyl_gk_block_geom_set_block(bgeom, 2, &(struct gkyl_gk_block_geom_info) {
      .lower = { params->x_min, params->y_min, params->z_min },
      .upper = { params->x_max, params->y_max, params->z_max },
      .cells = { Npsi_sol, nalpha, Ntheta_divertor},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP,
          .rclose = 2.5,           // Closest R to region of interest.
          .rright = 2.5,           // Closest R to outboard SOL.
          .rleft = 0.7,            // Closest R to inboard SOL.
          .rmax = 2.5,             // Largest R in machine.
          .rmin = 0.7,             // Smallest R in machine.
          .zmin = -1.3,            // Lower Z boundary.
          .zmax = 1.0,             // Upper Z boundary.
          .zmin_left = -1.0,       // Z of inboard divertor plate.
          .zmin_right = -1.0,      // Z of outboard divertor plate.
        }
      },
      
      .connections[0] = { // x-direction.
        { .bid = 2, .dir = 0, .edge = GKYL_PHYSICAL},
        { .bid = 2, .dir = 0, .edge = GKYL_PHYSICAL },
      },
      .connections[1] = { // x-direction.
        { .bid = 2, .dir = 1, .edge = GKYL_PHYSICAL},
        { .bid = 2, .dir = 1, .edge = GKYL_PHYSICAL },
      },
      .connections[2] = { // z-direction.
        { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 2, .dir = 1, .edge = GKYL_PHYSICAL},
      }
    }
  );

  return bgeom;
}


struct gk_app_ctx create_ctx(void)
{
  int cdim = 3, vdim = 2; // Dimensionality
  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0, eV = GKYL_ELEMENTARY_CHARGE;
  double proton_mass = GKYL_PROTON_MASS, electron_mass = GKYL_ELECTRON_MASS;

  // Location of the numerical equilibrium.

  // Get info from eqdsk file.
  struct gkyl_efit_inp efit_inp = {
    .filepath = "./gyrokinetic/data/eqdsk/asdex_33292_3.800.eqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };

  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);
  double R_axis = efit->rmaxis; // R of the magnetic axis.
  double Z_axis = efit->zmaxis; // Z of the magnetic axis.
  double psi_sep = efit->psisep; // psi at the separatrix.
  double psi_axis = efit->simag; // psi at the magnetic axis.
  double Rxpt = efit->Rxpt[0], Zxpt = efit->Zxpt[0];
  gkyl_efit_release(efit);

  // Species mass and charge.
  int num_species = 2;
  double me = electron_mass;
  double mi = proton_mass*2.01410177811; // Deuterium ions
  double qi =  eV; // ion charge
  double qe = -eV; // electron charge

  // Reference parameters.
  double n0  = 6.0e18;
  double Te0 = 11.0 * eV;
  double Ti0 = 29.0 * eV;

  // Parameters controlling initial conditions.
  double den_floor = 0.1*n0; // Min density in IC.
  double Te_floor = 4.0*eV; // Min Te in IC.
  double Ti_floor = 4.0*eV; // Min Ti in IC.

  double B0  = 2.5;
  double q0 = 5.0;
  double r0 = 0.5;
  double vte = sqrt(Te0/me);
  double vti = sqrt(Ti0/mi);
  double c_s = sqrt(Te0/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;

  // The radial extent gives [rhop_min, rhop_max] = [1.010, 1.046]
  double x_min = 0.336373;
  double x_max = 0.352308;
  double Lx = x_max - x_min;

  double z_min = -(M_PI-1e-10);
  double z_max =   M_PI-1e-10;
  double Lz = z_max - z_min;

  double Ly = 100*rho_s*q0/r0/3.0;
  double y_min = -Ly/2.;
  double y_max =  Ly/2.;

  double rho_min = rho_psi(x_min, psi_axis, psi_sep);
  double rho_max = rho_psi(x_max, psi_axis, psi_sep);

  // Source parameters
  int num_sources = 2;
  double P_exp = 0.235e6; // P_sol measured [W], half of 0.8MW
  double vol_frac = Ly/(2.0*M_PI); // Volume fraction of the simulation box
  double P_inj = P_exp * vol_frac / num_species;
  // Core source:
  // - The particles injection is only the one that are lost through the inner radial boundary.
  bool adapt_energy_srcCORE = true; // The source will compensate the losses in energy according to given boundaries.
  bool adapt_particle_srcCORE = true; // The source will compensate the losses in particle according to given boundaries.
  double energy_srcCORE = P_inj; // What the source must inject in energy [W]
  double particle_srcCORE = 0.0;// What the source must inject in particle [1/s]
  double center_srcCORE[3] = {x_min, 0.0, -Lz/4.0}; // This is the position of the ion source,
  double sigma_srcCORE[3] = {0.06*Lx, 0.0, Lz/8.0}; //  the electron source will be at +Lz/2.
  double floor_srcCORE = 1e-10;
  // Wall source:
  // - Injects energy only in the wall region.
  // - Reinjects particles that are absorbed by the wall.
  // - Energy is free to leave the system.
  bool adapt_energy_srcWALL = false; // The source will compensate the losses in energy according to given boundaries.
  bool adapt_particle_srcWALL = true; // The source will compensate the losses in particle according to given boundaries.
  double energy_srcWALL = 0.0; // What the source must inject in energy [W]
  double particle_srcWALL = 0.0;// What the source must inject in particle [1/s]
  double center_srcWALL[3] = {x_min, 0.0, -Lz/2.0}; // This is the position of the ion source,
  double sigma_srcWALL[3] = {0.25*Lx, 0.0, 0.1*Lz}; //  the electron source will be at +Lz/2.
  double floor_srcWALL = 1e-10;
  double temp_recycle_srcWALL = 4.0 * eV;

  // Grid parameters
  int num_cell_x = 12; 
  int num_cell_y = 12;
  int num_cell_z = 8;
  int num_cell_vpar = 8;
  int num_cell_mu = 4;
  int poly_order = 1;

  // Velocity box dimensions
  double vpar_max_elc = 6.*vte;
  double mu_max_elc   = me*pow(4*vte,2)/(2*B0);
  double vpar_max_ion = 6.*vti;
  double mu_max_ion   = mi*pow(4*vti,2)/(2*B0);

  double t_end = 1.0e-4;
  int num_frames = 100;
  double write_phase_freq = 0.01;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-3; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_app_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .psi_sep  = psi_sep ,
    .psi_axis = psi_axis,
    .x_min = x_min,  .x_max = x_max,
    .y_min = y_min,  .y_max = y_max,
    .z_min = z_min,  .z_max = z_max,
    .Lx = Lx,
    .Ly = Ly,
    .Lz = Lz,
    .rho_min = rho_min,  .rho_max = rho_max,
    .num_species = num_species,
    .me = me,  .qe = qe,
    .mi = mi,  .qi = qi,
    .den_floor        = den_floor       ,
    .Te_floor         = Te_floor        ,
    .Ti_floor         = Ti_floor        ,
    .n0 = n0,  .Te0 = Te0,  .Ti0 = Ti0,
    .B0 = B0,
    .num_sources = num_sources,
    .adapt_energy_srcCORE = adapt_energy_srcCORE,
    .adapt_particle_srcCORE = adapt_particle_srcCORE,
    .center_srcCORE = {center_srcCORE[0], center_srcCORE[1], center_srcCORE[2]},
    .sigma_srcCORE = {sigma_srcCORE[0], sigma_srcCORE[1], sigma_srcCORE[2]},
    .energy_srcCORE = energy_srcCORE,  .particle_srcCORE = particle_srcCORE,
    .floor_srcCORE = floor_srcCORE,
    .adapt_energy_srcWALL = adapt_energy_srcWALL,
    .adapt_particle_srcWALL = adapt_particle_srcWALL,
    .center_srcWALL = {center_srcWALL[0], center_srcWALL[1], center_srcWALL[2]},
    .sigma_srcWALL = {sigma_srcWALL[0], sigma_srcWALL[1], sigma_srcWALL[2]},
    .energy_srcWALL = energy_srcWALL,  .particle_srcWALL = particle_srcWALL,
    .floor_srcWALL = floor_srcWALL,
    .temp_recycle_srcWALL = temp_recycle_srcWALL,
    .num_cell_x     = num_cell_x,
    .num_cell_y     = num_cell_y,
    .num_cell_z     = num_cell_z,
    .num_cell_vpar  = num_cell_vpar,
    .num_cell_mu    = num_cell_mu,
    .cells = {num_cell_x, num_cell_y, num_cell_z, num_cell_vpar, num_cell_mu},
    .poly_order   = poly_order,
    .vpar_max_elc = vpar_max_elc,  .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion,  .mu_max_ion = mu_max_ion,
    .write_phase_freq = write_phase_freq,
    .t_end = t_end,  .num_frames = num_frames,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}



int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  struct gk_app_ctx ctx = create_ctx(); // Context for init functions.
                    
  // Construct block geometry
  struct gkyl_gk_block_geom *bgeom = create_gk_block_geom(&ctx);

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);


  // Electrons.
  // Electron core source:
  struct gkyl_gyrokinetic_projection proj_srcCORE_e = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
    .gaussian_mean = {ctx.center_srcCORE[0], ctx.center_srcCORE[1], ctx.center_srcCORE[2]},
    .gaussian_std_dev = {ctx.sigma_srcCORE[0], ctx.sigma_srcCORE[1], ctx.sigma_srcCORE[2]},
    .total_num_particles = ctx.particle_srcCORE,
    .total_kin_energy = ctx.energy_srcCORE,
    .temp_max = 5.0*ctx.Te0,
    .temp_min = 0.1*ctx.Te0,
    .f_floor = ctx.floor_srcCORE,
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcCORE_e ={
    .adapt_to_species = "elc",
    .adapt_particle = ctx.adapt_particle_srcCORE,
    .adapt_energy = ctx.adapt_energy_srcCORE,
    .num_boundaries = 1,
    .dir = {0},
    .edge = {GKYL_LOWER_EDGE},
  };
  // Electron wall source:
  struct gkyl_gyrokinetic_projection proj_srcWALL_e = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
    .gaussian_mean = {ctx.center_srcWALL[0], ctx.center_srcWALL[1], ctx.center_srcWALL[2]},
    .gaussian_std_dev = {ctx.sigma_srcWALL[0], ctx.sigma_srcWALL[1], ctx.sigma_srcWALL[2]},
    .total_num_particles = ctx.particle_srcWALL,
    .total_kin_energy = ctx.energy_srcWALL,
    .temp_max = 5.0*ctx.Te0,
    .temp_min = ctx.temp_recycle_srcWALL,
    .f_floor = ctx.floor_srcWALL,
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcWALL_e ={
    .adapt_to_species = "elc",
    .adapt_particle = ctx.adapt_particle_srcWALL,
    .adapt_energy = ctx.adapt_energy_srcWALL,
    .num_boundaries = 3,
    .dir = {0, 2, 2},
    .edge = {GKYL_UPPER_EDGE, GKYL_LOWER_EDGE, GKYL_UPPER_EDGE},
  };


  struct gkyl_gyrokinetic_projection elc_ic = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
    .ctx_density = &ctx,
    .density = eval_density,
    .ctx_upar = &ctx,
    .upar = eval_upar,
    .ctx_temp = &ctx,
    .temp = eval_temp_elc,
  };


  struct gkyl_gyrokinetic_multib_species_pb elc_blocks[3];

  elc_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb){
    .block_id = 0,
    .polarization_density = ctx.n0,
    .projection = elc_ic,
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = ctx.num_sources,
      .num_adapt_sources = ctx.num_sources,
      .projection[0] = proj_srcCORE_e,
      .adapt[0] = adapt_srcCORE_e,
      .projection[1] = proj_srcWALL_e,
      .adapt[1] = adapt_srcWALL_e,
      .diagnostics = {
        .num_diag_moments = 1,
        .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
      }
    },

  };

  elc_blocks[1] = (struct gkyl_gyrokinetic_multib_species_pb){
    .block_id = 1,
    .polarization_density = ctx.n0,
    .projection = elc_ic,
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = ctx.num_sources,
      .num_adapt_sources = ctx.num_sources,
      .projection[0] = proj_srcCORE_e,
      .adapt[0] = adapt_srcCORE_e,
      .projection[1] = proj_srcWALL_e,
      .adapt[1] = adapt_srcWALL_e,
      .diagnostics = {
        .num_diag_moments = 1,
        .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
      }
    },

  };
  elc_blocks[2] = (struct gkyl_gyrokinetic_multib_species_pb){
    .block_id = 2,
    .polarization_density = ctx.n0,
    .projection = elc_ic,
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = ctx.num_sources,
      .num_adapt_sources = ctx.num_sources,
      .projection[0] = proj_srcCORE_e,
      .adapt[0] = adapt_srcCORE_e,
      .projection[1] = proj_srcWALL_e,
      .adapt[1] = adapt_srcWALL_e,
      .diagnostics = {
        .num_diag_moments = 1,
        .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
      }
    },

  };



  struct gkyl_gyrokinetic_bc elc_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 2, .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
  };

  struct gkyl_gyrokinetic_multib_species elc = {
    .name = "elc",
    .charge = ctx.qe, .mass = ctx.me,
    .vdim = ctx.vdim,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0},
    .cells = { cells_v[0], cells_v[1] },

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm.
      .temp_ref = ctx.Te0, // Temperature used to calculate coulomb logarithm.
      .bmag_ref = ctx.B0,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },
    .num_physical_bcs = 8,
    .bcs = elc_phys_bcs,

    .blocks = elc_blocks,
    .duplicate_across_blocks = false,

    .num_diag_moments = 4,
    .diag_moments = { GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_BIMAXWELLIAN },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,
    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },
  };

  // Ions.
  // Ion core source:
  struct gkyl_gyrokinetic_projection proj_srcCORE_i = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN  ,
    .gaussian_mean = {ctx.center_srcCORE[0], ctx.center_srcCORE[1], ctx.center_srcCORE[2]},
    .gaussian_std_dev = {ctx.sigma_srcCORE[0], ctx.sigma_srcCORE[1], ctx.sigma_srcCORE[2]},
    .total_num_particles = ctx.particle_srcCORE,
    .total_kin_energy = ctx.energy_srcCORE,
    .temp_max = 5.0*ctx.Ti0,
    .temp_min = 0.1*ctx.Ti0,
    .f_floor = ctx.floor_srcCORE,
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcCORE_i ={
    .adapt_to_species = "ion",
    .adapt_particle = ctx.adapt_particle_srcCORE,
    .adapt_energy = ctx.adapt_energy_srcCORE,
    .num_boundaries = 1,
    .dir = {0},
    .edge = {GKYL_LOWER_EDGE},
  };
  // Ion wall source:
  struct gkyl_gyrokinetic_projection proj_srcWALL_i = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN  ,
    .gaussian_mean = {ctx.center_srcWALL[0], ctx.center_srcWALL[1], ctx.center_srcWALL[2]},
    .gaussian_std_dev = {ctx.sigma_srcWALL[0], ctx.sigma_srcWALL[1], ctx.sigma_srcWALL[2]},
    .total_num_particles = ctx.particle_srcWALL,
    .total_kin_energy = ctx.energy_srcWALL,
    .temp_max = 5.0*ctx.Ti0,
    .temp_min = ctx.temp_recycle_srcWALL,
    .f_floor = ctx.floor_srcWALL,
  };
  struct gkyl_gyrokinetic_adapt_source adapt_srcWALL_i ={
    .adapt_to_species = "ion",
    .adapt_particle = ctx.adapt_particle_srcWALL,
    .adapt_energy = ctx.adapt_energy_srcWALL,
    .num_boundaries = 3,
    .dir = {0, 2, 2},
    .edge = {GKYL_UPPER_EDGE, GKYL_LOWER_EDGE, GKYL_UPPER_EDGE},
  };

  struct gkyl_gyrokinetic_projection ion_ic = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
    .ctx_density = &ctx,
    .density = eval_density,
    .ctx_upar = &ctx,
    .upar = eval_upar,
    .ctx_temp = &ctx,
    .temp = eval_temp_ion,
  };

  struct gkyl_gyrokinetic_multib_species_pb ion_blocks[3];

  ion_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {
    .block_id = 0,
    .polarization_density = ctx.n0,
    .projection = ion_ic,
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = ctx.num_sources,
      .num_adapt_sources = ctx.num_sources,
      .projection[0] = proj_srcCORE_i,
      .adapt[0] = adapt_srcCORE_i,
      .projection[1] = proj_srcWALL_i,
      .adapt[1] = adapt_srcWALL_i,
      .diagnostics = {
        .num_diag_moments = 1,
        .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
      }
    },
  };

  ion_blocks[1] = (struct gkyl_gyrokinetic_multib_species_pb) {
    .block_id = 1,
    .polarization_density = ctx.n0,
    .projection = ion_ic,
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = ctx.num_sources,
      .num_adapt_sources = ctx.num_sources,
      .projection[0] = proj_srcCORE_i,
      .adapt[0] = adapt_srcCORE_i,
      .projection[1] = proj_srcWALL_i,
      .adapt[1] = adapt_srcWALL_i,
      .diagnostics = {
        .num_diag_moments = 1,
        .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
      }
    },
  };

  ion_blocks[2] = (struct gkyl_gyrokinetic_multib_species_pb) {
    .block_id = 2,
    .polarization_density = ctx.n0,
    .projection = ion_ic,
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = ctx.num_sources,
      .num_adapt_sources = ctx.num_sources,
      .projection[0] = proj_srcCORE_i,
      .adapt[0] = adapt_srcCORE_i,
      .projection[1] = proj_srcWALL_i,
      .adapt[1] = adapt_srcWALL_i,
      .diagnostics = {
        .num_diag_moments = 1,
        .diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = {GKYL_F_MOMENT_HAMILTONIAN},
      }
    },
  };

  struct gkyl_gyrokinetic_bc ion_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 2, .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
  };

  struct gkyl_gyrokinetic_multib_species ion = {
    .name = "ion",
    .charge = ctx.qi, .mass = ctx.mi,
    .vdim = ctx.vdim,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0},
    .cells = { cells_v[0], cells_v[1] },

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm.
      .temp_ref = ctx.Ti0, // Temperature used to calculate coulomb logarithm.
      .bmag_ref = ctx.B0,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
    },

    .num_physical_bcs = 8,
    .bcs = ion_phys_bcs,

    .blocks = ion_blocks,
    .duplicate_across_blocks = false,

    .num_diag_moments = 4,
    .diag_moments = { GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_BIMAXWELLIAN },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,
    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },
  };


  // Field object
  struct gkyl_gyrokinetic_multib_field_pb field_blocks[1];
  field_blocks[0] = (struct gkyl_gyrokinetic_multib_field_pb) {
    .polarization_bmag = 2.18,
  };

  struct gkyl_gyrokinetic_bc field_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 1, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
  };

  struct gkyl_gyrokinetic_multib_field field = {
    .duplicate_across_blocks = true,
    .blocks = field_blocks, 
    .num_physical_bcs = 6, 
    .bcs = field_phys_bcs,
    .time_rate_diagnostics = true,
  };

  struct gkyl_gyrokinetic_multib *app_inp = gkyl_calloc(1,sizeof(struct gkyl_gyrokinetic_multib));
  app_inp->cdim = ctx.cdim;
  app_inp->poly_order = 1;
  app_inp->basis_type = app_args.basis_type;
  app_inp->cfl_frac = 1.0;
  app_inp->gk_block_geom = bgeom;
  app_inp->num_species = 2;
  app_inp->species[0] = elc;
  app_inp->species[1] = ion;
  app_inp->num_periodic_dir=1;
  app_inp->periodic_dirs[0]=1;
  app_inp->field = field;
  app_inp->comm = comm;
  app_inp->use_gpu = app_args.use_gpu;

  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp->name, sizeof(app_inp->name), "%s", app_args.app_name);
  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_type = GKYL_GK_MULTIB,
    .multib_app_inp = *app_inp,
    .time_stepping = {
      .t_end = ctx.t_end,
      .num_frames = ctx.num_frames,
      .write_phase_freq = ctx.write_phase_freq,
      .int_diag_calc_num = ctx.int_diag_calc_num,
      .dt_failure_tol = ctx.dt_failure_tol,
      .num_failures_max = ctx.num_failures_max,
      .is_restart = app_args.is_restart,
      .restart_frame = app_args.restart_frame,
      .num_steps = app_args.num_steps,
    },
  };

  gkyl_gyrokinetic_run_simulation(&run_inp);

  gkyl_gk_block_geom_release(bgeom);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif
  
  return 0;

}
