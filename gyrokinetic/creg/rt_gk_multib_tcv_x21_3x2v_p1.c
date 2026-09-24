#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_efit.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_mpi_comm.h>
#include <gkyl_null_comm.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>

struct gk_tcv_ctx {
  int cdim, vdim; // Dimensionality.

  double charge_elc; // Electron charge.
  double charge_ion; // Ion charge.
  double mass_elc; // Electron mass.
  double mass_ion; // Ion mass.
  double Te_max; // Maximum electron temperature
  double Ti_max; // Maximum ion temperature
  double n_max; // Maximum density.
  double Te0; // Reference electron temperature
  double Ti0; // Reference ion temperature
  double n0; // Reference density.
  double B0; // Magnetic field.
  double c_s; // Sound speed
  double nu_frac; // Factor multiplying collision frequencies.

  // Source parameters
  double psi_src; // Location.
  double lambda_src; // Width.
  double ndot_src; // Particle source rate.
  double Te_src; // Electron source temperature.
  double Ti_src; // Ion source temperature.

  // Domain parameters.            
  char geqdsk_file[128]; // File with equilibrium.
  double psi_axis; // Psi at the magnetic axis.
  double psi_sep; // Psi at the separatrix.
  double psi_max_core, psi_min_sol, psi_max_pf; // Psi extents.
  double Lx_core; // Box size in y.
  double Ly; // Box size in y.
  double y_min, y_max; // Limits in y.

  // Grid.
  int Npsi_sol; // Number of cells in psi in the SOL.
  int Npsi_pf; // Number of cells in psi in the private flux.
  int Npsi_core; // Number of cells in psi in the core.
  int Ntheta_divertor_lfs; // Number of cells in theta in the divertor, LFS.
  int Ntheta_divertor_hfs; // Number of cells in theta in the divertor, HFS.
  int Ntheta_sol; // Number of cells in theta in the (upper) SOL (and core).
  int Nvpar, Nmu; // Number of cells in vpar,mu.
  int Ny; // Number of cells in y.
  int cells_v[2]; // Number of cells in all directions.
  int num_blocks; // Number of blocks.

  double parallel_compression_fac; // Compress cells near X-pt (0 to 1, 0 is full compression).
  double radial_compression_fac; // Compress cells separatrix (0 to 1, 0 is full compression).

  // Physical velocity space limits
  double vpar_max_elc; // Parallel velocity extents for electrons.
  double mu_max_elc; // Maximum magnetic moment for electrons.
  double vpar_max_ion; // Parallel velocity extents for ions.
  double mu_max_ion; // Maximum magnetic moment for ions.

  // Computational velocity space limits
  double vpar_min_elc_c, vpar_max_elc_c;
  double mu_min_elc_c, mu_max_elc_c;
  double vpar_min_ion_c, vpar_max_ion_c;
  double mu_min_ion_c, mu_max_ion_c;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

void divertor_plate_func_out(double s, double* RZ)
{
  // Straight bottom plate.
  double R_lo = 0.6662613750;
  double R_up = 0.9668828249;
  double Z = -0.75;

  RZ[0] = R_lo + (R_up-R_lo)*s;
  RZ[1] = Z;
}

void divertor_plate_func_in(double s, double* RZ)
{
  // Straight left plate.
  double R = 0.6240000129;
  double Z_lo = -0.7014500499;
  double Z_up =  -0.4014500499;

  RZ[0] = R;
  RZ[1] = Z_lo + (Z_up-Z_lo)*s;
}

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

struct gkyl_gk_block_geom*
create_asdex_lsn_gk_block_geom(void *ctx)
{
  struct gk_tcv_ctx *params = ctx;

  struct gkyl_gk_block_geom *bgeom = gkyl_gk_block_geom_new(params->cdim, params->num_blocks);

  /* Block layout and coordinates.

    z  
    ^  
    |

    |                     +------------------+------------------+
    |                     |                  |                  | 
    |                     |b4                |b3                |
    |                     |inner PF          |lower inner sol   |
    |                     |                  |                  | 
    |                     |%%%%%%%%%%%%%%%%%%|                  |
    |  +------------------+------------------+------------------+
    |  |                 $|                  |$                 |
    |  |                 $|                  |$                 |
    |  | b5              $|                  |$ b2              |
    |  + core            $+                  +$ upper sol       |
    |  |                 $|                  |$                 |
    |  |                 $|                  |$                 |
    |  +------------------+------------------+------------------+
    |                     |%%%%%%%%%%%%%%%%%%|                  |
    |                     |                  |                  | 
    |                     | b0               |b1                |
    |                     | outer PF         |lower outer sol   |
    |                     |                  |                  |
    0                     +------------------+------------------+

       0 -------------------------------------------------------- -> x

    Edges that touch coincide are physically connected unless
    otherwise indicated by a special symbol. Edges with a special
    symbol such as o,x,%, or % are instead connected to the other
    edge with the same symbol. Edges that do not coincide with
    another edge are a physical boundary.
  */  

  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };
  // Copy eqdsk file into efit_inp.
  memcpy(efit_inp.filepath, params->geqdsk_file, sizeof(params->geqdsk_file));

  // Theta limits are actually set internally by the code.
  double theta_min = -1.0, theta_max = 1.0;

  double psi_sep  = params->psi_sep; // Psi at the separatrix.
  double psi_axis = params->psi_axis; // Psi at the magnetic axis.
  double psi_max_core = params->psi_max_core; // Minimum psi the core.
  double psi_min_sol  = params->psi_min_sol ; // Maximum psi the SOL.
  double psi_max_pf   = params->psi_max_pf  ; // Minimum psi the private flux.

  double y_min = params->y_min; // Lower boundary in y.
  double y_max = params->y_max; // Upper boundary in y.

  // Number of cells.
  int Npsi_sol            = params->Npsi_sol       ;
  int Npsi_pf             = params->Npsi_pf        ;
  int Npsi_core           = params->Npsi_core      ;
  int Ntheta_divertor_lfs = params->Ntheta_divertor_lfs;
  int Ntheta_divertor_hfs = params->Ntheta_divertor_hfs;
  int Ntheta_sol          = params->Ntheta_sol     ;
  int Ny                  = params->Ny             ;

  // Block 0: outer private flux (PF) region.
  gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_sep, y_min, theta_min },
      .upper = { psi_max_pf, y_max, theta_max },
      .cells = { Npsi_pf, Ny, Ntheta_divertor_lfs },
      .cuts = { 1, 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_PF_LO_R,
          .rmin       = 0.618,
          .rmax       = 1.14,
          .zmin       = -0.75,
          .zmax       =  -0.38,
          .rleft      = 0.618,
          .rright     = 1.14,
          .zmin_left  = -0.45,
          .zmin_right = -0.75,
          .plate_spec = true,
          .plate_func_lower = divertor_plate_func_out,
          .plate_func_upper = divertor_plate_func_in,
        },
        .position_map_info = {
          .id = GKYL_PMAP_XPT_COMPRESSION,
          .compress_divertor = true,
          .radial_compression_factor = params->radial_compression_fac,
          .compression_factor = params->parallel_compression_fac
        },
      },

      .connections[0] = { // x-direction.
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL },  // Physical boundary.
        { .bid = 1, .dir = 0, .edge = GKYL_LOWER_POSITIVE },
      },
      .connections[1] = { // y-direction.
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL },
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
      },
      .connections[2] = { // z-direction.
        { .bid = 0, .dir = 2, .edge = GKYL_PHYSICAL}, // Physical boundary.
        { .bid = 4, .dir = 2, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // Block 1: lower outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 1, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_min_sol, y_min, theta_min },
      .upper = { psi_sep, y_max, theta_max },
      .cells = { Npsi_sol, Ny, Ntheta_divertor_lfs },
      .cuts = { 1, 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO,
          .rmin   = 0.618,
          .rmax   = 1.14,
          .zmin   = -0.75,
          .zmax   =  0.35,
          .rclose = 0.825,
          .rleft  = 0.62,
          .rright = 0.9,
          .zmin_left  = -0.45,
          .zmin_right = -0.75,
          .plate_spec = true,
          .plate_func_lower = divertor_plate_func_out,
          .plate_func_upper = divertor_plate_func_in,
        },
        .position_map_info = {
          .id = GKYL_PMAP_XPT_COMPRESSION,
          .compress_divertor = true,
          .radial_compression_factor = params->radial_compression_fac,
          .compression_factor = params->parallel_compression_fac
        },
      },
      
      .connections[0] = { // x-direction.
        { .bid = 0, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 1, .dir = 0, .edge = GKYL_PHYSICAL }, // Physical boundary.
      },
      .connections[1] = { // y-direction.
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL },
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
      },
      .connections[2] = { // z-direction.
        { .bid = 1, .dir = 2, .edge = GKYL_PHYSICAL}, // Physical boundary.
        { .bid = 2, .dir = 2, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // Block 2: mid SOL.
  gkyl_gk_block_geom_set_block(bgeom, 2, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_min_sol, y_min, theta_min },
      .upper = { psi_sep, y_max, theta_max },
      .cells = { Npsi_sol, Ny, Ntheta_sol },
      .cuts = { 1, 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID,
          .rmin   = 0.618,
          .rmax   = 1.14,
          .zmin   = -0.75,
          .zmax   =  0.35,
          .rclose = 0.825,
          .rleft  = 0.62,
          .rright = 0.9,
          .zmin_left  = -0.45,
          .zmin_right = -0.75,
          .plate_spec = true,
          .plate_func_lower = divertor_plate_func_out,
          .plate_func_upper = divertor_plate_func_in,
        },
        .position_map_info = {
          .id = GKYL_PMAP_XPT_COMPRESSION,
          .compress_divertor = true,
          .radial_compression_factor = params->radial_compression_fac,
          .compression_factor = params->parallel_compression_fac
        },
      },
      
      .connections[0] = { // x-direction.
        { .bid = 5, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 2, .dir = 0, .edge = GKYL_PHYSICAL }, // Physical boundary.
      },
      .connections[1] = { // y-direction.
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL },
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
      },
      .connections[2] = { // z-direction.
        { .bid = 1, .dir = 2, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 3, .dir = 2, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // Block 3: lower inner SOL.
  gkyl_gk_block_geom_set_block(bgeom, 3, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_min_sol, y_min, theta_min },
      .upper = { psi_sep, y_max, theta_max },
      .cells = { Npsi_sol, Ny, Ntheta_divertor_hfs },
      .cuts = { 1, 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP,
          .rmin   = 0.618,
          .rmax   = 1.14,
          .zmin   = -0.75,
          .zmax   =  0.35,
          .rclose = 0.825,
          .rleft  = 0.62,
          .rright = 0.9,
          .zmin_left  = -0.45,
          .zmin_right = -0.75,
          .plate_spec = true,
          .plate_func_lower = divertor_plate_func_out,
          .plate_func_upper = divertor_plate_func_in,
        },
        .position_map_info = {
          .id = GKYL_PMAP_XPT_COMPRESSION,
          .compress_divertor = true,
          .radial_compression_factor = params->radial_compression_fac,
          .compression_factor = params->parallel_compression_fac
        },
      },
      
      .connections[0] = { // x-direction.
        { .bid = 4, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 3, .dir = 0, .edge = GKYL_PHYSICAL }, // Physical boundary.
      },
      .connections[1] = { // y-direction.
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL },
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
      },
      .connections[2] = { // z-direction.
        { .bid = 2, .dir = 2, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 3, .dir = 2, .edge = GKYL_PHYSICAL}, // Physical boundary.
      }
    }
  );

  // Block 4: inner private flux (PF) region.
  gkyl_gk_block_geom_set_block(bgeom, 4, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_sep, y_min, theta_min },
      .upper = { psi_max_pf, y_max, theta_max },
      .cells = { Npsi_pf, Ny, Ntheta_divertor_hfs },
      .cuts = { 1, 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_PF_LO_L,
          .rmin       = 0.618,
          .rmax       = 1.14,
          .zmin       = -0.75,
          .zmax       =  -0.38,
          .rleft      = 0.618,
          .rright     = 1.14,
          .zmin_left  = -0.45,
          .zmin_right = -0.75,
          .plate_spec = true,
          .plate_func_lower = divertor_plate_func_out,
          .plate_func_upper = divertor_plate_func_in,
        },
        .position_map_info = {
          .id = GKYL_PMAP_XPT_COMPRESSION,
          .compress_divertor = true,
          .radial_compression_factor = params->radial_compression_fac,
          .compression_factor = params->parallel_compression_fac
        },
      },

      .connections[0] = { // x-direction.
        { .bid = 4, .dir = 0, .edge = GKYL_PHYSICAL },  // Physical boundary.
        { .bid = 3, .dir = 0, .edge = GKYL_LOWER_POSITIVE },
      },
      .connections[1] = { // y-direction.
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL },
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
      },
      .connections[2] = { // z-direction.
        { .bid = 0, .dir = 2, .edge = GKYL_UPPER_POSITIVE}, // Physical boundary.
        { .bid = 4, .dir = 2, .edge = GKYL_PHYSICAL},
      }
    }
  );

  // Block 5: core region.
  gkyl_gk_block_geom_set_block(bgeom, 5, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_sep, y_min, theta_min },
      .upper = { psi_max_core, y_max, theta_max },
      .cells = { Npsi_core, Ny, Ntheta_sol },
      .cuts = { 1, 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_CORE,
          .rmin   = 0.618,
          .rmax   = 1.14,
          .zmin   = -0.75,
          .zmax   =  0.35,
          .rclose = 1.14,
          .rleft  = 0.618,
          .rright = 1.14,
        },
        .position_map_info = {
          .id = GKYL_PMAP_XPT_COMPRESSION,
          .compress_divertor = true,
          .radial_compression_factor = params->radial_compression_fac,
          .compression_factor = params->parallel_compression_fac
        },
      },

      .connections[0] = { // x-direction.
        { .bid = 5, .dir = 0, .edge = GKYL_PHYSICAL },  // Physical boundary.
        { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE },
      },
      .connections[1] = { // y-direction.
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL },
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
      },
      .connections[2] = { // z-direction.
        { .bid = 5, .dir = 2, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 5, .dir = 2, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  return bgeom;
}

// Velocity space mappings.
void mapc2p_vel_elc(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_tcv_ctx *app = ctx;
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
  struct gk_tcv_ctx *app = ctx;
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

double
init_profile(double psi, double f_min, double f_max, void *ctx)
{
  // Profile in D. Michels, et al. Phys. Plasmas 29, 032307 (2022) eqn 17:
  struct gk_tcv_ctx *params = ctx;
  double psi_min = params->psi_max_core;
  double psi_max = params->psi_min_sol;
  double psi_axis = params->psi_axis;
  double psi_sep = params->psi_sep;

  double rho_min = rho_psi(psi_min, psi_axis, psi_sep);
  double rho_max = rho_psi(psi_max, psi_axis, psi_sep);
  double rho = rho_psi(psi, psi_axis, psi_sep);

  double c1 = (f_max-f_min)/2.0;
  double c2 = M_PI/(rho_max-rho_min);
  double c3 = M_PI/2 - c2*rho_min;
  double c4 = (f_max+f_min)/2.0;

  double f = -1.0;
  if (rho <= rho_min)
    f = f_max;
  else if (rho_min < rho && rho <= rho_max)
    f = c1*sin(c2*rho + c3) + c4;
  else
    f = f_min;

  return f;
}

void
init_dens(double t, const double *xn, double* restrict fout, void *ctx)
{
  double psi = xn[0], alpha = xn[1], theta = xn[2];

  struct gk_tcv_ctx *params = ctx;
  double psi_axis = params->psi_axis;
  double psi_sep = params->psi_sep;
  double n_max = params->n_max;

  double rho = rho_psi(psi, psi_axis, psi_sep);

  fout[0] = n_max * 1.165 * 0.5 * (1.0 + (0.815/1.165) * tanh((0.905 - rho)/0.048));
}

void
init_temp_elc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double psi = xn[0], alpha = xn[1], theta = xn[2];

  struct gk_tcv_ctx *params = ctx;
  double psi_axis = params->psi_axis;
  double psi_sep = params->psi_sep;
  double Te_max = params->Te_max;

  double rho = rho_psi(psi, psi_axis, psi_sep);

  fout[0] = Te_max * 1.1694 * 0.5 * ( 1.0 + (95.05/116.94) * tanh((0.87 - rho)/0.066));
}

void
init_temp_ion(double t, const double *xn, double* restrict fout, void *ctx)
{
  double psi = xn[0], alpha = xn[1], theta = xn[2];

  struct gk_tcv_ctx *params = ctx;
  double psi_axis = params->psi_axis;
  double psi_sep = params->psi_sep;
  double Ti_max = params->Ti_max;

  double rho = rho_psi(psi, psi_axis, psi_sep);

  fout[0] = Ti_max * 1.1694 * 0.5 * ( 1.0 + (95.05/116.94) * tanh((0.87 - rho)/0.066));
}

void
init_upar(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 0.0;
}

void
init_source_dens(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double psi = xn[0], alpha = xn[1], theta = xn[2];

  struct gk_tcv_ctx *params = ctx;
  double lambda_src = params->lambda_src;
  double psi_src = params->psi_src;
  double ndot_src = params->ndot_src;

  double source_floor = 1e-10;
  if (psi < psi_src + 3*lambda_src)
    source_floor = 1e-2;

  double src_prof = exp(-pow(psi-psi_src,2)/(2*pow(lambda_src,2)));
  fout[0] = ndot_src * fmax(src_prof, source_floor);
}

void
init_source_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
init_source_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double psi = xn[0], alpha = xn[1], theta = xn[2];

  struct gk_tcv_ctx *params = ctx;
  double lambda_src = params->lambda_src;
  double psi_src = params->psi_src;
  double Te_src = params->Te_src;
  double eV = GKYL_ELEMENTARY_CHARGE;

  if (psi < psi_src + 3*lambda_src)
    fout[0] = Te_src;
  else
    fout[0] = 2.0*eV;
}

void
init_source_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_tcv_ctx *params = ctx;
  double psi = xn[0], alpha = xn[1], theta = xn[2];

  double lambda_src = params->lambda_src;
  double psi_src = params->psi_src;
  double Ti_src = params->Ti_src;
  double eV = GKYL_ELEMENTARY_CHARGE;

  if (psi < psi_src + 3*lambda_src)
    fout[0] = Ti_src;
  else
    fout[0] = 2.0*eV;
}

struct gk_tcv_ctx
create_ctx(void)
{
  int cdim = 3, vdim = 2; // Dimensionality.

  double eps0 = GKYL_EPSILON0; // Permittivity of free space.
  double eV = GKYL_ELEMENTARY_CHARGE; // Elementary charge.
  double mi = 2.014*GKYL_PROTON_MASS; // Ion mass.
  double me = GKYL_ELECTRON_MASS; // Electron mass.
  double qi = eV; // Ion charge.
  double qe = -eV; // Electron charge.

  double Te_max = 200.0*eV; // Maximum electron temperature.
  double Ti_max = 200.0*eV; // Maximum ion temperature.
  double B_max  = 1.561831e+00; // Maximum B field amplitude.
  double n_max  = 2.0e19; // Maximum particle density.

  double Te_min = 20.0*eV; // Minimum electron temperature.
  double Ti_min = 20.0*eV; // Minimum ion temperature.
  double B_min  = 9.689716e-01; // Minimum B field amplitude.
  double n_min  = 0.2e19; // Minimum particle density.

  double Te0 = 0.5*(Te_min+Te_max); // Reference electron temperature.
  double Ti0 = 0.5*(Ti_min+Ti_max); // Reference ion temperature.
  double B0  = 0.5*(B_min +B_max ); // Reference B field amplitude.
  double n0  = 0.5*(n_min +n_max ); // Reference particle density.

  // Derived parameters.
  double vt_ion = sqrt(Ti0/mi);
  double vt_elc = sqrt(Te0/me);
  double c_s = sqrt(Te0/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;

  // Collision parameters.
  double nu_frac = 1.0;  
  double logLambda_elc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te0/eV);
  double nu_elc = nu_frac*logLambda_elc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*pow(eps0,2)*sqrt(me)*pow(Te0,3.0/2.0));

  double logLambda_ion = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti0/eV);
  double nu_ion = nu_frac*logLambda_ion*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*pow(eps0,2)*sqrt(mi)*pow(Ti0,3.0/2.0));

  // Location of the numerical equilibrium.
  char geqdsk_file[128] = "gyrokinetic/data/eqdsk/tcv_65402_t1.eqdsk";

  // Position space parameters.
  double num_blocks = 6;

  // Get the separatrix psi.
  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };
  // Copy eqdsk file into efit_inp.
  memcpy(efit_inp.filepath, geqdsk_file, sizeof(geqdsk_file));
  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);
  double psi_sep = efit->psisep;
  double psi_axis = efit->simag;
  double R_axis = efit->rmaxis;
  double Rxpt = efit->Rxpt[0], Zxpt = efit->Zxpt[0];
  gkyl_efit_release(efit);
  // Here rho = sqrt((psi-psi_axis) / (psi_sep - psi_axis)).
  double rho_min_core = 0.9;
  double rho_max_sol = 1.04;
  double psi_max_core = psi_rho(rho_min_core, psi_axis, psi_sep);
  double psi_min_sol = psi_rho(rho_max_sol, psi_axis, psi_sep);
  double psi_max_pf = psi_sep + 0.5*(psi_sep-psi_min_sol);

  // Number of cells.
  int Npsi_sol = 4;
  int Npsi_pf = 4;
  int Npsi_core = 8;
  int Ntheta_divertor_lfs = 8;
  int Ntheta_divertor_hfs = 4;
  int Ntheta_sol = 24;
  int Ny = 2;
  int Nvpar = 4; // Number of cells in vpar.
  int Nmu = 4; // Number of cells in mu.

  double parallel_compression_fac = 0.5; // Compress cells near X-pt (0 to 1, 0 is full compression).
  double radial_compression_fac = 0.5; // Compress cells separatrix (0 to 1, 0 is full compression).

  // Adjust psi so that psi_LCFS is at a cell boundary.
  double delta_psi = (psi_sep-psi_min_sol)/Npsi_sol; // Cell length in psi.
  psi_max_core = psi_sep + Npsi_core*delta_psi; // Adjust inner radial core boundary.
  rho_min_core = rho_psi(psi_max_core, psi_axis, psi_sep); // Minimum psi_N.

  double Lx_core = psi_sep - psi_max_core;

  double r0 = 1.093 - R_axis;
  double q0 = 3.2; // q95, see Oliveira 2022.
  double Ly = 100*rho_s*q0/r0/3.0;
  double y_min = -Ly/2.;
  double y_max =  Ly/2.;

  // Source parameters.
  double power_in = 150e3; // Input power [W].
  double ndot_src = 1.85e21; // Input particle source rate [particles/s].
  double Te_src = 0.5*power_in/((3.0/2.0)*ndot_src);
  double Ti_src = 0.5*power_in/((3.0/2.0)*ndot_src);
  double psi_src = psi_max_core;
  double lambda_src = psi_rho(0.915, psi_axis, psi_sep) - psi_max_core;

  // Physical velocity space limits
  double vpar_max_elc = 6.0*vt_elc;
  double mu_max_elc = me*pow(4.0*vt_elc,2)/(2.0*B0);

  double vpar_max_ion = 6.0*vt_ion;
  double mu_max_ion = mi*pow(4.0*vt_ion,2)/(2.0*B0);

  // Computational velocity space limits.
  double vpar_min_ion_c = -1.0/sqrt(2.0);
  double vpar_max_ion_c = 1.0/sqrt(2.0);
  double mu_min_ion_c = 0.;
  double mu_max_ion_c = 1.;
  // Computational velocity space limits.
  double vpar_min_elc_c = -1.0/sqrt(2.0);
  double vpar_max_elc_c = 1.0/sqrt(2.0);
  double mu_min_elc_c = 0.;
  double mu_max_elc_c = 1.;

  printf("  X-point @ (R,Z)     = (%.9e,%9e)\n",Rxpt,Zxpt);
  printf("  psi_axis            = %.13e\n",psi_axis);
  printf("  psi_sep             = %.13e\n",psi_sep);
  printf("  rho_min_core        = %.13e\n",rho_min_core);
  printf("  rho_max_sol         = %.13e\n",rho_max_sol);
  printf("  psi_max_core        = %.13e\n",psi_max_core);
  printf("  psi_min_sol         = %.13e\n",psi_min_sol);
  printf("  psi_max_pf          = %.13e\n",psi_max_pf);
  printf("  Npsi_sol            = %d\n",Npsi_sol       );
  printf("  Npsi_pf             = %d\n",Npsi_pf        );
  printf("  Npsi_core           = %d\n",Npsi_core      );
  printf("  Ntheta_divertor_lfs = %d\n",Ntheta_divertor_lfs);
  printf("  Ntheta_divertor_hfs = %d\n",Ntheta_divertor_hfs);
  printf("  Ntheta_sol          = %d\n",Ntheta_sol     );

  double t_end = 1.0e-7;
  double num_frames = 1;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_tcv_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .charge_elc = qe, 
    .charge_ion = qi, 
    .mass_elc = me, 
    .mass_ion = mi,
    .Te_max = Te_max, 
    .Ti_max = Ti_max, 
    .n_max = n_max, 
    .Te0 = Te0, 
    .Ti0 = Ti0, 
    .n0 = n0, 
    .B0 = B0, 
    .c_s = c_s, 
    .nu_frac = nu_frac, 
    .num_blocks = num_blocks,
    .psi_axis = psi_axis,
    .psi_sep = psi_sep,
    .psi_max_core = psi_max_core,
    .psi_min_sol = psi_min_sol,
    .psi_max_pf = psi_max_pf,
    .Lx_core = Lx_core, 
    .Ly = Ly,
    .y_min = y_min,
    .y_max = y_max,
    .lambda_src = lambda_src,
    .psi_src = psi_src,
    .ndot_src = ndot_src,
    .Te_src = Te_src,
    .Ti_src = Ti_src,
    // Physical velocity space limits
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
    // Computational velocity space limits
    .vpar_min_elc_c = vpar_min_elc_c,
    .vpar_max_elc_c = vpar_max_elc_c,
    .mu_min_elc_c = mu_min_elc_c,
    .mu_max_elc_c = mu_max_elc_c,
    .vpar_min_ion_c = vpar_min_ion_c,
    .vpar_max_ion_c = vpar_max_ion_c,
    .mu_min_ion_c = mu_min_ion_c,
    .mu_max_ion_c = mu_max_ion_c,
    .Npsi_sol        = Npsi_sol       ,
    .Npsi_pf         = Npsi_pf        ,
    .Npsi_core       = Npsi_core      ,
    .Ntheta_divertor_lfs = Ntheta_divertor_lfs,
    .Ntheta_divertor_hfs = Ntheta_divertor_hfs,
    .Ntheta_sol      = Ntheta_sol     ,
    .Ny              = Ny             ,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells_v = {Nvpar, Nmu},
    .parallel_compression_fac = parallel_compression_fac,
    .radial_compression_fac   = radial_compression_fac,
    .t_end = t_end, 
    .num_frames = num_frames, 
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  // Copy eqdsk file into ctx.
  memcpy(ctx.geqdsk_file, geqdsk_file, sizeof(geqdsk_file));
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

  struct gk_tcv_ctx ctx = create_ctx(); // Context for init functions.
                    
  // Construct block geometry
  struct gkyl_gk_block_geom *bgeom = create_asdex_lsn_gk_block_geom(&ctx);

  int cells_v[ctx.vdim];
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells_v[d]);

  struct gkyl_gyrokinetic_projection elc_ic = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
    .ctx_density = &ctx,
    .density = init_dens,
    .ctx_upar = &ctx,
    .upar = init_upar,
    .ctx_temp = &ctx,
    .temp = init_temp_elc,
  };

  struct gkyl_gyrokinetic_multib_species_pb elc_sol = {
    .polarization_density = ctx.n0,
    .projection = elc_ic,
  };

  // Electrons.
  struct gkyl_gyrokinetic_multib_species_pb elc_blocks[ctx.num_blocks];

  elc_sol.block_id = 0;
  elc_blocks[0] = elc_sol;

  elc_sol.block_id = 1;
  elc_blocks[1] = elc_sol;

  elc_sol.block_id = 2;
  elc_blocks[2] = elc_sol;

  elc_sol.block_id = 3;
  elc_blocks[3] = elc_sol;

  elc_sol.block_id = 4;
  elc_blocks[4] = elc_sol;

  elc_blocks[5] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 5,
    .polarization_density = ctx.n0,

    .projection = elc_ic,

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = init_source_dens,
        .ctx_upar = &ctx,
        .upar = init_source_upar,
        .ctx_temp = &ctx,
        .temp = init_source_temp_elc,
      }, 
      .diagnostics = {
        .num_diag_moments = 1,
        .diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP, },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//        .time_integrated = true,
      },
    },

  };

  struct gkyl_gyrokinetic_bc elc_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 1, .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 3, .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
  };

  struct gkyl_gyrokinetic_multib_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .vdim = ctx.vdim,
    .lower = { ctx.vpar_min_elc_c, ctx.mu_min_elc_c},
    .upper = { ctx.vpar_max_elc_c, ctx.mu_max_elc_c},
    .cells = { cells_v[0], cells_v[1] },

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES_NO_BY,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm.
      .temp_ref = ctx.Te0, // Temperature used to calculate coulomb logarithm.
      .bmag_ref = ctx.B0,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },

    .num_physical_bcs = 10,
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
  struct gkyl_gyrokinetic_projection ion_ic = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
    .ctx_density = &ctx,
    .density = init_dens,
    .ctx_upar = &ctx,
    .upar = init_upar,
    .ctx_temp = &ctx,
    .temp = init_temp_ion,
  };

  struct gkyl_gyrokinetic_multib_species_pb ion_sol = {
    .polarization_density = ctx.n0,
    .projection = ion_ic,
  };

  struct gkyl_gyrokinetic_multib_species_pb ion_blocks[ctx.num_blocks];
  ion_sol.block_id = 0;
  ion_blocks[0] = ion_sol;

  ion_sol.block_id = 1;
  ion_blocks[1] = ion_sol;

  ion_sol.block_id = 2;
  ion_blocks[2] = ion_sol;

  ion_sol.block_id = 3;
  ion_blocks[3] = ion_sol;

  ion_sol.block_id = 4;
  ion_blocks[4] = ion_sol;

  ion_blocks[5] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 5,
    .polarization_density = ctx.n0,

    .projection = ion_ic,

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = init_source_dens,
        .ctx_upar = &ctx,
        .upar = init_source_upar,
        .ctx_temp = &ctx,
        .temp = init_source_temp_ion,
      }, 
      .diagnostics = {
        .num_diag_moments = 1,
        .diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP, },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//        .time_integrated = true,
      },
    },

  };

  struct gkyl_gyrokinetic_bc ion_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 1, .dir = 2, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 3, .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 2, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
  };

  struct gkyl_gyrokinetic_multib_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .vdim = ctx.vdim,
    .lower = { ctx.vpar_min_ion_c, ctx.mu_min_ion_c},
    .upper = { ctx.vpar_max_ion_c, ctx.mu_max_ion_c},
    .cells = { cells_v[0], cells_v[1] },

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES_NO_BY,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm.
      .temp_ref = ctx.Ti0, // Temperature used to calculate coulomb logarithm.
      .bmag_ref = ctx.B0,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
    },

    .num_physical_bcs = 10,
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

  // Field object.
  struct gkyl_gyrokinetic_multib_field_pb field_blocks[1];
  field_blocks[0] = (struct gkyl_gyrokinetic_multib_field_pb) {
    // No block specific field info for this simulation
    .time_rate_diagnostics = true,
  };

  struct gkyl_gyrokinetic_bc field_phys_bcs[] = {
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 4, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 0, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, .value = {0.0} },
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, .value = {0.0} },
    { .bidx = 2, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, .value = {0.0} },
    { .bidx = 3, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, .value = {0.0} },
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, .value = {0.0} },
    { .bidx = 5, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, .value = {0.0} },
  };

  struct gkyl_gyrokinetic_multib_field field = {
    .blocks = field_blocks, 
    .duplicate_across_blocks = true,

    .num_physical_bcs = ctx.num_blocks,
    .bcs = field_phys_bcs,
  };

  struct gkyl_gyrokinetic_multib *app_inp = gkyl_calloc(1,sizeof(struct gkyl_gyrokinetic_multib));
  app_inp->cdim = ctx.cdim;
  app_inp->poly_order = 1;
  app_inp->basis_type = app_args.basis_type;
  app_inp->use_gpu = app_args.use_gpu;
  app_inp->cfl_frac = 1.0;
  app_inp->gk_block_geom = bgeom;
  app_inp->num_species = 2;
  app_inp->species[0] = elc;
  app_inp->species[1] = ion;
  app_inp->field = field;
  app_inp->comm = comm;

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
  gkyl_free(app_inp);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif
  
  return 0;
}
