#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_calc_bmag.h>
#include <gkyl_comm_io.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_position_map.h>
#include <gkyl_position_map_priv.h>

#include <float.h>
#include <assert.h>

// Remove with the print statements at the bottom
#include <gkyl_util.h>

void
gkyl_position_map_identity(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = xn[0];
}

void
gkyl_position_map_identity_slope(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = 1.0;
}


struct gkyl_position_map*
gkyl_position_map_null_new()
{
  struct gkyl_position_map *gpm = gkyl_calloc(1, sizeof(*gpm));
  gpm->id = GKYL_PMAP_USER_INPUT;

  gpm->use_map_derivs = true;
  gpm->to_optimize = false;
  gpm->mc2nu = gkyl_array_new(GKYL_DOUBLE, 1, 1);
  gpm->constB_ctx = gkyl_calloc(1, sizeof(struct gkyl_position_map_const_B_ctx));
  gpm->xpt_ctx = gkyl_calloc(1, sizeof(struct gkyl_position_map_xpt_ctx));
  gpm->bmag_ctx = gkyl_calloc(1, sizeof(struct gkyl_bmag_ctx));
  gpm->bmag_ctx->bmag = gkyl_array_new(GKYL_DOUBLE, 1, 1);
  gpm->ref_count = gkyl_ref_count_init(gkyl_position_map_free);
  
  for (int i = 0; i < 3; i++){
    gpm->maps[i] = gkyl_position_map_identity;
    gpm->map_derivs[i] = gkyl_position_map_identity_slope;
    gpm->ctxs[i] = 0;
    gpm->constB_ctx->maps_backup[i] = gkyl_position_map_identity;
    gpm->constB_ctx->ctxs_backup[i] = 0;
    gpm->xpt_ctx->maps_backup[i] = gkyl_position_map_identity;
    gpm->xpt_ctx->map_derivs_backup[i] = gkyl_position_map_identity_slope;
    gpm->xpt_ctx->ctxs_backup[i] = 0;
  }
  return gpm;
}

struct gkyl_position_map*
gkyl_position_map_inew(struct gkyl_position_map_inew_inp inp)
{
  return gkyl_position_map_new(inp.pmap_info, inp.grid, inp.local, inp.local_ext,
    inp.global, inp.global_ext, inp.basis);
}

struct gkyl_position_map*
gkyl_position_map_new(struct gkyl_position_map_inp pmap_info, struct gkyl_rect_grid grid,
  struct gkyl_range local, struct gkyl_range local_ext, struct gkyl_range global, struct gkyl_range global_ext,
  struct gkyl_basis basis)
{
  struct gkyl_position_map *gpm = gkyl_calloc(1, sizeof(*gpm));
  gpm->id = pmap_info.id;

  gpm->bmag_ctx = gkyl_calloc(1, sizeof(struct gkyl_bmag_ctx));
  gpm->bmag_ctx->bmag = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, global_ext.volume);
  gpm->to_optimize = false;
  gpm->use_map_derivs = (pmap_info.id == GKYL_PMAP_XPT_COMPRESSION || pmap_info.id == GKYL_PMAP_USER_INPUT_W_DERIVATIVE) ? true : false;


  gpm->constB_ctx = gkyl_calloc(1, sizeof(struct gkyl_position_map_const_B_ctx));
  gpm->xpt_ctx = gkyl_calloc(1, sizeof(struct gkyl_position_map_xpt_ctx));

  for (int i = 0; i < 3; i++){
    gpm->maps[i] = gkyl_position_map_identity;
    gpm->map_derivs[i] = gkyl_position_map_identity_slope;
    gpm->ctxs[i] = 0;
    gpm->constB_ctx->maps_backup[i] = gkyl_position_map_identity;
    gpm->constB_ctx->ctxs_backup[i] = 0;
    gpm->xpt_ctx->maps_backup[i] = gkyl_position_map_identity;
    gpm->xpt_ctx->map_derivs_backup[i] = gkyl_position_map_identity_slope;
    gpm->xpt_ctx->ctxs_backup[i] = 0;
  }

  switch (pmap_info.id)
  {
    case GKYL_PMAP_USER_INPUT:
    case GKYL_PMAP_USER_INPUT_W_DERIVATIVE:
      for (int i = 0; i < 3; i++) {
        if (pmap_info.maps[i]) {
          gpm->maps[i] = pmap_info.maps[i];
          gpm->map_derivs[i] = pmap_info.map_derivs[i];
          gpm->ctxs[i] = pmap_info.ctxs[i];
        }
      }
      break;

    case GKYL_PMAP_CONSTANT_DB_POLYNOMIAL:
    case GKYL_PMAP_CONSTANT_DB_NUMERIC:

      for (int i = 0; i < 2; i++){
        if (pmap_info.maps[i] != 0)
        { gpm->constB_ctx->maps_backup[i] = pmap_info.maps[i];
          gpm->constB_ctx->ctxs_backup[i] = pmap_info.ctxs[i];
        }
      }
      gpm->constB_ctx->map_strength = pmap_info.map_strength;

      if (pmap_info.maximum_slope_at_min_B == 0.)
      { gpm->constB_ctx->enable_maximum_slope_limits_at_min_B = false;  }
      else
      { gpm->constB_ctx->enable_maximum_slope_limits_at_min_B = true; }
      gpm->constB_ctx->maximum_slope_at_min_B = pmap_info.maximum_slope_at_min_B;

      if (pmap_info.maximum_slope_at_max_B == 0.)
      { gpm->constB_ctx->enable_maximum_slope_limits_at_max_B = false;  }
      else
      { gpm->constB_ctx->enable_maximum_slope_limits_at_max_B = true; }
      gpm->constB_ctx->maximum_slope_at_max_B = pmap_info.maximum_slope_at_max_B;
      gpm->constB_ctx->gaussian_std = pmap_info.gaussian_std;
      gpm->constB_ctx->gaussian_max_integration_width = pmap_info.gaussian_max_integration_width;

      break;

    case GKYL_PMAP_XPT_COMPRESSION:

      for (int i = 0; i < 3; i++){
        if (pmap_info.maps[i] != 0)
        { gpm->xpt_ctx->maps_backup[i] = pmap_info.maps[i];
          gpm->xpt_ctx->map_derivs_backup[i] = pmap_info.map_derivs[i];
          gpm->xpt_ctx->ctxs_backup[i] = pmap_info.ctxs[i];
        }
      }
      gpm->xpt_ctx->compression_factor = pmap_info.compression_factor;
      gpm->xpt_ctx->radial_compression_factor = pmap_info.radial_compression_factor;
      gpm->xpt_ctx->compress_divertor = pmap_info.compress_divertor;
  }

  gpm->grid = grid;
  gpm->local = local;
  gpm->local_ext = local_ext;
  gpm->global = global;
  gpm->global_ext = global_ext;
  gpm->basis = basis;
  gpm->cdim = grid.ndim; 
  gpm->mc2nu = gkyl_array_new(GKYL_DOUBLE, 3*gpm->basis.num_basis, gpm->local_ext.volume);
  gpm->ref_count = gkyl_ref_count_init(gkyl_position_map_free);

  struct gkyl_position_map *gpm_out = gpm;
  return gpm_out;
}

void
gkyl_position_map_set_mc2nu(struct gkyl_position_map* gpm, struct gkyl_array* mc2nu)
{
  gkyl_array_copy(gpm->mc2nu, mc2nu);
}

void
gkyl_position_map_set_bmag(struct gkyl_position_map* gpm, struct gkyl_comm* comm,
  struct gkyl_array* bmag)
{
  if (gpm->to_optimize) {
    gkyl_free(gpm->constB_ctx->theta_extrema);
    gkyl_free(gpm->constB_ctx->bmag_extrema);
    gkyl_free(gpm->constB_ctx->min_or_max);
  }
  gpm->to_optimize = true;
  // The extrema search samples 2*N boundaries, plus both endpoints.
  int N_boundaries = 2*gpm->constB_ctx->N_theta_boundaries + 1;
  gpm->constB_ctx->theta_extrema = gkyl_malloc(sizeof(double) * N_boundaries);
  gpm->constB_ctx->bmag_extrema = gkyl_malloc(sizeof(double) * N_boundaries);
  gpm->constB_ctx->min_or_max = gkyl_malloc(sizeof(bool) * N_boundaries);
  if (comm == NULL) {
    gkyl_array_release(gpm->bmag_ctx->bmag);
    gpm->bmag_ctx->bmag = gkyl_array_acquire(bmag);
    return;
  }
  else {
    gkyl_comm_array_allgather_host(comm, &gpm->local, \
    &gpm->global, bmag, (struct gkyl_array*) gpm->bmag_ctx->bmag);
  }
}

void
gkyl_position_map_set_compression(struct gkyl_position_map* gpm, double zcut, double zcenter, double w, double psisep)
{
  gpm->xpt_ctx->zcut = zcut;
  gpm->xpt_ctx->zcenter = zcenter;
  gpm->xpt_ctx->w = w;
  gpm->xpt_ctx->psisep = psisep;

  if (gpm->xpt_ctx->radial_compression_factor!=0.0) {
    gpm->maps[0] = position_map_sep_compression;
    gpm->map_derivs[0] = position_map_deriv_sep_compression;
    gpm->ctxs[0] = gpm->xpt_ctx;
  }
  else {
    gpm->maps[0] = gpm->xpt_ctx->maps_backup[0];
    gpm->map_derivs[0] = gpm->xpt_ctx->map_derivs_backup[0];
    gpm->ctxs[0] = gpm->xpt_ctx->ctxs_backup[0];
  }

  gpm->maps[1] = gpm->xpt_ctx->maps_backup[1];
  gpm->map_derivs[1] = gpm->xpt_ctx->map_derivs_backup[1];
  gpm->ctxs[1] = gpm->xpt_ctx->ctxs_backup[1];

  if (gpm->xpt_ctx->compression_factor!=0.0) {
    gpm->maps[2] = position_map_xpt_compression;
    gpm->map_derivs[2] = position_map_deriv_xpt_compression;
    gpm->ctxs[2] = gpm->xpt_ctx;
  }
  else {
    gpm->maps[2] = gpm->xpt_ctx->maps_backup[2];
    gpm->map_derivs[2] = gpm->xpt_ctx->map_derivs_backup[2];
    gpm->ctxs[2] = gpm->xpt_ctx->ctxs_backup[2];
  }
}

void 
gkyl_position_map_eval_mc2nu(const struct gkyl_position_map* gpm, const double *x_comp, double *x_fa)
{
  int cidx[GKYL_MAX_CDIM];
  for(int i = 0; i < gpm->grid.ndim; i++){
    int idxtemp = gpm->global.lower[i] + (int) floor((x_comp[i] - (gpm->grid.lower[i]) )/gpm->grid.dx[i]);
    idxtemp = GKYL_MAX2(GKYL_MIN2(idxtemp, gpm->local.upper[i]), gpm->local.lower[i]);
    cidx[i] = idxtemp;
  }
  long lidx = gkyl_range_idx(&gpm->local, cidx);
  const double *pmap_coeffs = gkyl_array_cfetch(gpm->mc2nu, lidx);
  double cxc[gpm->grid.ndim];
  double x_log[gpm->grid.ndim];
  gkyl_rect_grid_cell_center(&gpm->grid, cidx, cxc);
  for(int i = 0; i < gpm->grid.ndim; i++){
    x_log[i] = (x_comp[i]-cxc[i])/(gpm->grid.dx[i]*0.5);
  }
  double xyz_fa[3];
  for(int i = 0; i < 3; i++){
    xyz_fa[i] = gpm->basis.eval_expand(x_log, &pmap_coeffs[i*gpm->basis.num_basis]);
  }
  for (int i=0; i<gpm->grid.ndim; i++) {
    x_fa[i] = xyz_fa[i];
  }
  x_fa[gpm->grid.ndim-1] = xyz_fa[2];
}

void
gkyl_position_map_optimize(struct gkyl_position_map* gpm, struct gkyl_rect_grid grid,
  struct gkyl_range global)
{
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  gpm->constB_ctx->psi_max   = grid.upper[PSI_IDX];
  gpm->constB_ctx->psi_min   = grid.lower[PSI_IDX];
  gpm->constB_ctx->alpha_max = grid.upper[AL_IDX];
  gpm->constB_ctx->alpha_min = grid.lower[AL_IDX];
  gpm->constB_ctx->theta_max = grid.upper[TH_IDX];
  gpm->constB_ctx->theta_min = grid.lower[TH_IDX];
  gpm->constB_ctx->N_theta_boundaries = global.upper[TH_IDX] - global.lower[TH_IDX] + 2;

  if (gpm->id == GKYL_PMAP_CONSTANT_DB_POLYNOMIAL && gpm->to_optimize == true)
  {
    double psi_center = 0.5 * (gpm->constB_ctx->psi_min + gpm->constB_ctx->psi_max);
    double alpha_center = 0.5 * (gpm->constB_ctx->alpha_min + gpm->constB_ctx->alpha_max);

    gpm->maps[0] = gpm->constB_ctx->maps_backup[0];
    gpm->ctxs[0] = gpm->constB_ctx->ctxs_backup[0];
    gpm->maps[1] = gpm->constB_ctx->maps_backup[1];
    gpm->ctxs[1] = gpm->constB_ctx->ctxs_backup[1];
    gpm->maps[2] = position_map_constB_z_polynomial;
    gpm->ctxs[2] = gpm->constB_ctx;

    gpm->bmag_ctx->crange_global = &gpm->global;
    gpm->bmag_ctx->cbasis = &gpm->basis;
    gpm->bmag_ctx->cgrid = &gpm->grid;

    gpm->constB_ctx->psi    = psi_center;
    gpm->constB_ctx->alpha  = alpha_center;

    calculate_mirror_throat_location_polynomial(gpm->constB_ctx, gpm->bmag_ctx);
    calculate_optimal_mapping_polynomial(gpm->constB_ctx, gpm->bmag_ctx);
  }
  else if (gpm->id == GKYL_PMAP_CONSTANT_DB_NUMERIC && gpm->to_optimize == true)
  {
    double psi_center = 0.5 * (gpm->constB_ctx->psi_min + gpm->constB_ctx->psi_max);
    double alpha_center = 0.5 * (gpm->constB_ctx->alpha_min + gpm->constB_ctx->alpha_max);

    gpm->maps[0] = gpm->constB_ctx->maps_backup[0];
    gpm->ctxs[0] = gpm->constB_ctx->ctxs_backup[0];
    gpm->maps[1] = gpm->constB_ctx->maps_backup[1];
    gpm->ctxs[1] = gpm->constB_ctx->ctxs_backup[1];
    gpm->maps[2] = position_map_constB_z_numeric_moving_average;
    gpm->ctxs[2] = gpm;

    gpm->bmag_ctx->crange_global = &gpm->global;
    gpm->bmag_ctx->cbasis        = &gpm->basis;
    gpm->bmag_ctx->cgrid         = &gpm->grid;

    gpm->constB_ctx->psi    = psi_center;
    gpm->constB_ctx->alpha  = alpha_center;

    find_B_field_extrema(gpm);
    refine_B_field_extrema(gpm);
  }
}

double
gkyl_position_map_slope(const struct gkyl_position_map* gpm, int ix_map,
  double x, double dx, double lower, double upper)
{
  if (gpm->use_map_derivs && gpm->map_derivs[ix_map])
  {
    double slope;
    gpm->map_derivs[ix_map](0.0, &x, &slope, gpm->ctxs[ix_map]);
    return slope;
  }

  assert(dx > 0.0 && upper > lower);
  // Permit roundoff from coordinate construction, but never sample outside
  // the supported domain (including discarded stencil points).
  double tol = 16.0*DBL_EPSILON*fmax(1.0, fmax(fabs(lower), fabs(upper)));
  assert(x >= lower-tol && x <= upper+tol);
  x = fmin(upper, fmax(lower, x));
  double f_left, f_right;
  if (x-dx >= lower && x+dx <= upper) {
    double x_left = x-dx, x_right = x+dx;
    gpm->maps[ix_map](0.0, &x_left, &f_left, gpm->ctxs[ix_map]);
    gpm->maps[ix_map](0.0, &x_right, &f_right, gpm->ctxs[ix_map]);
    return (f_right-f_left)/(2.0*dx);
  }

  double step = upper-x >= x-lower ? fmin(dx, 0.5*(upper-x)) : -fmin(dx, 0.5*(x-lower));
  double x_one = fmin(upper, fmax(lower, x+step));
  double x_two = fmin(upper, fmax(lower, x+2.0*step));
  double f;
  gpm->maps[ix_map](0.0, &x, &f, gpm->ctxs[ix_map]);
  gpm->maps[ix_map](0.0, &x_one, &f_left, gpm->ctxs[ix_map]);
  gpm->maps[ix_map](0.0, &x_two, &f_right, gpm->ctxs[ix_map]);
  return (4.0*(f_left-f)-(f_right-f))/(2.0*step);
}

struct gkyl_position_map*
gkyl_position_map_acquire(const struct gkyl_position_map* gpm)
{
  gkyl_ref_count_inc(&gpm->ref_count);
  return (struct gkyl_position_map*) gpm;
}

void
gkyl_position_map_release(const struct gkyl_position_map *gpm)
{
  gkyl_ref_count_dec(&gpm->ref_count);
}

void
gkyl_position_map_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_position_map *gpm = container_of(ref, struct gkyl_position_map, ref_count);
  gkyl_array_release(gpm->mc2nu);
  gkyl_array_release(gpm->bmag_ctx->bmag);
  if (gpm->to_optimize == true)
  {
    gkyl_free(gpm->constB_ctx->theta_extrema);
    gkyl_free(gpm->constB_ctx->bmag_extrema);
    gkyl_free(gpm->constB_ctx->min_or_max);
  }
  gkyl_free(gpm->bmag_ctx);
  gkyl_free(gpm->constB_ctx);
  gkyl_free(gpm->xpt_ctx);
  gkyl_free(gpm);
}
