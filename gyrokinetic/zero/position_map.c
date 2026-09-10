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
#include <stddef.h>

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
  struct gkyl_position_map *gpm = gkyl_malloc(sizeof(*gpm));
  gpm->id = GKYL_PMAP_USER_INPUT;

  gpm->use_map_derivs = true;
  gpm->to_optimize = false;
  gpm->mc2nu = gkyl_array_new(GKYL_DOUBLE, 1, 1);
  gpm->constB_ctx = gkyl_malloc(sizeof(struct gkyl_position_map_const_B_ctx));
  gpm->xpt_ctx = gkyl_malloc(sizeof(struct gkyl_position_map_xpt_ctx));
  gpm->bmag_ctx = gkyl_malloc(sizeof(struct gkyl_bmag_ctx));
  gpm->bmag_ctx->bmag = gkyl_array_new(GKYL_DOUBLE, 1, 1);
  gpm->mc2nu_dev = 0;
  gpm->on_dev = gpm; // On the CPU the object points to itself.
  gpm->ref_count = gkyl_ref_count_init(gkyl_position_map_free);
  
  for (int i = 0; i < 3; i++){
    gpm->maps[i] = gkyl_position_map_identity;
    gpm->map_derivs[i] = gkyl_position_map_identity_slope;
    gpm->ctxs[i] = 0;
    gpm->constB_ctx->maps_backup[i] = gkyl_position_map_identity;
    gpm->constB_ctx->ctxs_backup[i] = 0;
    gpm->xpt_ctx->maps_backup[i] = gkyl_position_map_identity;
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
  struct gkyl_position_map *gpm = gkyl_malloc(sizeof(*gpm));
  gpm->id = pmap_info.id;

  gpm->bmag_ctx = gkyl_malloc(sizeof(struct gkyl_bmag_ctx));
  gpm->bmag_ctx->bmag = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, global_ext.volume);
  gpm->to_optimize = false;
  gpm->use_map_derivs = (pmap_info.id == GKYL_PMAP_XPT_COMPRESSION || pmap_info.id == GKYL_PMAP_USER_INPUT_W_DERIVATIVE) ? true : false;


  gpm->constB_ctx = gkyl_malloc(sizeof(struct gkyl_position_map_const_B_ctx));
  gpm->xpt_ctx = gkyl_malloc(sizeof(struct gkyl_position_map_xpt_ctx));

  for (int i = 0; i < 3; i++){
    gpm->maps[i] = gkyl_position_map_identity;
    gpm->map_derivs[i] = gkyl_position_map_identity_slope;
    gpm->ctxs[i] = 0;
    gpm->constB_ctx->maps_backup[i] = gkyl_position_map_identity;
    gpm->constB_ctx->ctxs_backup[i] = 0;
    gpm->xpt_ctx->maps_backup[i] = gkyl_position_map_identity;
    gpm->xpt_ctx->ctxs_backup[i] = 0;
  }

  switch (pmap_info.id)
  {
    case GKYL_PMAP_USER_INPUT:
      for (int i = 0; i < 3; i++){
        if (pmap_info.maps[i] != 0)
        { gpm->maps[i] = pmap_info.maps[i];
          gpm->ctxs[i] = pmap_info.ctxs[i];
        }
      }

    case GKYL_PMAP_USER_INPUT_W_DERIVATIVE:
      for (int i = 0; i < 3; i++){
        if (pmap_info.maps[i] != 0)
        { gpm->maps[i] = pmap_info.maps[i];
          gpm->map_derivs[i] = pmap_info.map_derivs[i];
          gpm->ctxs[i] = pmap_info.ctxs[i];
        }
      }

    case GKYL_PMAP_CONSTANT_DB_POLYNOMIAL:

      for (int i = 0; i < 2; i++){
        if (pmap_info.maps[i] != 0)
        { gpm->constB_ctx->maps_backup[i] = pmap_info.maps[i];
          gpm->constB_ctx->ctxs_backup[i] = pmap_info.ctxs[i];
        }
      }
      gpm->constB_ctx->map_strength = pmap_info.map_strength;

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

    case GKYL_PMAP_XPT_COMPRESSION:

      for (int i = 0; i < 2; i++){
        if (pmap_info.maps[i] != 0)
        { gpm->xpt_ctx->maps_backup[i] = pmap_info.maps[i];
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
  gpm->mc2nu_dev = 0;
  gpm->on_dev = gpm; // On the CPU the object points to itself.
  gpm->ref_count = gkyl_ref_count_init(gkyl_position_map_free);

  struct gkyl_position_map *gpm_out = gpm;
  return gpm_out;
}

void
gkyl_position_map_set_mc2nu(struct gkyl_position_map* gpm, struct gkyl_array* mc2nu)
{
  gkyl_array_copy(gpm->mc2nu, mc2nu);
  if (gpm->mc2nu_dev)
    gkyl_array_copy(gpm->mc2nu_dev, gpm->mc2nu);
}

void
gkyl_position_map_set_bmag(struct gkyl_position_map* gpm, struct gkyl_comm* comm,
  struct gkyl_array* bmag)
{
  gpm->to_optimize = true;
  int N_boundaries = gpm->constB_ctx->N_theta_boundaries;
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
    gpm->ctxs[0] = gpm->xpt_ctx->ctxs_backup[0];
  }

  gpm->maps[1] = gpm->xpt_ctx->maps_backup[1];
  gpm->ctxs[1] = gpm->xpt_ctx->ctxs_backup[1];

  if (gpm->xpt_ctx->compression_factor!=0.0) {
    gpm->maps[2] = position_map_xpt_compression;
    gpm->map_derivs[2] = position_map_deriv_xpt_compression;
    gpm->ctxs[2] = gpm->xpt_ctx;
  }
  else {
    gpm->maps[2] = gpm->xpt_ctx->maps_backup[2];
    gpm->ctxs[2] = gpm->xpt_ctx->ctxs_backup[2];
  }
}

void
gkyl_position_map_make_cu_dev(struct gkyl_position_map* gpm)
{
#ifdef GKYL_HAVE_CUDA
  // Device mirror of the mapping's DG coefficients, kept in sync by set_mc2nu.
  gpm->mc2nu_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, gpm->mc2nu->ncomp, gpm->mc2nu->size);
  gkyl_array_copy(gpm->mc2nu_dev, gpm->mc2nu);

  // Clone with device pointers; host-only members are nulled.
  struct gkyl_position_map gpm_ho = *gpm;
  gpm_ho.mc2nu = gpm->mc2nu_dev->on_dev;
  for (int i=0; i<3; i++) {
    gpm_ho.maps[i] = 0;
    gpm_ho.map_derivs[i] = 0;
    gpm_ho.ctxs[i] = 0;
  }
  gpm_ho.bmag_ctx = 0;
  gpm_ho.constB_ctx = 0;
  gpm_ho.xpt_ctx = 0;

  struct gkyl_position_map *gpm_cu = gkyl_cu_malloc(sizeof(struct gkyl_position_map));
  gpm_ho.on_dev = gpm_cu;
  gkyl_cu_memcpy(gpm_cu, &gpm_ho, sizeof(struct gkyl_position_map), GKYL_CU_MEMCPY_H2D);

  // Overwrite the basis in the clone with one whose function pointers are
  // device addresses (the temporary device container can be released after
  // copying its contents; the code addresses it holds remain valid).
  struct gkyl_basis *basis_cu;
  switch (gpm->basis.b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      basis_cu = gkyl_cart_modal_serendip_cu_dev_new(gpm->basis.ndim, gpm->basis.poly_order);
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      basis_cu = gkyl_cart_modal_tensor_cu_dev_new(gpm->basis.ndim, gpm->basis.poly_order);
      break;
    default:
      assert(false);
      break;
  }
  gkyl_cu_memcpy((char *)gpm_cu + offsetof(struct gkyl_position_map, basis), basis_cu,
    sizeof(struct gkyl_basis), GKYL_CU_MEMCPY_D2D);
  gkyl_cart_modal_basis_release_cu(basis_cu);

  gpm->on_dev = gpm_cu;
#else
  assert(false);
#endif
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
    double psi_center = pow(0.5 * (sqrt(gpm->constB_ctx->psi_min) + sqrt(gpm->constB_ctx->psi_max)), 2.0);
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
  double x, double dx, int ix_comp, const struct gkyl_range *nrange)
{
  if (gpm->use_map_derivs)
  {
    double slope;
    gpm->map_derivs[ix_map](0.0, &x, &slope, gpm->ctxs[ix_map]);
    return slope;
  }

  double x_left = x - dx;
  double x_right = x + dx;
  double f_left, f, f_right;
  gpm->maps[ix_map](0.0, &x_left, &f_left, gpm->ctxs[ix_map]);
  gpm->maps[ix_map](0.0, &x_right, &f_right, gpm->ctxs[ix_map]);
  double slope;
  if (ix_comp == nrange->lower[ix_map])
  {
    gpm->maps[ix_map](0.0, &x, &f, gpm->ctxs[ix_map]);
    slope = (f_right - f) / dx;
  }
  else if (ix_comp == nrange->upper[ix_map])
  {
    gpm->maps[ix_map](0.0, &x, &f, gpm->ctxs[ix_map]);
    slope = (f - f_left) / dx;
  }
  else
  {
    slope = (f_right - f_left) / (2.0 * dx);
  }
  return slope;
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
  if (gpm->mc2nu_dev)
    gkyl_array_release(gpm->mc2nu_dev);
  if (gpm->on_dev != gpm)
    gkyl_cu_free(gpm->on_dev);
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
