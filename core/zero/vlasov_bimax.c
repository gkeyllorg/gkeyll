#include <math.h>
#include <stddef.h>
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_dist_type.h>
#include <gkyl_dist_type_priv.h>


struct gkyl_vlasov_bimax {
  struct gkyl_dist_proj_type dist;

  enum gkyl_model_id model_id;
  int vdim;
};

/**
 * Free the BiMax distribution projection object.
 */
static void
vlasov_bimax_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_vlasov_bimax *bimax =
    (struct gkyl_vlasov_bimax *)
    ((char *) ref- offsetof(struct gkyl_vlasov_bimax, dist.ref_count));

  gkyl_free(bimax);
}


/**
 * Non-relativistic bi-Maxwellian distribution.
 *
 * moms layout:
 *
 *   moms[0]              = n
 *   moms[1 ... vdim]     = V_drift
 *   moms[vdim+1]         = T_parallel/m
 *   moms[vdim+2]         = T_perpendicular/m
 */
static void
vlasov_bimax_default(
  const struct gkyl_dist_proj_type *dist_proj,
  const double *xc,
  const double *dx,
  const int *idx,
  const double *moms,
  double *out)
{
  const struct gkyl_vlasov_bimax *bimax = (const struct gkyl_vlasov_bimax *) dist_proj;

  int cdim = dist_proj->cdim;
  int vdim = bimax->vdim;

  double n = moms[0];
  const double *V_drift = &moms[1];
  double Tpar_over_m = moms[vdim+1];
  double Tperp_over_m = moms[vdim+2];

  if (n <= 0.0 || Tpar_over_m <= 0.0 || Tperp_over_m <= 0.0) {
    out[0] = 0.0;
    return;
  }

  double vpar = xc[cdim] - V_drift[0];

  double vperp_sq = 0.0;

  for (int d=1; d<vdim; ++d) {
    double v =
      xc[cdim+d] - V_drift[d];

    vperp_sq += v*v;
  }

  double amplitude = n/(sqrt(2.0*GKYL_PI*Tpar_over_m) * pow( 2.0*GKYL_PI*Tperp_over_m, 0.5*(vdim-1)));

  out[0] = amplitude* exp( -vpar*vpar/(2.0*Tpar_over_m) -vperp_sq/(2.0*Tperp_over_m));
}


struct gkyl_dist_proj_type*
gkyl_vlasov_bimax_new(
  int cdim,
  int pdim,
  int poly_order,
  int num_config,
  int num_phase,
  enum gkyl_model_id model_id)
{
  struct gkyl_vlasov_bimax *bimax = gkyl_malloc(sizeof(*bimax));

  bimax->model_id = model_id;
  bimax->vdim = pdim-cdim;

  bimax->dist.dist_id = GKYL_DIST_TYPE_BIMAX;

  bimax->dist.cdim = cdim;

  bimax->dist.pdim = pdim;

  bimax->dist.poly_order = poly_order;

  bimax->dist.num_config =  num_config;

  bimax->dist.num_phase = num_phase;

  // n + V_drift[vdim] + T_parallel/m + T_perpendicular/m
  bimax->dist.num_mom = bimax->vdim + 3;

  switch (model_id) {
    default:
      bimax->dist.proj_kernel = vlasov_bimax_default;
      break;
  }

  bimax->dist.ref_count = gkyl_ref_count_init(vlasov_bimax_free);

  bimax->dist.flags = 0;
  bimax->dist.on_dev = &bimax->dist;

  return &bimax->dist;
}