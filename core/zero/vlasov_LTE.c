#include <math.h>
#include <stddef.h>
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_dist_type.h>
#include <gkyl_dist_type_priv.h>


struct gkyl_vlasov_lte {
  struct gkyl_dist_proj_type dist;

  enum gkyl_model_id model_id;
  int vdim;
};

/**
 * Free the LTE distribution projection object.
 */
static void
vlasov_lte_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_vlasov_lte *lte =
    (struct gkyl_vlasov_lte *)
    ((char *) ref - offsetof(struct gkyl_vlasov_lte, dist.ref_count));

  gkyl_free(lte);
}


/**
 * Non-relativistic LTE distribution.
 *
 * moms layout:
 *
 *   moms[0]              = n
 *   moms[1 ... vdim]     = V_drift
 *   moms[vdim+1]         = T/m
 */
static void
vlasov_lte_default(
  const struct gkyl_dist_proj_type *dist_proj,
  const double *xc,
  const double *dx,
  const int *idx,
  const double *moms,
  double *out)

{
  const struct gkyl_vlasov_lte *lte =
    (const struct gkyl_vlasov_lte *) dist_proj;

  int cdim = dist_proj->cdim;
  int vdim = lte->vdim;

  double n = moms[0];
  const double *V_drift = &moms[1];
  double T_over_m = moms[vdim+1];

  if (n <= 0.0 || T_over_m <= 0.0) {
    out[0] = 0.0;
    return;
  }

  double efact = 0.0;

  for (int d=0; d<vdim; ++d) {
    double v = xc[cdim+d] - V_drift[d];

    efact += v*v;
  }

  double amplitude = n/sqrt(pow(2.0*GKYL_PI*T_over_m, vdim));

  out[0] = amplitude * exp(-efact/(2.0*T_over_m));
}


/**
 * Special-relativistic LTE distribution.
 *
 * Evaluates the Maxwell-Juttner distribution using the stationary-frame
 * density, the bulk four-velocity, and T/m.
 *
 * moms layout:
 *
 *   moms[0]              = n
 *   moms[1 ... vdim]     = spatial bulk four-velocity
 *   moms[vdim+1]         = T/m
 */
static void
vlasov_lte_sr(
  const struct gkyl_dist_proj_type *dist_proj,
  const double *xc,
  const double *dx,
  const int *idx,
  const double *moms,
  double *out)
{
  const struct gkyl_vlasov_lte *lte =
    (const struct gkyl_vlasov_lte *) dist_proj;

  int cdim = dist_proj->cdim;
  int vdim = lte->vdim;

  double n = moms[0];
  const double *V_drift = &moms[1];
  double T_over_m = moms[vdim+1];

  if (n <= 0.0 || T_over_m <= 0.0) {
    out[0] = 0.0;
    return;
  }

  double vv = 0.0;
  double vu = 0.0;
  double uu = 0.0;

  for (int d=0; d<vdim; ++d) {
    double u = V_drift[d];
    double p = xc[cdim+d];

    vv += u*u;
    vu += u*p;
    uu += p*p;
  }

  double GammaV =
    sqrt(1.0 + vv);

  double amplitude = n * (1.0/(4.0*GKYL_PI*T_over_m)) * sqrt(2.0*T_over_m/GKYL_PI);

  out[0] = amplitude * exp((1.0/T_over_m) - (1.0/T_over_m) * (GammaV*sqrt(1.0+uu) - vu));
}


/**
 * Canonical-PB LTE distribution.
 *
 * Evaluates the LTE distribution in canonical momentum coordinates.
 *
 * moms layout:
 *
 *   moms[0]              = n
 *   moms[1 ... vdim]     = V_drift
 *   moms[vdim+1]         = T/m
 *
 * The required geometric quantities are not yet available through the
 * generic dist_proj_t interface.
 */
static void
vlasov_lte_canonical_pb(
  const struct gkyl_dist_proj_type *dist_proj,
  const double *xc,
  const double *dx,
  const int *idx,
  const double *moms,
  double *out)
{
  assert(false && "Canonical-PB LTE missing params");
}


struct gkyl_dist_proj_type*
gkyl_vlasov_lte_new(
  int cdim,
  int pdim,
  int poly_order,
  int num_config,
  int num_phase,
  enum gkyl_model_id model_id)
{
  struct gkyl_vlasov_lte *lte = gkyl_malloc(sizeof(*lte));

  lte->model_id = model_id;
  lte->vdim = pdim-cdim;

  lte->dist.dist_id = GKYL_DIST_TYPE_LTE;

  lte->dist.cdim = cdim;

  lte->dist.pdim = pdim;

  lte->dist.poly_order = poly_order;

  lte->dist.num_config = num_config;

  lte->dist.num_phase = num_phase;

  // n + V_drift[vdim] + T/m
  lte->dist.num_mom = lte->vdim + 2;

  switch (model_id) {
    case GKYL_MODEL_SR:
      lte->dist.proj_kernel = vlasov_lte_sr;
      break;

    case GKYL_MODEL_CANONICAL_PB:
    case GKYL_MODEL_CANONICAL_PB_GR:
      lte->dist.proj_kernel = vlasov_lte_canonical_pb;
      break;

    default:
      lte->dist.proj_kernel = vlasov_lte_default;
      break;
  }

  lte->dist.ref_count = gkyl_ref_count_init(vlasov_lte_free);

  lte->dist.flags = 0;
  lte->dist.on_dev = &lte->dist;

  return &lte->dist;
}