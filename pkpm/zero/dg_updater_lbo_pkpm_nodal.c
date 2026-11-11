#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_updater_lbo_pkpm_nodal.h>
#include <gkyl_hyper_dg.h>

struct gkyl_dg_updater_lbo_pkpm_nodal {
  struct gkyl_dg_eqn *eqn; // the operator
  struct gkyl_hyper_dg *slvr; // solver over the velocity directions
};

struct gkyl_dg_updater_lbo_pkpm_nodal *
gkyl_dg_updater_lbo_pkpm_nodal_new(
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_basis *conf_basis,
  const struct gkyl_basis *phase_basis, const struct gkyl_range *conf_range,
  struct gkyl_dg_lbo_pkpm_nodal_auxfields *auxin, bool use_gpu
)
{
  struct gkyl_dg_updater_lbo_pkpm_nodal *up = gkyl_malloc(sizeof(*up));

  up->eqn = gkyl_dg_lbo_pkpm_nodal_new(conf_basis, phase_basis, conf_range, phase_grid, use_gpu);
  gkyl_lbo_pkpm_nodal_set_auxfields(up->eqn, *auxin);

  int cdim = conf_basis->ndim, pdim = phase_basis->ndim;
  int vdim = pdim - cdim;
  int num_up_dirs = vdim;
  int up_dirs[GKYL_MAX_DIM] = {0};
  int zero_flux_flags[2 * GKYL_MAX_DIM] = {0};
  for (int d = 0; d < vdim; ++d) {
    up_dirs[d] = cdim + d;
    zero_flux_flags[cdim + d] = zero_flux_flags[cdim + d + pdim] = 1;
  }

  up->slvr = gkyl_hyper_dg_new(
    phase_grid, phase_basis, up->eqn, num_up_dirs, up_dirs, zero_flux_flags, 1, use_gpu
  );
  return up;
}

void
gkyl_dg_updater_lbo_pkpm_nodal_advance(
  struct gkyl_dg_updater_lbo_pkpm_nodal *up, const struct gkyl_range *update_rng,
  const struct gkyl_array *GKYL_RESTRICT fnodal, struct gkyl_array *GKYL_RESTRICT cflrate,
  struct gkyl_array *GKYL_RESTRICT rhs_nodal
)
{
  gkyl_hyper_dg_advance(up->slvr, update_rng, fnodal, cflrate, rhs_nodal);
}

void
gkyl_dg_updater_lbo_pkpm_nodal_release(struct gkyl_dg_updater_lbo_pkpm_nodal *up)
{
  gkyl_dg_eqn_release(up->eqn);
  gkyl_hyper_dg_release(up->slvr);
  gkyl_free(up);
}
