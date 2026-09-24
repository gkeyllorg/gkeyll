#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>

#include <assert.h>
#include <math.h>

void
gk_field_flr_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  assert(app->cdim > 1);
  f->invert_flr = gk_field_invert_flr;

  double flr_weight = 0.0;
  for (int i = 0; i < app->num_species; ++i){
    struct gk_species *s = &app->species[i];
    double gyroradius_bmag = s->info.flr.bmag ? s->info.flr.bmag : app->bmag_ref;
    flr_weight += s->info.flr.Tperp * s->info.mass / (pow(s->info.charge * gyroradius_bmag, 2.0));
  }
  // Initialize the weight in the Laplacian operator.
  f->flr_rhoSq_sum = mkarr(app->use_gpu, (2 * (app->cdim - 1) - 1) * app->basis.num_basis, app->local_ext.volume);
  gkyl_array_set_offset(f->flr_rhoSq_sum, flr_weight, app->gk_geom->geo_int.gxxj, 0 * app->basis.num_basis);
  if (app->cdim > 2) {
    gkyl_array_set_offset(f->flr_rhoSq_sum, flr_weight, app->gk_geom->geo_int.gxyj, 1 * app->basis.num_basis);
    gkyl_array_set_offset(f->flr_rhoSq_sum, flr_weight, app->gk_geom->geo_int.gyyj, 2 * app->basis.num_basis);
  }
  // Initialize the factor multiplying the field in the FLR operator.
  f->flr_kSq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkyl_array_shiftc(f->flr_kSq, -pow(sqrt(2.0), app->cdim), 0); // Sets kSq=-1.

  // If domain is not periodic use Dirichlet BCs.
  struct gkyl_poisson_bc flr_bc = {};

  bool bc_is_np[GKYL_MAX_CDIM]; // Is the BC in this direction non-periodic?
  for (int d = 0; d < app->cdim; ++d) {
    bc_is_np[d] = true;
  }
  for (int d = 0; d < app->num_periodic_dir; ++d) {
    bc_is_np[app->periodic_dirs[d]] = false;
  }

  for (int d = 0; d < app->cdim - 1; d++) {
    if (bc_is_np[d]) {
      struct gkyl_gyrokinetic_bc *bc_lo = gk_fetch_bc_with_dir_edge(f->info.poisson_bcs, 2 * app->cdim, d, GKYL_LOWER_EDGE);
      if (bc_lo != 0) {
        flr_bc.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_DIRICHLET_VARYING);
      }

      struct gkyl_gyrokinetic_bc *bc_up = gk_fetch_bc_with_dir_edge(f->info.poisson_bcs, 2 * app->cdim, d, GKYL_UPPER_EDGE);
      if (bc_up != 0) {
        flr_bc.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_DIRICHLET_VARYING);
      }
    } else {
      flr_bc.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
      flr_bc.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
    }
  }
  // Deflated Poisson solve is performed on range assuming decomposition is *only* in z.
  f->flr_op = gkyl_deflated_fem_poisson_new(app->grid, app->basis_on_dev, app->basis,
                                            app->local, app->local, f->flr_rhoSq_sum, f->flr_kSq, flr_bc, NULL, app->use_gpu);
}

void
gk_field_invert_flr(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi)
{
  gkyl_deflated_fem_poisson_advance(field->flr_op, phi, phi, phi);
}

void
gk_field_invert_flr_none(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi)
{
}

void
gk_field_flr_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  gkyl_array_release(f->flr_rhoSq_sum);
  gkyl_array_release(f->flr_kSq);
  gkyl_deflated_fem_poisson_release(f->flr_op);
}
