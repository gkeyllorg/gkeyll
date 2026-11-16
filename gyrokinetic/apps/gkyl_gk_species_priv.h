#pragma once

#include <gkyl_alloc.h>
#include <gkyl_app.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gyrokinetic_pol_density.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_dg_interpolate.h>
#include <gkyl_translate_dim.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_gk_maxwellian_moments.h>

#include <assert.h>
#include <time.h>

void gk_species_gyroaverage_disabled(
  gkyl_gyrokinetic_app *app, struct gk_species *species, struct gkyl_array *field_in,
  struct gkyl_array *field_gyroavg
);

void gk_species_gyroaverage_enabled(
  gkyl_gyrokinetic_app *app, struct gk_species *species, struct gkyl_array *field_in,
  struct gkyl_array *field_gyroavg
);

double gk_species_rhs_static(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  const struct gkyl_array *fbar_in, struct gkyl_array *rhs, struct gkyl_array **bflux_moms
);

double gk_species_rhs_implicit_static(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt
);

void gk_species_apply_bc_static(
  gkyl_gyrokinetic_app *app, struct gk_species *species, double tm, struct gkyl_array *f
);

void gk_species_step_f_dynamic(struct gkyl_array *out, double dt, const struct gkyl_array *inp);

void gk_species_step_f_static(struct gkyl_array *out, double dt, const struct gkyl_array *inp);

void gk_species_combine_dynamic(
  struct gkyl_array *out, double c1, const struct gkyl_array *arr1, double c2,
  const struct gkyl_array *arr2, const struct gkyl_range *rng
);

void gk_species_combine_static(
  struct gkyl_array *out, double c1, const struct gkyl_array *arr1, double c2,
  const struct gkyl_array *arr2, const struct gkyl_range *rng
);

void gk_species_copy_range_dynamic(
  struct gkyl_array *out, const struct gkyl_array *inp, const struct gkyl_range *range
);

void gk_species_copy_range_static(
  struct gkyl_array *out, const struct gkyl_array *inp, const struct gkyl_range *range
);

void gk_species_write_vcutsq(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gkyl_bc_sheath_gyrokinetic *bc_sheath,
  const char *dir_str, const char *edge_str, double tm, int frame
);

void gk_species_write_dynamic(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame
);

void gk_species_write_static(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame
);

void gk_species_write_cfl_enabled(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame
);

void gk_species_write_cfl_disabled(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame
);

void gk_species_write_mom_dynamic(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame
);

void gk_species_write_mom_static(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame
);

void gk_species_calc_int_mom_dt_enabled(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, double dt, struct gkyl_array *fdot_int_mom,
  struct gkyl_array *fbardot_int_mom
);

void gk_species_calc_int_mom_dt_disabled(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, double dt, struct gkyl_array *fdot_int_mom,
  struct gkyl_array *fbardot_int_mom
);

void gk_species_calc_int_mom_dt(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, double dt, struct gkyl_array *fdot_int_mom,
  struct gkyl_array *fbardot_int_mom
);

void gk_species_calc_integrated_mom_dynamic(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm
);

void gk_species_calc_integrated_mom_static(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm
);

void gk_species_write_integrated_mom_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *gks);

void gk_species_write_integrated_mom_static(gkyl_gyrokinetic_app *app, struct gk_species *gks);

void gk_species_calc_L2norm_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm);

void gk_species_calc_L2norm_static(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm);

void gk_species_write_L2norm_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *gks);

void gk_species_write_L2norm_static(gkyl_gyrokinetic_app *app, struct gk_species *gks);

void gk_species_release_static(const gkyl_gyrokinetic_app *app, const struct gk_species *gks);

void gk_species_init_static(
  struct gkyl_gk *gk_app_inp, struct gkyl_gyrokinetic_app *app, struct gk_species *gks
);

void gk_species_apply_ic(gkyl_gyrokinetic_app *app, struct gk_species *gks, double t0);

/**
 * Import an initial state during species construction.
 * @param app Gyrokinetic app.
 * @param gks Initialized species state and moment objects.
 * @param inp Initial-condition import options.
 */
void gk_species_file_import_init(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gkyl_gyrokinetic_ic_import inp
);

/**
 * Import a kinetic initial state, including requested donor-grid transformations.
 * @param app Gyrokinetic app.
 * @param gks Kinetic species receiving the imported distribution.
 * @param inp Initial-condition import options.
 */
void gk_species_kinetic_file_import_init(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gkyl_gyrokinetic_ic_import inp
);

/**
 * Import an unmodified, matching-grid conservative fluid state.
 * @param app Gyrokinetic app.
 * @param gks Fluid species receiving the imported moments.
 * @param inp Initial-condition import options.
 */
void gk_species_fluid_file_import_init(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gkyl_gyrokinetic_ic_import inp
);

void gk_species_apply_ic_cross(gkyl_gyrokinetic_app *app, struct gk_species *gks_self, double t0);

double gk_species_rhs(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  const struct gkyl_array *fbar_in, struct gkyl_array *rhs, struct gkyl_array **bflux_moms
);

double gk_species_rhs_implicit(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt
);

void gk_species_step_f(
  struct gk_species *species, struct gkyl_array *out, double a, const struct gkyl_array *inp
);

void gk_species_combine(
  struct gk_species *species, struct gkyl_array *out, double c1, const struct gkyl_array *arr1,
  double c2, const struct gkyl_array *arr2, const struct gkyl_range *rng
);

void gk_species_copy_range(
  struct gk_species *species, struct gkyl_array *out, const struct gkyl_array *inp,
  const struct gkyl_range *range
);

void gk_species_apply_bc(
  gkyl_gyrokinetic_app *app, struct gk_species *species, double tm, struct gkyl_array *f
);

void gk_species_n_iter_corr(gkyl_gyrokinetic_app *app, const struct gk_species *gks, int sidx);

void gk_species_write(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame);

void gk_species_write_mom(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame);

void gk_species_calc_integrated_mom(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm);

void gk_species_write_integrated_mom(gkyl_gyrokinetic_app *app, struct gk_species *gks);

void gk_species_calc_L2norm(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm);

void gk_species_write_L2norm(gkyl_gyrokinetic_app *app, struct gk_species *gks);

void gk_species_init_io(struct gkyl_gyrokinetic_app *app, struct gk_species *gks);

void gk_species_init_integrated_mom(struct gkyl_gyrokinetic_app *app, struct gk_species *gks);

void gk_species_release_integrated_mom(
  const struct gkyl_gyrokinetic_app *app, const struct gk_species *gks
);

void gk_species_kinetic_init(struct gkyl_gk *, struct gkyl_gyrokinetic_app *, struct gk_species *);
void gk_species_fluid_init(struct gkyl_gk *, struct gkyl_gyrokinetic_app *, struct gk_species *);
void gk_species_kinetic_release(const gkyl_gyrokinetic_app *, const struct gk_species *);
