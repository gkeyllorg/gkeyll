#include <gkyl_alloc.h>
#include <gkyl_app.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_dg_interpolate.h>
#include <gkyl_translate_dim.h>

#include <assert.h>
#include <time.h>

/**
 * Initialize fluid neutral species.
 * 
 * @param gk Input parameters from input files.
 * @param app Gyrokinetic app.
 * @param ns Neutral species to initialize.
 */ 
void
gk_neut_species_fluid_init(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns);

/**
 * Initialize kinetic neutral species.
 * 
 * @param gk Input parameters from input files.
 * @param app Gyrokinetic app.
 * @param ns Neutral species to initialize.
 */ 
void
gk_neut_species_kinetic_init(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns);

/**
 *
 * Methods common to fluid and kinetic neutrals.
 *
 */

double
gk_neut_species_rhs_static(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms);

double
gk_neut_species_rhs_implicit_static(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt);

void
gk_neut_species_apply_bc_static(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, struct gkyl_array *f);

void
gk_neut_species_step_f_dynamic(struct gkyl_array* out, double dt,
  const struct gkyl_array* inp);

void
gk_neut_species_step_f_static(struct gkyl_array* out, double dt,
  const struct gkyl_array* inp);

void
gk_neut_species_combine_dynamic(struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng);

void
gk_neut_species_combine_static(struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2,
  const struct gkyl_range *rng);

void
gk_neut_species_copy_range_dynamic(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range);

void
gk_neut_species_copy_range_static(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range);

void
gk_neut_species_write_dynamic(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);

void
gk_neut_species_write_mom_dynamic(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);

void
gk_neut_species_calc_integrated_mom_dynamic(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm);

void
gk_neut_species_write_integrated_mom_dynamic(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns);

void
gk_neut_species_write_static(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);

void
gk_neut_species_write_init_only(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);

void
gk_neut_species_write_mom_static(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);

void
gk_neut_species_write_mom_init_only(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm, int frame);

void
gk_neut_species_calc_integrated_mom_static(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, double tm);

void
gk_neut_species_write_integrated_mom_static(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns);

void
gk_neut_species_n_iter_corr_enabled(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns, int sidx);

void
gk_neut_species_n_iter_corr_disabled(gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns, int sidx);

void
gk_neut_species_release_static(const gkyl_gyrokinetic_app* app, const struct gk_neut_species *s);

