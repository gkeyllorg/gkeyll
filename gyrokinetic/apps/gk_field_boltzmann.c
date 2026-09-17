#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_position_map.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_comm_io.h>

#include <assert.h>
#include <float.h>
#include <time.h>

static void gk_field_accumulate_rho_c_boltzmann(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gk_species *s, struct gkyl_array **bflux)
{
  if (s->info.charge > 0.0) {
    // For Boltzmann electrons, we only need ion density (and the ion density
    // times the conf-space Jacobian), not charge density.
    // Rescale moment by inverse of Jacobian.
    gkyl_dg_div_op_range(s->m0.mem_geo, &app->basis, 0, field->rho_c, 0, s->m0.marr, 0,
                         app->gk_geom->geo_int.jacobgeo, &app->local);

    // We also need the M0 flux of the boundary flux through the z
    // boundaries. Put it in the ghost cells of s->m0.marr.
    gk_species_bflux_get_flux_mom(&s->bflux, app->cdim-1, GKYL_LOWER_EDGE,
      GKYL_F_MOMENT_M0, bflux, s->m0.marr, &app->local_lower_ghost[app->cdim-1]);
    gk_species_bflux_get_flux_mom(&s->bflux, app->cdim-1, GKYL_UPPER_EDGE,
      GKYL_F_MOMENT_M0, bflux, s->m0.marr, &app->local_upper_ghost[app->cdim-1]);
  }
}

static void
gk_field_calc_ambi_pot_sheath_vals(gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  // Note that the M0 moment of boundary fluxes along z should
  // be stored in the ghost cells of m0.marr at this point.
  int idx_par = app->cdim-1;
  int off = 2*idx_par;

  int comm_sz;
  gkyl_comm_get_size(app->comm, &comm_sz);

  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];

    if (s->info.charge > 0.0) {
      // Assumes symmetric sheath BCs for now only in 1D
      // NOTE: this relies on the accumulate_rho_c calling gk_species_moment_calc(s->m0)
      // to calculate the particle flux and place it in the ghost cells of s->m0.marr.
      gkyl_ambi_bolt_potential_sheath_calc(field->ambi_pot, GKYL_LOWER_EDGE, 
        &app->local_lower_skin[idx_par], &app->local_lower_ghost[idx_par], app->gk_geom->geo_int.cmag, 
        app->gk_geom->geo_int.jacobtot_inv, s->m0.marr, field->rho_c, s->m0.marr, field->sheath_vals[off]);
      gkyl_ambi_bolt_potential_sheath_calc(field->ambi_pot, GKYL_UPPER_EDGE, 
        &app->local_upper_skin[idx_par], &app->local_upper_ghost[idx_par], app->gk_geom->geo_int.cmag,
        app->gk_geom->geo_int.jacobtot_inv, s->m0.marr, field->rho_c, s->m0.marr, field->sheath_vals[off+1]);

      // Broadcast the sheath values from skin processes to other processes.
      gkyl_comm_array_bcast(app->comm, field->sheath_vals[off]  , field->sheath_vals[off], 0);
      gkyl_comm_array_bcast(app->comm, field->sheath_vals[off+1], field->sheath_vals[off+1], comm_sz-1);

      // Copy upper sheath values into lower ghost & add to lower sheath values for averaging.
      gkyl_array_copy_range_to_range(field->sheath_vals[off+1], field->sheath_vals[off+1],
        &app->local_lower_ghost[idx_par], &app->local_upper_ghost[idx_par]);
      gkyl_array_accumulate(field->sheath_vals[off], 1., field->sheath_vals[off+1]);
      gkyl_array_scale(field->sheath_vals[off], 0.5);
    }
  }
}

static void
gk_field_rhs_boltzmann(struct gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  // Compute sheath density n_i,s and potential phi_s = (Te/e)*ln(n_i,s*v_te/(sqrt(2*pi)*Gamma_i)).
  gk_field_calc_ambi_pot_sheath_vals(app, app->field);

  // Solve phi = phi_s + (Te/e)*ln(n_i/n_i,s).
  gkyl_ambi_bolt_potential_phi_calc(field->ambi_pot, &app->local, &app->local_ext,
    field->rho_c, field->sheath_vals[2*(app->cdim-1)], field->phi_smooth);

  // Smooth the potential along z.
  gk_field_fem_projection_par(app, field, field->phi_smooth, field->phi_smooth);
}

static void
gk_field_fem_release_boltzmann(const gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  gkyl_array_release(f->rho_c);
  gkyl_array_release(f->rho_c_global_dg);
  gkyl_array_release(f->phi_fem);
  gkyl_array_release(f->phi_smooth);

  if (app->use_gpu) {
    gkyl_array_release(f->phi_host);
  }

  gkyl_ambi_bolt_potential_release(f->ambi_pot);
  for (int i = 0; i < 2*app->cdim; ++i) {
    gkyl_array_release(f->sheath_vals[i]);
  }
  gkyl_fem_parproj_release(f->fem_parproj);
  gkyl_array_integrate_release(f->calc_em_energy);
}

void
gk_field_fem_new_boltzmann(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  // Create global subrange we'll copy the field solver solution from (into local).
  gkyl_sub_range_intersect(&f->global_sub_range, &app->global, &app->local);

  // Allocate arrays for charge density.
  f->rho_c = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->rho_c_global_dg = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);

  // Allocate arrays for electrostatic potential.
  f->phi_fem = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  f->phi_smooth = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  // Allocate phi_host for I/O.
  f->phi_host = f->phi_smooth;
  if (app->use_gpu) {
    f->phi_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  }

  double polarization_weight = 1.0;
  f->rhs_phi_func = gk_field_rhs_boltzmann;
  f->accumulate_rhoc_func = gk_field_accumulate_rho_c_boltzmann;

  f->ambi_pot = gkyl_ambi_bolt_potential_new(&app->grid, &app->basis, 
    f->info.electron_mass, f->info.electron_charge, f->info.electron_temp, app->use_gpu);
  
  // Sheath_vals contains both the density and potential sheath values.
  for (int j=0; j<app->cdim; ++j) {
    f->sheath_vals[2*j]   = mkarr(app->use_gpu, 2*app->basis.num_basis, app->local_ext.volume);
    f->sheath_vals[2*j+1] = mkarr(app->use_gpu, 2*app->basis.num_basis, app->local_ext.volume);
  }

  enum gkyl_fem_parproj_bc_type fem_parproj_bc = GKYL_FEM_PARPROJ_NONE;
  for (int d=0; d<app->num_periodic_dir; ++d) {
    if (app->periodic_dirs[d] == app->cdim-1) {
      fem_parproj_bc = GKYL_FEM_PARPROJ_PERIODIC;
    }
  }

  f->fem_parproj = gkyl_fem_parproj_new(&app->global, &app->grid, &app->basis,
    fem_parproj_bc, 0, 0, 0, app->use_gpu);

  if (app->cdim == 1) {
    f->es_energy_fac_1d = polarization_weight;
    f->calc_em_energy = gkyl_array_integrate_new(&app->grid, &app->basis, 
      1, GKYL_ARRAY_INTEGRATE_OP_SQ, app->use_gpu);
  } else {
    f->calc_em_energy = gkyl_array_integrate_new(&app->grid, &app->basis, 
      1, GKYL_ARRAY_INTEGRATE_OP_EPS_GRADPERP_SQ, app->use_gpu);
  }

  // Create operator needed for FLR effects.
  f->use_flr = false;
  f->invert_flr = gk_field_invert_flr_none;
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    if (s->info.flr.type) {
      f->use_flr = f->use_flr || s->info.flr.type;
    }
  }

  f->release_func = gk_field_fem_release_boltzmann;
}
