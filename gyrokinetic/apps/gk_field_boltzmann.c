#include <gkyl_alloc.h>
#include <gkyl_boltz_pot_0curr.h>
#include <gkyl_array_ops.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>
#include <gkyl_gk_particle_flux_parallel.h>

struct gk_field_boltzmann {
  struct gkyl_boltz_pot_0curr *pot_op; // Boltzmann potential updater.
  struct gkyl_array *nphi_ref; // Reference density and potential in the local lower ghost slab.
  int num_ions; // Number of positive kinetic species contributing to the closure.
  int ion_idx[GKYL_MAX_SPECIES]; // Indices of those species in the app.
  struct gkyl_comm
    *comm_par; // App-owned parallel line communicator; NULL if the closure is disabled.
  int rank_par; // This rank's index in comm_par.
  int boundary_rank[2]; // Lower/upper boundary ranks in comm_par.
  bool is_boundary_rank[2]; // Whether this rank owns each physical parallel boundary.
  struct gkyl_range skin[2]; // Local configuration-space skin ranges at the parallel ends.
  struct gk_species_wall_potential *wall[2]; // Wall profiles used at the lower/upper ends.
  struct gk_species_wall_potential absorbing_wall[2]; // Owned profiles for absorbing ion BCs.
  struct gkyl_gk_particle_flux_parallel
    *m0flux_op[GKYL_MAX_SPECIES][2]; // Particle-flux updaters per ion and end.
  struct gkyl_array *m0_buff; // Ion density after removing the configuration Jacobian.
  struct gkyl_array *m0flux_buff; // Instantaneous particle flux from one ion species.
  struct gkyl_array *m0flux_tot[2]; // Summed outward ion m0 fluxes weighted by q_i/|q_e|.
  struct gkyl_array *end_vals[2]; // Packed surface flux, area factor, density, and wall potential.
};

const struct gkyl_array *
gk_field_boltzmann_nphi_ref(const struct gk_field *field)
{
  return field->boltzmann->nphi_ref;
}

void
gk_field_boltzmann_init(gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  struct gk_field_boltzmann *boltz = field->boltzmann;
  int par_dir = app->cdim - 1;

  // The constructor leaves comm_par unset when the sheath closure is disabled.
  if (!boltz->comm_par) {
    return;
  }

  // Identify ions.
  for (int k = 0; k < field->num_rhoc_species; ++k) {
    int i = field->rhoc_species_idx[k];
    struct gk_species *gks = &app->species[i];
    if (gks->info.charge <= 0.0) {
      continue;
    }
    if (gks->is_fluid) {
      gkyl_exit("The Boltzmann parallel-loss closure requires kinetic ions.");
    }

    int ion = boltz->num_ions++;
    boltz->ion_idx[ion] = i;

    for (int edge = 0; edge < 2; ++edge) {
      const struct gkyl_gyrokinetic_bc *bc = edge ? &gks->upper_bc[par_dir] :
                                                    &gks->lower_bc[par_dir];
      if (bc->type != GKYL_BC_GK_SPECIES_ABSORB &&
          bc->type != GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING) {
        gkyl_exit(
          "The Boltzmann parallel-loss closure requires absorbing or conducting-sheath ion boundaries."
        );
      }
      // All ion species see the same prescribed wall potential at a given end.
      if (ion > 0) {
        struct gk_species *first = &app->species[boltz->ion_idx[0]];
        const struct gkyl_gyrokinetic_bc *ref = edge ? &first->upper_bc[par_dir] :
                                                       &first->lower_bc[par_dir];
        if (bc->aux_profile != ref->aux_profile ||
            (bc->aux_profile && bc->aux_ctx != ref->aux_ctx)) {
          gkyl_exit("Boltzmann ions must share the prescribed wall profile and context at each end."
          );
        }
      }
      if (!boltz->is_boundary_rank[edge]) {
        continue;
      }

      // Reuse the first ion's conducting-wall profile when available.
      if (ion == 0) {
        if (bc->type == GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING) {
          boltz->wall[edge] = edge ? &gks->phi_wall_up : &gks->phi_wall_lo;
        } else {
          // Absorbing species do not own wall data; the field owns this profile.
          struct gkyl_gyrokinetic_bc wall_bc = *bc;
          wall_bc.type = GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING;
          boltz->wall[edge] = &boltz->absorbing_wall[edge];
          gk_species_phi_wall_init(app, &wall_bc, boltz->wall[edge]);
        }
      }

      // Each end_vals owner computes parallel particle losses from its phase-space skin cells.
      struct gkyl_gk_particle_flux_parallel_inp inp = {
        .edge = edge ? GKYL_UPPER_EDGE : GKYL_LOWER_EDGE,
        .cbasis = &app->basis,
        .pbasis = &gks->basis,
        .phase_grid = &gks->grid,
        .conf_skin = &boltz->skin[edge],
        .phase_skin = edge ? &gks->local_upper_skin[par_dir] : &gks->local_lower_skin[par_dir],
        .charge = gks->info.charge,
        .mass = gks->info.mass,
        .gk_geom = app->gk_geom,
        .dg_geom = app->dg_geom,
        .gk_dg_geom = app->gk_dg_geom,
        .vel_map = gks->vel_map,
        .use_gpu = app->use_gpu,
      };
      boltz->m0flux_op[ion][edge] = gkyl_gk_particle_flux_parallel_inew(&inp);
    }
  }
  if (!boltz->num_ions) {
    gkyl_exit("The Boltzmann field requires a positive kinetic species.");
  }
}

static void
gk_field_accumulate_rho_c_boltzmann(
  gkyl_gyrokinetic_app *app, struct gk_field *field, struct gk_species *gks
)
{
  if (gks->info.charge > 0.0) {
    struct gkyl_array *density = field->boltzmann->m0_buff;
    gkyl_dg_div_op_range(
      gks->m0.mem_geo, &app->basis, 0, density, 0, gks->m0.marr, 0, app->gk_geom->geo_int.jacobgeo,
      &app->local
    );
    gkyl_array_accumulate_range(
      field->rho_c, -gks->info.charge / field->info.electron_charge, density, &app->local
    );
  }
}

static void
gk_field_rhs_boltzmann(
  gkyl_gyrokinetic_app *app, struct gk_field *field, double tcurr, const struct gkyl_array *fin[]
)
{
  struct gk_field_boltzmann *boltz = field->boltzmann;
  if (!boltz->num_ions) {
    gkyl_array_clear(field->phi_smooth, 0.0);
    return;
  }

  int par_dir = app->cdim - 1;
  for (int edge = 0; edge < 2; ++edge) {
    if (boltz->is_boundary_rank[edge]) {
      gkyl_array_clear(boltz->m0flux_tot[edge], 0.0);

      // Sum instantaneous outward ion fluxes, weighted by charge number.
      for (int ion = 0; ion < boltz->num_ions; ++ion) {
        int i = boltz->ion_idx[ion];
        gkyl_gk_particle_flux_parallel_advance(
          boltz->m0flux_op[ion][edge], fin[i], boltz->m0flux_buff
        );
        gkyl_array_accumulate(
          boltz->m0flux_tot[edge], -app->species[i].info.charge / field->info.electron_charge,
          boltz->m0flux_buff
        );
      }

      struct gk_species_wall_potential *wall = boltz->wall[edge];
      gk_species_phi_wall_advance(app, wall, tcurr);

      // Form end_vals data from the ion losses, density, magnetic incidence, and wall bias.
      gkyl_boltz_pot_0curr_boundary(
        boltz->pot_op, edge ? GKYL_UPPER_EDGE : GKYL_LOWER_EDGE, &boltz->skin[edge],
        boltz->m0flux_tot[edge], gkyl_gk_particle_flux_parallel_area(boltz->m0flux_op[0][edge]),
        field->rho_c, wall->phi, boltz->end_vals[edge]
      );
    }

    // Make each end_vals's data available to every rank.
    gkyl_comm_array_bcast(
      boltz->comm_par, boltz->end_vals[edge], boltz->end_vals[edge], boltz->boundary_rank[edge]
    );
  }

  // Balance the combined current at both ends to determine the Boltzmann normalization.
  gkyl_boltz_pot_0curr_reference(
    boltz->pot_op, &app->local_lower_ghost[par_dir], boltz->end_vals[0], boltz->end_vals[1],
    boltz->nphi_ref
  );

  // Recover the potential from the local density, then enforce parallel continuity.
  gkyl_boltz_pot_0curr_phi_calc(
    boltz->pot_op, &app->local, &app->local_ext, field->rho_c, boltz->nphi_ref, field->phi_smooth
  );

  gk_field_fem_projection_par(app, field, field->phi_smooth, field->phi_smooth);
}

static void
gk_field_boltzmann_release(const gkyl_gyrokinetic_app *app, struct gk_field_boltzmann *boltz)
{
  for (int ion = 0; ion < boltz->num_ions; ++ion) {
    for (int edge = 0; edge < 2; ++edge) {
      if (boltz->m0flux_op[ion][edge]) {
        gkyl_gk_particle_flux_parallel_release(boltz->m0flux_op[ion][edge]);
      }
    }
  }
  for (int edge = 0; edge < 2; ++edge) {
    if (boltz->wall[edge] == &boltz->absorbing_wall[edge]) {
      gk_species_phi_wall_release(app, boltz->wall[edge]);
    }
    if (boltz->m0flux_tot[edge]) {
      gkyl_array_release(boltz->m0flux_tot[edge]);
    }
    if (boltz->end_vals[edge]) {
      gkyl_array_release(boltz->end_vals[edge]);
    }
  }
  if (boltz->m0flux_buff) {
    gkyl_array_release(boltz->m0flux_buff);
  }
  gkyl_array_release(boltz->m0_buff);
  gkyl_boltz_pot_0curr_release(boltz->pot_op);
  gkyl_array_release(boltz->nphi_ref);
  gkyl_free(boltz);
}

static void
gk_field_fem_release_boltzmann(const gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  gk_field_boltzmann_release(app, gkf->boltzmann);
  gkyl_array_release(gkf->rho_c);
  gkyl_array_release(gkf->rho_c_global_dg);
  gkyl_array_release(gkf->phi_fem);
  gkyl_array_release(gkf->phi_smooth);

  if (app->use_gpu) {
    gkyl_array_release(gkf->phi_host);
  }

  gkyl_fem_parproj_release(gkf->fem_parproj);
  gkyl_array_integrate_release(gkf->calc_em_energy);
}

void
gk_field_fem_new_boltzmann(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  gkf->boltzmann = gkyl_calloc(1, sizeof(struct gk_field_boltzmann));
  gkf->boltzmann->m0_buff = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  // Create global subrange we'll copy the field solver solution from (into local).
  gkyl_sub_range_intersect(&gkf->global_sub_range, &app->global, &app->local);

  // Allocate arrays for charge density.
  gkf->rho_c = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkf->rho_c_global_dg = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);

  // Allocate arrays for electrostatic potential.
  gkf->phi_fem = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  gkf->phi_smooth = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  // Allocate phi_host for I/O.
  gkf->phi_host = gkf->phi_smooth;
  if (app->use_gpu) {
    gkf->phi_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  }

  double polarization_weight = 1.0;
  gkf->rhs_phi_func = gk_field_rhs_boltzmann;
  gkf->accumulate_rhoc_func = gk_field_accumulate_rho_c_boltzmann;

  gkf->boltzmann->pot_op = gkyl_boltz_pot_0curr_new(
    &app->basis, gkf->info.electron_mass, gkf->info.electron_charge, gkf->info.electron_temp,
    app->use_gpu
  );

  gkf->boltzmann->nphi_ref = mkarr(app->use_gpu, 2 * app->basis.num_basis, app->local_ext.volume);

  enum gkyl_fem_parproj_bc_type fem_parproj_bc = GKYL_FEM_PARPROJ_NONE;
  for (int d = 0; d < app->num_periodic_dir; ++d) {
    if (app->periodic_dirs[d] == app->cdim - 1) {
      fem_parproj_bc = GKYL_FEM_PARPROJ_PERIODIC;
    }
  }

  gkf->fem_parproj = gkyl_fem_parproj_new(
    &app->global, &app->grid, &app->basis, fem_parproj_bc, 0, 0, 0, app->use_gpu
  );

  if (app->cdim == 1) {
    gkf->es_energy_fac_1d = polarization_weight;
    gkf->calc_em_energy = gkyl_array_integrate_new(
      &app->grid, &app->basis, 1, GKYL_ARRAY_INTEGRATE_OP_SQ, app->use_gpu
    );
  } else {
    gkf->calc_em_energy = gkyl_array_integrate_new(
      &app->grid, &app->basis, 1, GKYL_ARRAY_INTEGRATE_OP_EPS_GRADPERP_SQ, app->use_gpu
    );
  }

  // Create operator needed for FLR effects.
  gkf->use_flr = false;
  gkf->invert_flr = gk_field_invert_flr_none;
  for (int i = 0; i < app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    if (gks->info.flr.type) {
      gkf->use_flr = gkf->use_flr || gks->info.flr.type;
    }
  }

  gkf->release_func = gk_field_fem_release_boltzmann;

  struct gk_field_boltzmann *boltz = gkf->boltzmann;
  int par_dir = app->cdim - 1;

  // Set up the sheath closure only if a field solve is needed, initially or during evolution.
  bool needs_init_solve = gkf->calc_init_field && !gkf->info.init_field_profile &&
                          gkf->info.init_from_file.type == 0;
  if (!gkf->update_field && !needs_init_solve) {
    return;
  }
  for (int d = 0; d < app->num_periodic_dir; ++d) {
    if (app->periodic_dirs[d] == par_dir) {
      if (gkf->update_field) {
        gkyl_exit("An evolving Boltzmann field requires non-periodic boundaries.");
      }
      return; // Preserve the zero field in static periodic collision tests.
    }
  }

  // Borrow the app's parallel line communicator and identify its boundary ranks.
  boltz->comm_par = app->comm_line[par_dir];
  const struct gkyl_rect_decomp *decomp = app->decomp_line[par_dir];
  for (int r = 0; r < decomp->ndecomp; ++r) {
    if (decomp->ranges[r].lower[par_dir] == app->global.lower[par_dir]) {
      boltz->boundary_rank[0] = r;
    }
    if (decomp->ranges[r].upper[par_dir] == app->global.upper[par_dir]) {
      boltz->boundary_rank[1] = r;
    }
  }

  gkyl_comm_get_rank(boltz->comm_par, &boltz->rank_par);

  // Allocate surface data on every rank for broadcasting values from both ends along the group.
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  for (int edge = 0; edge < 2; ++edge) {
    gkyl_copy_int_arr(app->cdim, app->local.lower, lower);
    gkyl_copy_int_arr(app->cdim, app->local.upper, upper);
    lower[par_dir] = upper[par_dir] = edge ? app->local.upper[par_dir] : app->local.lower[par_dir];
    gkyl_sub_range_init(&boltz->skin[edge], &app->local_ext, lower, upper);

    boltz->is_boundary_rank[edge] = boltz->rank_par == boltz->boundary_rank[edge];

    int nq = 1 << par_dir;
    boltz->m0flux_tot[edge] = mkarr(app->use_gpu, nq, boltz->skin[edge].volume);
    boltz->end_vals[edge] = mkarr(app->use_gpu, 4 * nq, boltz->skin[edge].volume);
  }

  boltz->m0flux_buff = mkarr(app->use_gpu, 1 << par_dir, boltz->skin[0].volume);
}
