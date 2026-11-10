#include <gkyl_gk_species_priv.h>
#include <mpack.h>

// Validate supplied metadata for either species model; missing metadata is allowed.
static bool
gk_species_read_metadata(const struct gk_species *gks, const struct gkyl_array_header_info *hdr)
{
  if (hdr->meta_size == 0) {
    return true;
  }
  mpack_tree_t tree;
  mpack_tree_init_data(&tree, hdr->meta, hdr->meta_size);
  mpack_tree_parse(&tree);
  mpack_node_t root = mpack_tree_root(&tree);

  bool matches = mpack_node_type(root) == mpack_type_map;
  if (matches && mpack_node_map_contains_cstr(root, "basis_type")) {
    mpack_node_t node = mpack_node_map_cstr(root, "basis_type");
    matches = mpack_node_type(node) == mpack_type_str &&
              mpack_node_strlen(node) == strlen(gks->basis.id) &&
              memcmp(mpack_node_str(node), gks->basis.id, strlen(gks->basis.id)) == 0;
  }

  const char *int_keys[] = {"vdim", "poly_order"};
  unsigned int int_values[] = {gks->info.vdim, gks->basis.poly_order};
  for (int i = 0; matches && i < 2; ++i) {
    if (mpack_node_map_contains_cstr(root, int_keys[i])) {
      mpack_node_t node = mpack_node_map_cstr(root, int_keys[i]);
      matches = mpack_node_u64(node) == int_values[i];
    }
  }

  mpack_error_t error = mpack_tree_destroy(&tree);

  return matches && error == mpack_ok;
}

enum gkyl_array_rio_status
gk_species_read(gkyl_gyrokinetic_app *app, struct gk_species *gks, const char *fname)
{
  struct gkyl_rect_grid grid;
  struct gkyl_array_header_info hdr = {0};
  enum gkyl_array_rio_status status = gkyl_grid_sub_array_header_read(&grid, &hdr, fname);
  if (status != GKYL_ARRAY_RIO_SUCCESS) {
    return status;
  }

  bool matches = grid.ndim == gks->grid.ndim && hdr.etype == GKYL_DOUBLE &&
                 hdr.esznc == gks->f->esznc;
  for (int d = 0; matches && d < grid.ndim; ++d) {
    matches = grid.cells[d] == gks->grid.cells[d] &&
              gkyl_compare(grid.lower[d], gks->grid.lower[d], 1e-12) &&
              gkyl_compare(grid.upper[d], gks->grid.upper[d], 1e-12);
  }
  matches = matches && gk_species_read_metadata(gks, &hdr);
  gkyl_array_header_info_release(&hdr);
  if (!matches) {
    return GKYL_ARRAY_RIO_DATA_MISMATCH;
  }

  status = gkyl_comm_array_read(gks->comm, &gks->grid, &gks->local, gks->f_host, fname);
  if (status == GKYL_ARRAY_RIO_SUCCESS && app->use_gpu) {
    gkyl_array_copy(gks->f, gks->f_host);
  }

  return status;
}

void
gk_species_file_import_init(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gkyl_gyrokinetic_ic_import inp
)
{
  if (gks->is_fluid) {
    gk_species_fluid_file_import_init(app, gks, inp);
  } else {
    gk_species_kinetic_file_import_init(app, gks, inp);
  }
}

void
gk_species_gyroaverage_disabled(
  gkyl_gyrokinetic_app *app, struct gk_species *species, struct gkyl_array *field_in,
  struct gkyl_array *field_gyroavg
)
{
  struct timespec wst = gkyl_wall_clock();
  // Don't perform gyroaveraging.
  app->stat.species_gyroavg_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_gyroaverage_enabled(
  gkyl_gyrokinetic_app *app, struct gk_species *species, struct gkyl_array *field_in,
  struct gkyl_array *field_gyroavg
)
{
  struct timespec wst = gkyl_wall_clock();
  // Gyroaverage input field.
  gkyl_deflated_fem_poisson_advance(species->flr_op, field_in, field_in, field_gyroavg);
  app->stat.species_gyroavg_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_rhs_static(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms
)
{
  // Do nothing.
}

double
gk_species_get_cfl_static(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms
)
{
  double omega_cfl = 1 / DBL_MAX;
  return app->cfl / omega_cfl;
}

double
gk_species_rhs_implicit_static(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt
)
{
  double omega_cfl = 1 / DBL_MAX;
  return app->cfl / omega_cfl;
}

static void
gk_species_phi_wall_advance_disabled(
  gkyl_gyrokinetic_app *app, const struct gk_species_wall_potential *wall, double tm
)
{
}

static void
gk_species_phi_wall_advance_enabled(
  gkyl_gyrokinetic_app *app, const struct gk_species_wall_potential *wall, double tm
)
{
  gkyl_eval_on_nodes_advance(wall->projector, tm, &app->local_ext, wall->phi_host);
  if (app->use_gpu) {
    gkyl_array_copy(wall->phi, wall->phi_host);
  }
}

void
gk_species_phi_wall_init(
  gkyl_gyrokinetic_app *app, const struct gkyl_gyrokinetic_bc *bc,
  struct gk_species_wall_potential *wall
)
{
  *wall = (struct gk_species_wall_potential){.advance_func = gk_species_phi_wall_advance_disabled};
  if (bc->type != GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING &&
      bc->type != GKYL_BC_GK_SPECIES_SHEATH_SURROGATE) {
    return;
  }

  wall->phi = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkyl_array_clear(wall->phi, 0.0);

  wall->phi_host = wall->phi;
  if (bc->aux_profile) {
    if (app->use_gpu) {
      wall->phi_host = mkarr(false, wall->phi->ncomp, wall->phi->size);
    }

    wall->projector =
      gkyl_eval_on_nodes_new(&app->grid, &app->basis, 1, bc->aux_profile, bc->aux_ctx);
    wall->advance_func = gk_species_phi_wall_advance_enabled;
    gk_species_phi_wall_advance(app, wall, 0.0);
  }
}

void
gk_species_phi_wall_advance(
  gkyl_gyrokinetic_app *app, const struct gk_species_wall_potential *wall, double tm
)
{
  wall->advance_func(app, wall, tm);
}

void
gk_species_phi_wall_release(
  const gkyl_gyrokinetic_app *app, const struct gk_species_wall_potential *wall
)
{
  if (!wall->phi) {
    return;
  }

  gkyl_array_release(wall->phi);
  if (wall->projector) {
    gkyl_eval_on_nodes_release(wall->projector);
    if (app->use_gpu) {
      gkyl_array_release(wall->phi_host);
    }
  }
}

void
gk_species_apply_bc_static(
  gkyl_gyrokinetic_app *app, struct gk_species *species, double tm, struct gkyl_array *f
)
{
  // do nothing
}

void
gk_species_step_f_dynamic(struct gkyl_array *out, double dt, const struct gkyl_array *inp)
{
  gkyl_array_accumulate(gkyl_array_scale(out, dt), 1.0, inp);
}

void
gk_species_step_f_static(struct gkyl_array *out, double dt, const struct gkyl_array *inp)
{
  // do nothing
}

void
gk_species_combine_dynamic(
  struct gkyl_array *out, double c1, const struct gkyl_array *arr1, double c2,
  const struct gkyl_array *arr2, const struct gkyl_range *rng
)
{
  gkyl_array_accumulate_range(gkyl_array_set_range(out, c1, arr1, rng), c2, arr2, rng);
}

void
gk_species_combine_static(
  struct gkyl_array *out, double c1, const struct gkyl_array *arr1, double c2,
  const struct gkyl_array *arr2, const struct gkyl_range *rng
)
{
  // do nothing
}

void
gk_species_copy_range_dynamic(
  struct gkyl_array *out, const struct gkyl_array *inp, const struct gkyl_range *range
)
{
  gkyl_array_copy_range(out, inp, range);
}

void
gk_species_copy_range_static(
  struct gkyl_array *out, const struct gkyl_array *inp, const struct gkyl_range *range
)
{
  // do nothing
}

void
gk_species_write_vcutsq(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gkyl_bc_sheath_gyrokinetic *bc_sheath,
  const char *dir_str, const char *edge_str, double tm, int frame
)
{
  struct timespec wst = gkyl_wall_clock();

  // Package metadata (time/frame + DG basis info).
  gkyl_msgpack_map_elem_set_double(app->io_meta_dg_len, app->io_meta_dg, "time", tm);
  gkyl_msgpack_map_elem_set_uint(app->io_meta_dg_len, app->io_meta_dg, "frame", frame);
  struct gkyl_msgpack_map_elem mpe_vcutsq[] = {
    {.key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = app->basis.poly_order},
    {.key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = "serendipity"},
    {
      .key = "Description",
      .elem_type = GKYL_MP_STRING,
      .cval = "Square of the velocity cutoff used in the gyrokinetic sheath BC.",
    }
  };
  int mpe_vcutsq_len = sizeof(mpe_vcutsq) / sizeof(mpe_vcutsq[0]);
  int io_meta_len[] = {app->io_meta_dg_len, mpe_vcutsq_len};
  const struct gkyl_msgpack_map_elem *io_meta[] = {app->io_meta_dg, mpe_vcutsq};
  struct gkyl_msgpack_data *mt =
    gkyl_msgpack_create_union(sizeof(io_meta_len) / sizeof(int), io_meta_len, io_meta);

  const char *fmt = "%s-%s_bc_%s%s_vcutsq_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, dir_str, edge_str, frame);
  char fileNm[sz + 1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, dir_str, edge_str, frame);

  gkyl_bc_sheath_gyrokinetic_write_vcutsq(bc_sheath, mt, fileNm);

  gkyl_msgpack_data_release(mt);
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag_io += 1;
}

void
gk_species_write_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  gks->write_cfl_func(app, gks, tm, frame);

  // Metadata from app and geo object.
  gkyl_msgpack_map_elem_set_double(gks->io_meta_phase_len, gks->io_meta_phase, "time", tm);
  gkyl_msgpack_map_elem_set_uint(gks->io_meta_phase_len, gks->io_meta_phase, "frame", frame);
  struct gkyl_msgpack_map_elem io_meta_f[] = {
    {
      .key = "Description",
      .elem_type = GKYL_MP_STRING,
      .cval = gks->is_fluid ? "Fluid moments: J*rho, J*Mpar, J*E, J*pperp/B." :
                              "Distribution function times Jacobians.",
    },
    {.key = "is_fluid", .elem_type = GKYL_MP_BOOL, .bval = gks->is_fluid}
  };
  int io_meta_f_len = sizeof(io_meta_f) / sizeof(io_meta_f[0]);
  int io_meta_len[] = {gks->io_meta_phase_len, app->gk_geom->io_meta_basic_len, io_meta_f_len};
  const struct gkyl_msgpack_map_elem *io_meta[] = {
    gks->io_meta_phase, app->gk_geom->io_meta_basic, io_meta_f
  };
  struct gkyl_msgpack_data *mt =
    gkyl_msgpack_create_union(sizeof(io_meta_len) / sizeof(int), io_meta_len, io_meta);

  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz + 1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);

  // copy data from device to host before writing it out
  if (app->use_gpu) {
    gkyl_array_copy(gks->f_host, gks->f);
  }
  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->f_host, fileNm);

  gkyl_msgpack_data_release(mt);

  app->stat.species_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_io += 1;

  // Write out the sheath BC velocity cutoff (vcutsq).
  // Every rank has local skin ranges, but only physical boundary ranks may
  // write these serial files. Decomposition is only along the parallel direction.
  int par_dir = app->cdim - 1; // Sheath BC acts in the parallel direction.
  if (gks->local.lower[par_dir] == gks->global.lower[par_dir] &&
      (gks->lower_bc[par_dir].type == GKYL_BC_GK_SPECIES_SHEATH_SURROGATE ||
       gks->lower_bc[par_dir].type == GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING)) {
    gk_species_write_vcutsq(app, gks, gks->bc_sheath_lo, "z", "lower", tm, frame);
  }
  if (gks->local.upper[par_dir] == gks->global.upper[par_dir] &&
      (gks->upper_bc[par_dir].type == GKYL_BC_GK_SPECIES_SHEATH_SURROGATE ||
       gks->upper_bc[par_dir].type == GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING)) {
    gk_species_write_vcutsq(app, gks, gks->bc_sheath_up, "z", "upper", tm, frame);
  }
}

void
gk_species_write_static(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame)
{
  // do nothing
}

void
gk_species_write_cfl_enabled(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame)
{
  // The 0th frame is often wrong for this diagnostic because cflrate is computed AFTER
  // the first time step
  struct timespec wst = gkyl_wall_clock();
  // DG metadata for cflrate.
  struct gkyl_msgpack_map_elem mpe_cfl[] = {
    {.key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0},
    {.key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = "serendipity"},
    {.key = "Description", .elem_type = GKYL_MP_STRING, .cval = "CFL frequency."},
    {.key = "value_form", .elem_type = GKYL_MP_STRING, .cval = "nodal"}
  };
  int mpe_cfl_len = sizeof(mpe_cfl) / sizeof(mpe_cfl[0]);
  // Update app basic metada with time/frame.
  gkyl_msgpack_map_elem_set_double(gks->io_meta_phase_len, gks->io_meta_phase, "time", tm);
  gkyl_msgpack_map_elem_set_uint(gks->io_meta_phase_len, gks->io_meta_phase, "frame", frame);
  // Unionize metadata.
  int io_meta_len[] = {gks->io_meta_phase_len, mpe_cfl_len, app->gk_geom->io_meta_basic_len};
  const struct gkyl_msgpack_map_elem *io_meta[] = {
    gks->io_meta_phase, mpe_cfl, app->gk_geom->io_meta_basic
  };
  struct gkyl_msgpack_data *mt =
    gkyl_msgpack_create_union(sizeof(io_meta_len) / sizeof(int), io_meta_len, io_meta);

  const char *fmt = "%s-%s_cflrate_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz + 1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);
  gkyl_array_copy(gks->cflrate_ho, gks->cflrate);
  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->cflrate_ho, fileNm);

  gkyl_msgpack_data_release(mt);

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_io += 1;
}

void
gk_species_write_cfl_disabled(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame
)
{
  // do nothing
}

void
gk_species_write_mom_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame)
{
  // Update time/frame metadata (shared fields, set once before the loop).
  gkyl_msgpack_map_elem_set_double(gks->io_meta_conf_len, gks->io_meta_conf, "time", tm);
  gkyl_msgpack_map_elem_set_uint(gks->io_meta_conf_len, gks->io_meta_conf, "frame", frame);

  for (int m = 0; m < gks->info.num_diag_moments; ++m) {
    struct timespec wtm = gkyl_wall_clock();
    gk_species_moment_calc(app, &gks->moms[m], gks->local, app->local, gks->f);
    app->stat.n_mom += 1;

    // Rescale moment by inverse of Jacobian if necessary.
    gk_species_moment_diag_jacobgeo_div(app, &gks->moms[m], gks->moms[m].marr, gks->moms[m].marr);
    app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wtm);

    struct timespec wst = gkyl_wall_clock();
    if (app->use_gpu) {
      gkyl_array_copy(gks->moms[m].marr_host, gks->moms[m].marr);
    }

    struct gkyl_msgpack_map_elem io_meta_mom[] = {{
      .key = "Description",
      .elem_type = GKYL_MP_STRING,
      .cval = (char *)gkyl_distribution_moments_descriptions[gks->info.diag_moments[m]],
    }};
    int io_meta_len[] = {gks->io_meta_conf_len, app->gk_geom->io_meta_basic_len, 1};
    const struct gkyl_msgpack_map_elem *io_meta[] = {
      gks->io_meta_conf, app->gk_geom->io_meta_basic, io_meta_mom
    };
    struct gkyl_msgpack_data *mt =
      gkyl_msgpack_create_union(sizeof(io_meta_len) / sizeof(int), io_meta_len, io_meta);

    const char *fmt = "%s-%s_%s_%d.gkyl";
    int sz = gkyl_calc_strlen(
      fmt, app->name, gks->info.name, gkyl_distribution_moments_strs[gks->info.diag_moments[m]],
      frame
    );
    char fileNm[sz + 1]; // ensures no buffer overflow
    snprintf(
      fileNm, sizeof fileNm, fmt, app->name, gks->info.name,
      gkyl_distribution_moments_strs[gks->info.diag_moments[m]], frame
    );

    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gks->moms[m].marr_host, fileNm);
    gkyl_msgpack_data_release(mt);
    app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
    app->stat.n_diag_io += 1;
  }

  app->stat.n_diag += 1;
}

void
gk_species_write_mom_static(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame)
{
  // do nothing
}

void
gk_species_calc_int_mom_dt_enabled(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, double dt, struct gkyl_array *fdot_int_mom
)
{
  struct timespec wst = gkyl_wall_clock();
  // Compute moment of f_new to compute moment of df/dt.
  // Need to do it after the fields are updated.
  gk_species_moment_calc(app, &gks->integ_moms, gks->local, app->local, gks->f);
  gkyl_array_set(fdot_int_mom, 1.0 / dt, gks->integ_moms.marr);
  app->stat.fdot_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_calc_int_mom_dt_disabled(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, double dt, struct gkyl_array *fdot_int_mom
)
{
}

void
gk_species_calc_int_mom_dt(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, double dt, struct gkyl_array *fdot_int_mom
)
{
  gks->calc_int_mom_dt_func(app, gks, dt, fdot_int_mom);
}

void
gk_species_calc_integrated_mom_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm)
{
  struct timespec wst = gkyl_wall_clock();

  int num_mom = gks->integ_moms.num_mom;
  double avals_global[num_mom];

  gk_species_moment_calc(app, &gks->integ_moms, gks->local, app->local, gks->f);
  app->stat.n_mom += 1;

  // Reduce (sum) over whole domain, append to diagnostics.
  gkyl_array_reduce_range(gks->red_integ_diag, gks->integ_moms.marr, GKYL_SUM, &app->local);
  gkyl_comm_allreduce(
    app->comm, GKYL_DOUBLE, GKYL_SUM, num_mom, gks->red_integ_diag, gks->red_integ_diag_global
  );
  if (app->use_gpu) {
    gkyl_cu_memcpy(
      avals_global, gks->red_integ_diag_global, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H
    );
  } else {
    memcpy(avals_global, gks->red_integ_diag_global, sizeof(double[num_mom]));
  }
  gkyl_dynvec_append(gks->integ_diag, tm, avals_global);

  if (gks->info.time_rate_diagnostics) {
    // Reduce (sum) over whole domain, append to diagnostics.
    gkyl_array_accumulate(gks->fdot_mom_new, -1.0, gks->fdot_mom_old);
    gkyl_array_reduce_range(gks->red_integ_diag, gks->fdot_mom_new, GKYL_SUM, &app->local);
    gkyl_comm_allreduce(
      app->comm, GKYL_DOUBLE, GKYL_SUM, num_mom, gks->red_integ_diag, gks->red_integ_diag_global
    );
    if (app->use_gpu) {
      gkyl_cu_memcpy(
        avals_global, gks->red_integ_diag_global, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H
      );
    } else {
      memcpy(avals_global, gks->red_integ_diag_global, sizeof(double[num_mom]));
    }
    gkyl_dynvec_append(gks->fdot_integ_diag, tm, avals_global);
  }

  app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag += 1;
}

void
gk_species_calc_integrated_mom_static(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm)
{
  // do nothing
}

void
gk_species_write_integrated_mom_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  struct timespec wst = gkyl_wall_clock();

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    // Write integrated diagnostic moments.
    const char *fmt = "%s-%s_%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "integrated_moms");
    char fileNm[sz + 1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "integrated_moms");

    if (gks->is_first_integ_write_call) {
      struct gkyl_msgpack_map_elem io_meta_phi[] = {
        {.key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Volume integrated moments."}
      };
      int io_meta_len[] = {gks->io_meta_basic_len, app->gk_geom->io_meta_basic_len, 1};
      const struct gkyl_msgpack_map_elem *io_meta[] = {
        gks->io_meta_basic, app->gk_geom->io_meta_basic, io_meta_phi
      };
      struct gkyl_msgpack_data *mt =
        gkyl_msgpack_create_union(sizeof(io_meta_len) / sizeof(int), io_meta_len, io_meta);

      gkyl_dynvec_write_wmeta(gks->integ_diag, fileNm, mt);
      gks->is_first_integ_write_call = false;
      gkyl_msgpack_data_release(mt);
    } else {
      gkyl_dynvec_awrite(gks->integ_diag, fileNm);
    }
  }
  gkyl_dynvec_clear(gks->integ_diag);
  app->stat.n_diag_io += 1;

  if (gks->info.time_rate_diagnostics) {
    if (rank == 0) {
      // Write integrated diagnostic moments.
      const char *fmt = "%s-%s_fdot_%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "integrated_moms");
      char fileNm[sz + 1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "integrated_moms");

      if (gks->is_first_fdot_integ_write_call) {
        struct gkyl_msgpack_map_elem io_meta_phi[] = {{
          .key = "Description",
          .elem_type = GKYL_MP_STRING,
          .cval = "Volume integrated moments of time rate of change.",
        }};
        int io_meta_len[] = {gks->io_meta_basic_len, app->gk_geom->io_meta_basic_len, 1};
        const struct gkyl_msgpack_map_elem *io_meta[] = {
          gks->io_meta_basic, app->gk_geom->io_meta_basic, io_meta_phi
        };
        struct gkyl_msgpack_data *mt =
          gkyl_msgpack_create_union(sizeof(io_meta_len) / sizeof(int), io_meta_len, io_meta);

        gkyl_dynvec_write_wmeta(gks->fdot_integ_diag, fileNm, mt);
        gks->is_first_fdot_integ_write_call = false;
        gkyl_msgpack_data_release(mt);
      } else {
        gkyl_dynvec_awrite(gks->fdot_integ_diag, fileNm);
      }
    }
    gkyl_dynvec_clear(gks->fdot_integ_diag);
    app->stat.n_diag_io += 1;
  }

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_write_integrated_mom_static(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // do nothing
}

void
gk_species_calc_L2norm_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_array_integrate_advance(
    gks->integ_wfsq_op, gks->f, (2.0 * M_PI) / gks->info.mass, app->jacobtot_inv_weak, &gks->local,
    &app->local, gks->L2norm_local
  );
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, gks->L2norm_local, gks->L2norm_global);

  double L2norm_global[] = {0.0};
  if (app->use_gpu) {
    gkyl_cu_memcpy(L2norm_global, gks->L2norm_global, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    memcpy(L2norm_global, gks->L2norm_global, sizeof(double));
  }
  gkyl_dynvec_append(gks->L2norm, tm, L2norm_global);

  app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag += 1;
}

void
gk_species_calc_L2norm_static(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm)
{
  // do nothing
}

void
gk_species_write_L2norm_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  struct timespec wst = gkyl_wall_clock();

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    // Write the L2 norm.
    const char *fmt = "%s-%s_%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "L2norm");
    char fileNm[sz + 1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "L2norm");

    if (gks->is_first_L2norm_write_call) {
      struct gkyl_msgpack_map_elem io_meta_phi[] = {
        {.key = "Description", .elem_type = GKYL_MP_STRING, .cval = "L2 norm of the distribution."}
      };
      int io_meta_len[] = {gks->io_meta_basic_len, app->gk_geom->io_meta_basic_len, 1};
      const struct gkyl_msgpack_map_elem *io_meta[] = {
        gks->io_meta_basic, app->gk_geom->io_meta_basic, io_meta_phi
      };
      struct gkyl_msgpack_data *mt =
        gkyl_msgpack_create_union(sizeof(io_meta_len) / sizeof(int), io_meta_len, io_meta);

      gkyl_dynvec_write_wmeta(gks->L2norm, fileNm, mt);
      gks->is_first_L2norm_write_call = false;
      gkyl_msgpack_data_release(mt);
    } else {
      gkyl_dynvec_awrite(gks->L2norm, fileNm);
    }
  }
  gkyl_dynvec_clear(gks->L2norm);

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag_io += 1;
}

void
gk_species_write_L2norm_static(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // do nothing
}

void
gk_species_release_static(const gkyl_gyrokinetic_app *app, const struct gk_species *gks)
{
}

void
gk_species_init_static(
  struct gkyl_gk *gk_app_inp, struct gkyl_gyrokinetic_app *app, struct gk_species *gks
)
{
  // Allocate distribution function arrays.
  gks->f1 = gks->f;
  gks->fnew = gks->f;

  // Set function pointers.
  gks->rhs_func = gk_species_rhs_static;
  gks->rhs_em_complete_func = gk_species_rhs_static;
  gks->bflux_update = gk_species_rhs_static;
  gks->get_cfl = gk_species_get_cfl_static;
  gks->rhs_implicit_func = gk_species_rhs_implicit_static;
  gks->bc_func = gk_species_apply_bc_static;
  gks->release_func = gk_species_release_static;
  gks->step_f_func = gk_species_step_f_static;
  gks->combine_func = gk_species_combine_static;
  gks->copy_func = gk_species_copy_range_static;
  gks->write_func = gk_species_write_static;
  gks->write_cfl_func = gk_species_write_cfl_disabled;
  gks->write_mom_func = gk_species_write_mom_static;
  gks->calc_integrated_mom_func = gk_species_calc_integrated_mom_static;
  gks->write_integrated_mom_func = gk_species_write_integrated_mom_static;
  gks->calc_L2norm_func = gk_species_calc_L2norm_static;
  gks->write_L2norm_func = gk_species_write_L2norm_static;
  gks->calc_int_mom_dt_func = gk_species_calc_int_mom_dt_disabled;
}

void
gk_species_apply_ic(gkyl_gyrokinetic_app *app, struct gk_species *gks, double t0)
{
  if (gks->info.init_from_file.type == 0) {
    gk_species_projection_calc(app, gks, &gks->proj_init, gks->f, t0);
  }

  // We are pre-computing source for now as it is time-independent.
  gk_species_source_calc(app, gks, &gks->src, gks->lte.f_lte, t0);
}

void
gk_species_apply_ic_cross(gkyl_gyrokinetic_app *app, struct gk_species *gks_self, double t0)
{
  // IC setup step that depends on the IC of other species.
  if (app->field->init_phi_pol && gks_self->info.scale_with_polarization) {
    // Scale the distribution function so its guiding center density is computed
    // given the polarization density and the guiding center density of other species.

    // Compute the polarization density.
    struct gkyl_array *npol = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

    struct gkyl_gyrokinetic_pol_density *npol_op =
      gkyl_gyrokinetic_pol_density_new(app->basis, app->grid, app->use_gpu);
    gkyl_gyrokinetic_pol_density_advance(
      npol_op, &app->local, app->field->epsilon, app->field->phi_pol, npol
    );
    gkyl_gyrokinetic_pol_density_release(npol_op);

    // Calculate the guiding center density of this species: (npol - q_other*n^G_other)/q_self.
    for (int i = 0; i < app->num_species; ++i) {
      struct gk_species *gks = &app->species[i];
      gk_species_moment_calc(app, &gks->m0, gks->local, app->local, gks->f);
      if (strcmp(gks->info.name, gks_self->info.name)) {
        gkyl_array_accumulate(npol, -gks->info.charge, gks->m0.marr);
      }
    }
    gkyl_array_scale(npol, 1. / gks_self->info.charge);

    // Scale the distribution function so it has this guiding center density.
    struct gkyl_dg_bin_op_mem *div_mem =
      app->use_gpu ? gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->basis.num_basis) :
                     gkyl_dg_bin_op_mem_new(app->local.volume, app->basis.num_basis);
    struct gkyl_array *den_mod = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

    gkyl_dg_div_op_range(
      div_mem, &app->basis, 0, den_mod, 0, npol, 0, gks_self->m0.marr, &app->local
    );
    gkyl_dg_mul_conf_phase_op_range(
      &app->basis, &gks_self->basis, gks_self->f, den_mod, gks_self->f, &app->local_ext,
      &gks_self->local_ext
    );

    gkyl_array_release(den_mod);
    gkyl_dg_bin_op_mem_release(div_mem);
    gkyl_array_release(npol);
  }

  // Store initial density in scaling operator.
  gk_species_scaling_apply_ic_cross(app, gks_self, &gks_self->sca);
}

void
gk_species_rhs(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms
)
{
  species->rhs_func(app, species, fin, rhs, bflux_moms);
}

void
gk_species_rhs_em_complete(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms
)
{
  species->rhs_em_complete_func(app, species, fin, rhs, bflux_moms);
}

void
gk_species_update_bflux(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms
)
{
  species->bflux_update(app, species, fin, rhs, bflux_moms);
}

double
gk_species_get_cfl(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms
)
{
  return species->get_cfl(app, species, fin, rhs, bflux_moms);
}

double
gk_species_rhs_implicit(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt
)
{
  return species->rhs_implicit_func(app, species, fin, rhs, bflux_moms, dt);
}

void
gk_species_step_f(
  struct gk_species *species, struct gkyl_array *out, double a, const struct gkyl_array *inp
)
{
  species->step_f_func(out, a, inp);
}

void
gk_species_combine(
  struct gk_species *species, struct gkyl_array *out, double c1, const struct gkyl_array *arr1,
  double c2, const struct gkyl_array *arr2, const struct gkyl_range *rng
)
{
  species->combine_func(out, c1, arr1, c2, arr2, rng);
}

void
gk_species_copy_range(
  struct gk_species *species, struct gkyl_array *out, const struct gkyl_array *inp,
  const struct gkyl_range *range
)
{
  species->copy_func(out, inp, range);
}

void
gk_species_apply_bc(
  gkyl_gyrokinetic_app *app, struct gk_species *species, double tm, struct gkyl_array *f
)
{
  species->bc_func(app, species, tm, f);
}

void
gk_species_n_iter_corr(gkyl_gyrokinetic_app *app, const struct gk_species *gks, int sidx)
{
  app->stat.num_corr[sidx] = gks->lte.num_corr;
  app->stat.n_iter_corr[sidx] = gks->lte.n_iter;
}

void
gk_species_write(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame)
{
  if (frame == 0) {
    gk_species_write_dynamic(app, gks, tm, frame);
  } else {
    gks->write_func(app, gks, tm, frame);
  }
}

void
gk_species_write_mom(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm, int frame)
{
  if (frame == 0) {
    gk_species_write_mom_dynamic(app, gks, tm, frame);
  } else {
    gks->write_mom_func(app, gks, tm, frame);
  }
}

void
gk_species_calc_integrated_mom(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm)
{
  gks->calc_integrated_mom_func(app, gks, tm);
}

void
gk_species_write_integrated_mom(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  gks->write_integrated_mom_func(app, gks);
}

void
gk_species_calc_L2norm(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm)
{
  gks->calc_L2norm_func(app, gks, tm);
}

void
gk_species_write_L2norm(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  gks->write_L2norm_func(app, gks);
}

void
gk_species_init_io(struct gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // Species properties metadata.
  struct gkyl_msgpack_map_elem io_meta_sprop[] = {
    {.key = "mass", .elem_type = GKYL_MP_DOUBLE, .dval = gks->info.mass},
    {.key = "charge", .elem_type = GKYL_MP_DOUBLE, .dval = gks->info.charge},
    {.key = "vdim", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = gks->info.vdim}
  };

  // Metadata for integrated quantities.
  const struct gkyl_msgpack_map_elem *io_meta_basic_union[] = {app->io_meta_basic, io_meta_sprop};
  int io_meta_basic_union_len[] = {
    app->io_meta_basic_len, sizeof(io_meta_sprop) / sizeof(io_meta_sprop[0])
  };
  gks->io_meta_basic = gkyl_msgpack_map_elem_union(
    sizeof(io_meta_basic_union) / sizeof(io_meta_basic_union[0]), io_meta_basic_union_len,
    io_meta_basic_union, &gks->io_meta_basic_len
  );

  // Metadata for conf-space quantities.
  struct gkyl_msgpack_map_elem io_meta_conf[] = {
    {.key = "value_form", .elem_type = GKYL_MP_STRING, .cval = "modal"},
    {.key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = app->basis.poly_order},
    {.key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = app->basis.id},
    {.key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0},
    {.key = "frame", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0}
  };
  const struct gkyl_msgpack_map_elem *io_meta_conf_union[] = {
    app->io_meta_basic, io_meta_sprop, io_meta_conf
  };
  int io_meta_conf_union_len[] = {
    app->io_meta_basic_len, sizeof(io_meta_sprop) / sizeof(io_meta_sprop[0]),
    sizeof(io_meta_conf) / sizeof(io_meta_conf[0])
  };
  gks->io_meta_conf = gkyl_msgpack_map_elem_union(
    sizeof(io_meta_conf_union) / sizeof(io_meta_conf_union[0]), io_meta_conf_union_len,
    io_meta_conf_union, &gks->io_meta_conf_len
  );

  // Metadata for phase-space quantities.
  struct gkyl_msgpack_map_elem io_meta_phase[] = {
    {.key = "value_form", .elem_type = GKYL_MP_STRING, .cval = "modal"},
    {.key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = gks->basis.poly_order},
    {.key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = gks->basis.id},
    {.key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0},
    {.key = "frame", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0}
  };
  const struct gkyl_msgpack_map_elem *io_meta_phase_union[] = {
    app->io_meta_basic, io_meta_sprop, io_meta_phase
  };
  int io_meta_phase_union_len[] = {
    app->io_meta_basic_len, sizeof(io_meta_sprop) / sizeof(io_meta_sprop[0]),
    sizeof(io_meta_phase) / sizeof(io_meta_phase[0])
  };
  gks->io_meta_phase = gkyl_msgpack_map_elem_union(
    sizeof(io_meta_phase_union) / sizeof(io_meta_phase_union[0]), io_meta_phase_union_len,
    io_meta_phase_union, &gks->io_meta_phase_len
  );
}

void
gk_species_init_integrated_mom(struct gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // Allocate data for integrated moments.
  int num_diag_int_moms = gks->info.num_integrated_diag_moments;
  assert(num_diag_int_moms < 2); // 1 int moment allowed now.
  gk_species_moment_init(
    app, gks, &gks->integ_moms,
    num_diag_int_moms == 0 ? GKYL_F_MOMENT_M0M1M2PARM2PERP : gks->info.integrated_diag_moments[0],
    true
  );

  if (app->use_gpu) {
    gks->red_integ_diag = gkyl_cu_malloc(sizeof(double[gks->integ_moms.num_mom]));
    gks->red_integ_diag_global = gkyl_cu_malloc(sizeof(double[gks->integ_moms.num_mom]));
  } else {
    gks->red_integ_diag = gkyl_malloc(sizeof(double[gks->integ_moms.num_mom]));
    gks->red_integ_diag_global = gkyl_malloc(sizeof(double[gks->integ_moms.num_mom]));
  }
  // Allocate dynamic-vector to store all-reduced integrated moments.
  gks->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, gks->integ_moms.num_mom);
  gks->is_first_integ_write_call = true;

  // Allocate dynamic-vector to store Delta f integrated moments.
  if (gks->info.time_rate_diagnostics) {
    gks->fdot_mom_old =
      mkarr(app->use_gpu, gks->integ_moms.marr->ncomp, gks->integ_moms.marr->size);
    gks->fdot_mom_new =
      mkarr(app->use_gpu, gks->integ_moms.marr->ncomp, gks->integ_moms.marr->size);
    gks->fdot_integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, gks->integ_moms.num_mom);
    gks->is_first_fdot_integ_write_call = true;
  }
}

void
gk_species_release_integrated_mom(
  const struct gkyl_gyrokinetic_app *app, const struct gk_species *gks
)
{
  // Release integrated moment memory.
  gk_species_moment_release(app, &gks->integ_moms);

  // Release integrated diag memory.
  gkyl_dynvec_release(gks->integ_diag);
  if (app->use_gpu) {
    gkyl_cu_free(gks->red_integ_diag);
    gkyl_cu_free(gks->red_integ_diag_global);
  } else {
    gkyl_free(gks->red_integ_diag);
    gkyl_free(gks->red_integ_diag_global);
  }

  if (gks->info.time_rate_diagnostics) {
    // Free df/dt diagnostics memory.
    gkyl_array_release(gks->fdot_mom_old);
    gkyl_array_release(gks->fdot_mom_new);
    gkyl_dynvec_release(gks->fdot_integ_diag);
  }
}

void
gk_species_init(struct gkyl_gk *inp, struct gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  struct gkyl_gyrokinetic_species info = gks->info;
  if (info.is_static && info.scaling.type) {
    fprintf(stderr, "Static species '%s' cannot use a scaling operation.\n", info.name);
    gkyl_exit("Static species must remain unchanged in time.");
  }

  memset(gks, 0, sizeof *gks);

  gks->info = info;
  gks->is_fluid = info.cells[0] == 0;

  if (gks->is_fluid) {
    gk_species_fluid_init(inp, app, gks);
  } else {
    gk_species_kinetic_init(inp, app, gks);
  }
}

void
gk_species_release(const gkyl_gyrokinetic_app *app, const struct gk_species *gks)
{
  gks->release_species_func(app, gks);
}
