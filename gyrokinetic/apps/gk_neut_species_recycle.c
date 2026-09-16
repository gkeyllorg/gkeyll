#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_neut_species_recycle_write_flux_enabled(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_recycle_wall *recyc, double tm, int frame)
{
  // Output boundary flux from ions and neutral ghost cells

  // Package metadata.
  gkyl_msgpack_map_elem_set_double(s->io_meta_conf_len, s->io_meta_conf, "time", tm);
  gkyl_msgpack_map_elem_set_uint(s->io_meta_conf_len, s->io_meta_conf, "frame", frame);
  int io_meta_len[] = {s->io_meta_conf_len, app->gk_geom->io_meta_basic_len, 1};

  const char *vars[] = {"x","y","z"};
  const char *edge[] = {"lower","upper"};

  int dir = recyc->dir;
  int edi = recyc->edge == GKYL_LOWER_EDGE? 0 : 1;

  struct gkyl_range *cskin_r = edi ==0 ? &app->local_lower_skin[recyc->dir] : &app->local_upper_skin[recyc->dir];
  
  for (int i=0; i<recyc->num_species; ++i) {
    struct timespec wst = gkyl_wall_clock();
    // Write out the particle flux of the impacting species.
    struct gk_species *gks = recyc->impact_species[i];

    gk_species_bflux_get_flux(&gks->bflux, recyc->dir, recyc->edge, recyc->phase_flux_gk[i], &recyc->impact_buff_r[i]);
    gkyl_dg_updater_moment_gyrokinetic_advance(recyc->m0op_gk[i], &recyc->impact_normal_r[i],
      &recyc->impact_cbuff_r[i], recyc->phase_flux_gk[i], recyc->m0_flux_gk[i]);
    app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);

    struct timespec wtm = gkyl_wall_clock();
    // Copy to skin to write it out.
    gkyl_array_clear(recyc->diag_out, 0.0);
    gkyl_array_copy_range_to_range(recyc->diag_out, recyc->m0_flux_gk[i], cskin_r, &recyc->impact_cbuff_r[i]);

    if (app->use_gpu)
      gkyl_array_copy(recyc->diag_out_ho, recyc->diag_out);

    struct gkyl_msgpack_map_elem desc0[] = {
      { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Impacting boundary particle flux." }
    };
    const struct gkyl_msgpack_map_elem* io_meta[] = {s->io_meta_conf, app->gk_geom->io_meta_basic, desc0};
    struct gkyl_msgpack_data *mt0 = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

    const char *fmt = "%s-%s_recycling_%s%s_%s_flux_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, vars[dir], edge[edi], gks->info.name, frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, vars[dir], edge[edi], gks->info.name, frame);
    
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt0, recyc->diag_out_ho, fileNm);
    app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wtm);

    gkyl_msgpack_data_release(mt0); 
  }

  struct timespec wst = gkyl_wall_clock();
  // Write out the particle flux of the emitting neutral species.
  gkyl_array_clear(recyc->f_diag, 0.0);
  gkyl_array_copy_range_to_range(recyc->f_diag, recyc->f_emit, recyc->emit_skin_r, &recyc->emit_buff_r);
  // The from-flux Vlasov equation's boundary term reads the stored conf-space
  // surface fluxes, so compute them for this distribution first (cflrate is
  // scratch here; it is cleared at the start of every RHS evaluation).
  gkyl_dg_vlasov_conf_flux_surf_advance(s->collisionless.calc_conf_flux, &app->local,
    &s->local, &s->local_ext, s->conf_poisson_tensor, s->hamil, recyc->f_diag,
    s->cflrate, s->collisionless.conf_flux_surf);
  gkyl_boundary_flux_advance(recyc->f0_flux_slvr, recyc->f_diag, recyc->f_diag);
  gkyl_array_copy_range_to_range(recyc->unit_phase_flux_neut, recyc->f_diag, &recyc->emit_buff_r, recyc->emit_ghost_r);
  gkyl_mom_calc_advance(recyc->m0op_neut, &recyc->emit_normal_r, &recyc->emit_cbuff_r,
    recyc->unit_phase_flux_neut, recyc->emit_flux);
  app->stat.neut_species_diag_calc_tm += gkyl_time_diff_now_sec(wst);

  struct timespec wtm = gkyl_wall_clock();
  // Copy to skin to write it out.
  gkyl_array_clear(recyc->diag_out, 0.0);
  gkyl_array_copy_range_to_range(recyc->diag_out, recyc->emit_flux, cskin_r, &recyc->emit_cbuff_r);
  if (app->use_gpu)
    gkyl_array_copy(recyc->diag_out_ho, recyc->diag_out);
  
  struct gkyl_msgpack_map_elem desc1[] = {
    { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Emitted boundary particle flux." }
  };
  const struct gkyl_msgpack_map_elem* io_meta[] = {s->io_meta_conf, app->gk_geom->io_meta_basic, desc1};
  struct gkyl_msgpack_data *mt1 = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

  const char *fmt = "%s-%s_recycling_%s%s_%s_flux_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, s->info.name, vars[dir], edge[edi], s->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, s->info.name, vars[dir], edge[edi], s->info.name, frame);
  
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt1, recyc->diag_out_ho, fileNm);
  app->stat.neut_species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
  app->stat.n_diag_io += 1;

  gkyl_msgpack_data_release(mt1); 
}

static void
gk_neut_species_recycle_write_flux_disabled(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_recycle_wall *recyc, double tm, int frame)
{
}

void
gk_neut_species_recycle_write_flux(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_recycle_wall *recyc, double tm, int frame)
{
  recyc->write_flux_func(app, s, recyc, tm, frame);
}

struct gk_neut_recycling_maxwellian_params {
  double temp; // Temperature of the neutral species emitted during recycling.
};

static void
gk_neut_recycling_maxwellian_den(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  fout[0] = 1.0;
}

static void
gk_neut_recycling_maxwellian_udrift(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

static void
gk_neut_recycling_maxwellian_temp(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct gk_neut_recycling_maxwellian_params *params = ctx;
  fout[0] = params->temp;
}

void
gk_neut_species_recycle_init(struct gkyl_gyrokinetic_app *app, struct gk_recycle_wall *recyc,
  int dir, enum gkyl_edge_loc edge, struct gkyl_gyrokinetic_emission_inp *params,
  struct gk_neut_species *s, bool use_gpu)
{
  recyc->params = params;
  recyc->num_species = params->num_species;
  recyc->edge = edge;
  recyc->dir = dir;
  recyc->write_diagnostics = edge == GKYL_LOWER_EDGE? s->lower_bc[dir].write_diagnostics : s->upper_bc[dir].write_diagnostics;

  int cdim = app->cdim;
  int ndim = app->cdim + s->info.vdim;

  int e = recyc->edge == GKYL_LOWER_EDGE? 0 : 1;

  // Create boundary grids and ranges.
  int cells[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  for (int i=0; i<ndim; ++i) {
    cells[i] = s->grid.cells[i];
    lower[i] = s->grid.lower[i];
    upper[i] = s->grid.upper[i];
  }
  cells[dir] = 1;
  lower[dir] = e==0? s->grid.lower[dir] - s->grid.dx[dir] : s->grid.upper[dir];
  upper[dir] = e==0? s->grid.lower[dir] : s->grid.upper[dir] + s->grid.dx[dir];
  gkyl_rect_grid_init(&recyc->emit_grid, ndim, lower, upper, cells);

  recyc->emit_ghost_r = e==0? &s->local_lower_ghost[dir] : &s->local_upper_ghost[dir];
  recyc->emit_skin_r = e==0? &s->local_lower_skin[dir] : &s->local_upper_skin[dir];
  gkyl_range_init(&recyc->emit_buff_r, ndim, recyc->emit_ghost_r->lower, recyc->emit_ghost_r->upper);
  gkyl_range_init(&recyc->emit_cbuff_r, cdim, recyc->emit_ghost_r->lower, recyc->emit_ghost_r->upper);

  // Buffer for scaled Maxwellian in ghost.
  recyc->bc_buffer = mkarr(app->use_gpu, s->basis.num_basis, recyc->emit_skin_r->volume);
  // Initialize fixed func bc object to project the unit Maxwellian in ghost
  struct gkyl_bc_basic *bc_basic_op = gkyl_bc_basic_new(dir, edge, GKYL_BC_FIXED_FUNC, s->basis_on_dev,
    recyc->emit_skin_r, recyc->emit_ghost_r, s->f->ncomp, app->cdim, app->use_gpu);
  // Project unit Maxwellian.
  struct gk_neut_recycling_maxwellian_params neut_max_pars = {
    .temp = e==0? s->lower_bc[dir].emission.emission_temp : s->upper_bc[dir].emission.emission_temp,
  };
  struct gkyl_gyrokinetic_projection recyc_proj_inp = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
    .ctx_density = &neut_max_pars,
    .density = gk_neut_recycling_maxwellian_den,
    .ctx_upar = &neut_max_pars,
    .udrift = gk_neut_recycling_maxwellian_udrift,
    .ctx_temp = &neut_max_pars,
    .temp = gk_neut_recycling_maxwellian_temp,
  };
  struct gk_proj proj_unit_maxwellian;
  gk_neut_species_projection_init(app, s, recyc_proj_inp, &proj_unit_maxwellian);
  gk_neut_species_projection_calc(app, s, &proj_unit_maxwellian, s->f1, 0.0); // Temporarily use f1.

  // Calculate flux associated with unit Maxwellian projected in f0.
  int num_eqns = 0;
  if (s->collisionless.collisionless_id == GKYL_GK_COLLISIONLESS_NEUTRAL)
    num_eqns += 1; // Collisionless terms.

  const struct gkyl_dg_eqn **eqns = gkyl_malloc(num_eqns*sizeof(struct gkyl_dg_eqn *));

  int eqc = 0;
  if (s->collisionless.collisionless_id == GKYL_GK_COLLISIONLESS_NEUTRAL)
    eqns[eqc++] = gkyl_dg_eqn_acquire(s->collisionless.eqn_vlasov);

  recyc->f0_flux_slvr = gkyl_boundary_flux_new(recyc->dir, recyc->edge, &s->grid,
    recyc->emit_skin_r, recyc->emit_ghost_r, num_eqns, eqns, app->use_gpu);  

  for (int i=0; i<num_eqns; i++)
    gkyl_dg_eqn_release(eqns[i]);

  gkyl_free(eqns);

  // The from-flux Vlasov equation's boundary term reads the stored conf-space
  // surface fluxes, so compute them for the unit Maxwellian first (cflrate is
  // scratch here; it is cleared at the start of every RHS evaluation).
  gkyl_dg_vlasov_conf_flux_surf_advance(s->collisionless.calc_conf_flux, &app->local,
    &s->local, &s->local_ext, s->conf_poisson_tensor, s->hamil, s->f1,
    s->cflrate, s->collisionless.conf_flux_surf);
  gkyl_boundary_flux_advance(recyc->f0_flux_slvr, s->f1, s->f1);
  recyc->unit_phase_flux_neut = mkarr(app->use_gpu, s->basis.num_basis, recyc->emit_buff_r.volume);
  gkyl_array_copy_range_to_range(recyc->unit_phase_flux_neut, s->f1, &recyc->emit_buff_r,
    recyc->emit_ghost_r);

  recyc->write_flux_func = gk_neut_species_recycle_write_flux_disabled;
  if (recyc->write_diagnostics) {
    recyc->write_flux_func = gk_neut_species_recycle_write_flux_enabled;
  }

  gkyl_bc_basic_buffer_fixed_func(bc_basic_op, recyc->bc_buffer, s->f1);
  gkyl_array_clear(s->f1, 0.0);

  gkyl_bc_basic_release(bc_basic_op);
  gk_neut_species_projection_release(app, &proj_unit_maxwellian);
}

void
gk_neut_species_recycle_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_recycle_wall *recyc)
{
  int cdim = app->cdim;
  int vdim = s->info.vdim;
 
  // Define necessary grid, ranges, and array for calculating the desired Maxwellian for ghost.
  recyc->f_emit = mkarr(app->use_gpu, s->basis.num_basis, recyc->emit_buff_r.volume);

  int ghost[GKYL_MAX_DIM]; // Number of ghost cells in each direction.
  for (int d=0; d<cdim; ++d)
    ghost[d] = 1;
  for (int d=0; d<vdim; ++d)
    ghost[cdim+d] = 0;

  // Calculate the flux associated with unit-density Maxwellian.
  gkyl_bc_emission_flux_ranges(&recyc->emit_normal_r, recyc->dir + cdim, &recyc->emit_buff_r,
    ghost, recyc->edge);
  
  recyc->unit_m0_flux_neut = mkarr(app->use_gpu, app->basis.num_basis, recyc->emit_cbuff_r.volume);

  struct gkyl_mom_vlasov_inp inp_mom = {
    .conf_basis = &app->basis,
    .phase_basis = &s->basis,
    .vel_range = &s->local_vel,
    .vel_map = s->vlasov_vel_map,
    .hamil_range = &s->hamil_range,
    .hamil = s->hamil,
    .model_id = s->model_id,
    .hamil_id = s->hamil_id,
    .mom_type = GKYL_F_MOMENT_M0,
    .use_gpu = app->use_gpu,
  };
  recyc->m0op_neut_type = gkyl_mom_vlasov_inew(&inp_mom);
  recyc->m0op_neut = gkyl_mom_calc_new(&recyc->emit_grid, recyc->m0op_neut_type, app->use_gpu);

  gkyl_mom_calc_advance(recyc->m0op_neut, &recyc->emit_normal_r,
    &recyc->emit_cbuff_r, recyc->unit_phase_flux_neut, recyc->unit_m0_flux_neut);

  // Define memory for div bin op for calculating correct scaling factor.
  if (app->use_gpu)
    recyc->mem_geo = gkyl_dg_bin_op_mem_cu_dev_new(recyc->emit_cbuff_r.volume, app->basis.num_basis);
  else
    recyc->mem_geo = gkyl_dg_bin_op_mem_new(recyc->emit_cbuff_r.volume, app->basis.num_basis);
  
  for (int i=0; i<recyc->num_species; ++i) {
    recyc->impact_species[i] = gk_find_species(app, recyc->params->in_species[i]);
    struct gk_species *gks = recyc->impact_species[i];

    int e = recyc->edge == GKYL_LOWER_EDGE? 0 : 1;

    // Create boundary grids and ranges for impacting species.
    int cells[GKYL_MAX_DIM];
    double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
    for (int i=0; i<cdim+gks->info.vdim; ++i) {
      cells[i] = gks->grid.cells[i];
      lower[i] = gks->grid.lower[i];
      upper[i] = gks->grid.upper[i];
    }
    cells[recyc->dir] = 1;
    lower[recyc->dir] = e==0? gks->grid.lower[recyc->dir] - gks->grid.dx[recyc->dir] : gks->grid.upper[recyc->dir];
    upper[recyc->dir] = e==0? gks->grid.lower[recyc->dir] : gks->grid.upper[recyc->dir] + gks->grid.dx[recyc->dir];
    gkyl_rect_grid_init(&recyc->impact_grid[i], cdim+gks->info.vdim, lower, upper, cells);

    struct gkyl_range *phase_skin_r = e==0? &gks->local_lower_skin[recyc->dir] : &gks->local_upper_skin[recyc->dir];
    struct gkyl_range *phase_ghost_r = e==0? &gks->local_lower_ghost[recyc->dir] : &gks->local_upper_ghost[recyc->dir];
    recyc->impact_ghost_r[i] = phase_ghost_r;
    gkyl_range_init(&recyc->impact_buff_r[i], cdim+gks->info.vdim, phase_ghost_r->lower, phase_ghost_r->upper);
    gkyl_range_init(&recyc->impact_cbuff_r[i], cdim, phase_ghost_r->lower, phase_ghost_r->upper);

    recyc->phase_flux_gk[i] = mkarr(app->use_gpu, gks->basis.num_basis, phase_ghost_r->volume);
    recyc->m0_flux_gk[i] = mkarr(app->use_gpu, app->basis.num_basis, recyc->impact_cbuff_r[i].volume);    
    
    recyc->spectrum[i] = mkarr(app->use_gpu, s->basis.num_basis, recyc->emit_buff_r.volume);

    recyc->m0op_gk[i] = gkyl_dg_updater_moment_gyrokinetic_new(&recyc->impact_grid[i], &app->basis,
      &gks->basis, &recyc->emit_cbuff_r, gks->info.mass, gks->info.charge, gks->vel_map,
      app->gk_geom, 0, GKYL_F_MOMENT_M0, false, app->use_gpu);
    
    // Create a phase-space range that only includes velocities towards the boundary.
    gkyl_bc_emission_flux_ranges(&recyc->impact_normal_r[i], recyc->dir + cdim, &recyc->impact_buff_r[i],
      ghost, recyc->edge);
  }

  // For writing diagnostics, if needed.
  if (recyc->write_diagnostics) {
    recyc->f_diag = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);
    recyc->emit_flux = mkarr(app->use_gpu, app->basis.num_basis, recyc->emit_cbuff_r.volume);
    recyc->diag_out = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    recyc->diag_out_ho = recyc->diag_out;
    if (app->use_gpu)
      recyc->diag_out_ho = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  }
   
}

void
gk_neut_species_recycle_apply_bc(struct gkyl_gyrokinetic_app *app, const struct gk_recycle_wall *recyc,
  const struct gk_neut_species *s, struct gkyl_array *fout)
{
  gkyl_array_clear(recyc->f_emit, 0.0); // Zero emitted distribution before beginning accumulate

  // Inelastic emission contribution.
  // This relies on the calculation of the ion flux (phase_flux_gk).
  // from gk_species_bflux_rhs_solver.
  for (int i=0; i<recyc->num_species; ++i) {
    // Copy unit-density Maxwellian from buffer.
    gkyl_array_set(recyc->spectrum[i], recyc->params->recycling_frac, recyc->bc_buffer);

    struct gk_species *gks = recyc->impact_species[i];

    // Calculate M0 moment of ion flux.
    gk_species_bflux_get_flux(&gks->bflux, recyc->dir, recyc->edge, recyc->phase_flux_gk[i], &recyc->impact_buff_r[i]);
    gkyl_dg_updater_moment_gyrokinetic_advance(recyc->m0op_gk[i], &recyc->impact_normal_r[i],
      &recyc->impact_cbuff_r[i], recyc->phase_flux_gk[i], recyc->m0_flux_gk[i]);

    // Calculate scaling factor from ratio of ion flux to unit-density flux.
    gkyl_dg_div_op_range(recyc->mem_geo, &app->basis, 0, recyc->m0_flux_gk[i], 0, recyc->m0_flux_gk[i],
      0, recyc->unit_m0_flux_neut, &recyc->emit_cbuff_r);

    gkyl_dg_mul_conf_phase_op_accumulate_range(&app->basis, &s->basis, recyc->f_emit, 1.0,
      recyc->m0_flux_gk[i], recyc->spectrum[i], &recyc->impact_cbuff_r[i], &recyc->emit_buff_r);    
  }

  gkyl_array_set_range_to_range(fout, 1.0, recyc->f_emit, recyc->emit_ghost_r,
    &recyc->emit_buff_r);
}

void
gk_neut_species_recycle_release(const struct gkyl_gyrokinetic_app *app, const struct gk_recycle_wall *recyc)
{
  gkyl_boundary_flux_release(recyc->f0_flux_slvr);
  gkyl_array_release(recyc->unit_phase_flux_neut);
  gkyl_array_release(recyc->f_emit);

  if (recyc->write_diagnostics) {
    gkyl_array_release(recyc->f_diag);
    gkyl_array_release(recyc->emit_flux);
    gkyl_array_release(recyc->diag_out);
    if (app->use_gpu)
      gkyl_array_release(recyc->diag_out_ho);
  }

  gkyl_array_release(recyc->unit_m0_flux_neut);
  gkyl_mom_type_release(recyc->m0op_neut_type);
  gkyl_mom_calc_release(recyc->m0op_neut);
  gkyl_dg_bin_op_mem_release(recyc->mem_geo);

  for (int i=0; i<recyc->num_species; ++i) {
    gkyl_array_release(recyc->phase_flux_gk[i]);
    gkyl_array_release(recyc->m0_flux_gk[i]);
    gkyl_array_release(recyc->spectrum[i]);
    gkyl_dg_updater_moment_gyrokinetic_release(recyc->m0op_gk[i]);
  }

  gkyl_array_release(recyc->bc_buffer);
}
