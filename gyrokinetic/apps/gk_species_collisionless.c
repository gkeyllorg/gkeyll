#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_species_collisionless_flux_disabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin)
{
  // Do nothing.
}

static void
gk_species_collisionless_flux_enabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin)
{
  // Compute the surface expansion of the phase space flux
  // Note: Each cell stores the *lower* surface expansions of the 
  // phase space flux, so local_ext range needed to index the output
  // values of flux_surf even though we only loop over local ranges
  // to avoid evaluating quantities such as geometry in ghost cells
  // where they are not defined.
  gkyl_gk_collisionless_flux_surf(gkcls->surf_flux_op, 
    &app->local, &species->local, &app->local_ext, &species->local_ext, 
    species->gyro_phi, fin, gkcls->yfield, gkcls->flux_surf, species->cflrate);
}

static void
gk_species_collisionless_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
}

static void
gk_species_collisionless_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  gkcls->flux_func(app, species, gkcls, fin);

  gkyl_dg_updater_gyrokinetic_advance(gkcls->slvr, &species->local, 
    fin, species->cflrate, rhs);

  gkcls->fdot_scaling(app, species, gkcls, rhs, species->cflrate, &species->local);

  app->stat.species_collisionless_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_collisionless_fdot_scaling_disabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_collisionless *gkcls, struct gkyl_array *rhs, struct gkyl_array *cflrate, struct gkyl_range *rng)
{
  // Do nothing.
}

static void
gk_species_collisionless_fdot_scaling_enabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_collisionless *gkcls, struct gkyl_array *rhs, struct gkyl_array *cflrate, struct gkyl_range *rng)
{
  gkyl_array_scale_range(rhs, gkcls->scale_fac, rng);
  gkyl_array_scale_range(cflrate, gkcls->scale_fac, rng);
}

static void
gk_species_collisionless_write_diags_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_collisionless *gkcls, double tm, int frame)
{
  // Do nothing.
}

static void
gk_species_collisionless_write_diags_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_collisionless *gkcls, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();

  // Package metadata.
  gkyl_msgpack_map_elem_set_double(gks->io_meta_phase_len, gks->io_meta_phase, "time", tm);
  gkyl_msgpack_map_elem_set_uint(gks->io_meta_phase_len, gks->io_meta_phase, "frame", frame);
  struct gkyl_msgpack_map_elem desc[] = {
    { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Collisionless flux at cell surface." }
  };
  int io_meta_len[] = {gks->io_meta_phase_len, app->gk_geom->io_meta_basic_len, 1};
  const struct gkyl_msgpack_map_elem* io_meta[] = {gks->io_meta_phase, app->gk_geom->io_meta_basic, desc};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

  // Write gkcls->flux_surf
  const char *fmt = "%s-%s_collisionless_surf_flux_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);
  gkyl_array_copy(gkcls->flux_surf_ho, gkcls->flux_surf);
  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt,
    gkcls->flux_surf_ho, fileNm);

  gkyl_msgpack_data_release(mt);
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_io += 1;
}

void 
gk_species_collisionless_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, 
  struct gk_collisionless *gkcls)
{
  gkcls->collisionless_id = gks->info.collisionless.type;
  gkcls->write_diagnostics = gks->info.collisionless.write_diagnostics;

  gkcls->write_diags_func = gk_species_collisionless_write_diags_disabled;
  gkcls->flux_func = gk_species_collisionless_flux_disabled;
  gkcls->rhs_func = gk_species_collisionless_rhs_disabled;

  if (gkcls->collisionless_id) {

    int cdim = app->cdim, vdim = gks->info.vdim;
    int pdim = cdim+vdim;

    // Determine which directions are zero-flux. By default
    // we do not have zero-flux boundary conditions in any direction.
    bool is_zero_flux[2*GKYL_MAX_DIM] = {false};
    for (int dir=0; dir<app->cdim; ++dir) {
      if (gks->lower_bc[dir].type == GKYL_BC_GK_SPECIES_ZERO_FLUX)
        is_zero_flux[dir] = true;
      if (gks->upper_bc[dir].type == GKYL_BC_GK_SPECIES_ZERO_FLUX)
        is_zero_flux[dir+pdim] = true;
    }

    // Need to figure out size of flux_surf by finding size of surface basis set 
    struct gkyl_basis surf_basis;
    struct gkyl_basis surf_vpar_basis;
    // Define surface bases
    if (app->poly_order > 1) {
      gkyl_cart_modal_serendip(&surf_basis, pdim-1, app->poly_order);
    }
    else {
      gkyl_cart_modal_serendip(&surf_vpar_basis, pdim-1, app->poly_order);
      if (vdim>1) {
        gkyl_cart_modal_gkhybrid(&surf_basis, cdim-1, vdim); // p=2 in vparallel
      }
      else {
        gkyl_cart_modal_serendip(&surf_basis, pdim-1, 2); // p=2 in vparallel
      }
    }
    int flux_surf_sz = (cdim)*surf_basis.num_basis + surf_vpar_basis.num_basis;

    // Allocate arrays to store surface phase space flux.
    gkcls->yfield = mkarr(app->use_gpu, cdim*gks->basis.num_basis, gks->local_ext.volume);
    gkcls->flux_surf = mkarr(app->use_gpu, flux_surf_sz, gks->local_ext.volume);

    if (gkcls->collisionless_id == GKYL_GK_COLLISIONLESS_EM_BPERP) {
      // Parallel component of magnetic vector potential.
      gkcls->apar = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      gkcls->apardot = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);    
    }
    else {
      gkcls->apar    = gkyl_array_acquire(app->field->phi_smooth); // Not used.
      gkcls->apardot = gkyl_array_acquire(app->field->phi_smooth); // Not used.
    }

    enum gkyl_gyrokinetic_bc_type bctype_conf[2*GKYL_MAX_CDIM];
    bool is_mpi_edge[2*GKYL_MAX_CDIM] = { false };
    for (int d=0; d<app->cdim; d++) {
      // Only global domain edges use the configured physical boundary
      // stencil. Mark local edges inside the domain as MPI interfaces so the
      // surface-flux updater uses the ordinary two-sided interior stencil.
      is_mpi_edge[d] = app->local.lower[d] != app->global.lower[d];
      is_mpi_edge[GKYL_MAX_CDIM+d] = app->local.upper[d] != app->global.upper[d];
      bctype_conf[d] = !is_mpi_edge[d]
        ? gks->lower_bc[d].type : GKYL_BC_GK_SPECIES_PERIODIC;
      bctype_conf[GKYL_MAX_CDIM+d] = !is_mpi_edge[GKYL_MAX_CDIM+d]
        ? gks->upper_bc[d].type : GKYL_BC_GK_SPECIES_PERIODIC;
    }

    gkcls->surf_flux_op = gkyl_gk_collisionless_flux_new(&gks->grid, &app->basis, &gks->basis, 
      gks->info.charge, gks->info.mass,
      gkcls->collisionless_id, app->gk_geom, 
      app->dg_geom, app->gk_dg_geom, gks->vel_map, bctype_conf, is_mpi_edge, app->use_gpu);

    struct gkyl_dg_gyrokinetic_auxfields aux_inp = { .flux_surf = gkcls->flux_surf, 
      .phi = gks->gyro_phi, .apar = gkcls->apar, .apardot = gkcls->apardot,
      .yfield = gkcls->yfield, };
    // Create solver.
    gkcls->slvr = gkyl_dg_updater_gyrokinetic_new(&gks->grid, &app->basis, &gks->basis, 
      &app->local, &gks->local, is_zero_flux, gks->info.charge, gks->info.mass,
      gkcls->collisionless_id, app->gk_geom, gks->vel_map, 
      &aux_inp, app->use_gpu);

    gkcls->scale_fac = 1.0; // Not used if scale_factor in input file is not given.
    gkcls->fdot_scaling = gk_species_collisionless_fdot_scaling_disabled;
    if (1.0e-16 < fabs(gks->info.collisionless.scale_factor)) {
      gkcls->scale_fac = gks->info.collisionless.scale_factor;
      gkcls->fdot_scaling = gk_species_collisionless_fdot_scaling_enabled;
    }

    // Other methods chosen at runtime.
    gkcls->flux_func = gk_species_collisionless_flux_enabled;
    gkcls->rhs_func = gk_species_collisionless_rhs_enabled;
    if (gkcls->write_diagnostics) {
      gkcls->flux_surf_ho = mkarr(false, gkcls->flux_surf->ncomp, gkcls->flux_surf->size);
      gkcls->write_diags_func = gk_species_collisionless_write_diags_enabled;
    }
  }
}

void
gk_species_collisionless_flux(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin)
{
  gkcls->flux_func(app, species, gkcls, fin);
}

void
gk_species_collisionless_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkcls->rhs_func(app, species, gkcls, fin, rhs);
}

void
gk_species_collisionless_write_diags(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_collisionless *gkcls, double tm, int frame)
{
  gkcls->write_diags_func(app, gks, gkcls, tm, frame);
}

void
gk_species_collisionless_release(const struct gkyl_gyrokinetic_app *app, const struct gk_collisionless *gkcls)
{
  if (gkcls->collisionless_id) {

    gkyl_array_release(gkcls->yfield);
    gkyl_array_release(gkcls->flux_surf);
    gkyl_array_release(gkcls->apar);
    gkyl_array_release(gkcls->apardot);
  
    gkyl_gk_collisionless_flux_release(gkcls->surf_flux_op);
    gkyl_dg_updater_gyrokinetic_release(gkcls->slvr);

    if (gkcls->write_diagnostics) {
      gkyl_array_release(gkcls->flux_surf_ho);
    }
  }
}

void
gk_species_collisionless_reset(gkyl_gyrokinetic_app* app, double tm, struct gk_species *gks,
  struct gk_collisionless *gkcls, struct gkyl_gyrokinetic_collisionless gkcls_inp)
{
  gkcls->scale_fac = 1.0;
  gkcls->fdot_scaling = gk_species_collisionless_fdot_scaling_disabled;
  if (1.0e-16 < fabs(gkcls_inp.scale_factor)) {
    gks->info.collisionless.scale_factor = gkcls_inp.scale_factor;

    gkcls->scale_fac = gkcls_inp.scale_factor;
    gkcls->fdot_scaling = gk_species_collisionless_fdot_scaling_enabled;
  }
}
