#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_neut_species_scaling_cross_moms_enabled(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_scaling *sca, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  struct timespec wst = gkyl_wall_clock();    

  struct gk_species *gks_elc = &app->species[sca->elc_idx]; 

  // Compute electron Maxwellian moments (J*n, u_par, T/m).
  gk_species_moment_calc(&gks_elc->lte.moms, 
    gks_elc->local, app->local, fin[sca->elc_idx]);

  // Divide the electron density by the Jacobian.
  gkyl_dg_div_op_range(gks_elc->lte.moms.mem_geo, &app->basis, 0, gks_elc->lte.moms.marr,
    0, gks_elc->lte.moms.marr, 0, app->gk_geom->geo_int.jacobgeo, &app->local); 

  // Compute ionization reactivity <sigma v>_iz.
  gkyl_dg_iz_coll(sca->iz_react_calc, gks_elc->lte.moms.marr, 
    sca->dndt_react, sca->dndt_react, sca->reactivity, 0);
  
  app->stat.neut_species_react_mom_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_neut_species_scaling_cross_moms_disabled(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_scaling *sca, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  // Do nothing.
}

static void
gk_neut_species_scaling_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_scaling *sca, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();  

  struct gk_species *gks_elc = &app->species[sca->elc_idx];

  // Compute (J*n_neut)*n_elc*<sigma v>_iz.
  gkyl_dg_mul_op_range(&app->basis, 0, sca->dndt_react,
    0, sca->Jm0_init, 0, sca->reactivity, &app->local);  
  gkyl_dg_mul_op_range(&app->basis, 0, sca->dndt_react,
    0, gks_elc->lte.moms.marr, 0, sca->dndt_react, &app->local);  

  // Volume integrate the reaction contribution.
  gkyl_array_integrate_advance(sca->integrate_op, sca->dndt_react, 1.0, 0, &app->local, 0, sca->react_vol_integ_local);
  // Reduce over MPI processes.
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, 
    sca->react_vol_integ_local, sca->react_vol_integ_global);
  if (app->use_gpu)
    gkyl_cu_memcpy(&sca->react_vol_integ, sca->react_vol_integ_global, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(&sca->react_vol_integ, sca->react_vol_integ_global, sizeof(double));

  app->stat.neut_species_react_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_neut_species_scaling_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_scaling *sca, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  // Do nothing.
}

static void
gk_neut_species_scaling_apply_enabled(gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
  struct gk_scaling *sca, struct gkyl_array *fin, struct gkyl_array **bflux[])
{
  struct gk_species *gks_ion = &app->species[sca->ion_idx];
  struct gkyl_array **bflux_ion = bflux[sca->ion_idx];

  if (sca->react_vol_integ > 0.0) {
    gkyl_array_clear(sca->dndt_react, 0.0);

    double bflux_intm0_local_ho = 0.0;
    for (int j=0; j<sca->num_boundaries; ++j) {
      // Add integrated M0 moments of boundary fluxes.
      gk_species_bflux_get_flux_mom(&gks_ion->bflux, sca->boundaries_dir[j], sca->boundaries_edge[j],
        GKYL_F_MOMENT_M0, bflux_ion, sca->dndt_react, &sca->boundaries_conf_ghost[j]);
      gkyl_array_integrate_advance(sca->integrate_op, sca->dndt_react, 1.0, 0,
        &sca->boundaries_conf_ghost[j], 0, sca->bflux_m0_vol_integ_local);

      double bflux_m0_vol_integ_local_ho;
      if (app->use_gpu)
        gkyl_cu_memcpy(&bflux_m0_vol_integ_local_ho, sca->bflux_m0_vol_integ_local, sizeof(double), GKYL_CU_MEMCPY_D2H);
      else
        memcpy(&bflux_m0_vol_integ_local_ho, sca->bflux_m0_vol_integ_local, sizeof(double));

      bflux_intm0_local_ho += bflux_m0_vol_integ_local_ho;
    }

    double bflux_intm0_global_ho;
    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &bflux_intm0_local_ho, &bflux_intm0_global_ho);

    double neut_scaling_fac = sca->recycling_coeff * bflux_intm0_global_ho / sca->react_vol_integ;

    // Divide by the present J*rho, and multiply by neut_scaling_fac*mass*Jm0_init.
    gkyl_array_set_offset_range(sca->dndt_react, 1.0, fin, 0, &app->local);
    for (int i=0; i<ns->num_moments; ++i)
      gkyl_dg_div_op_range(gks_ion->lte.moms.mem_geo, &app->basis, i, fin,
        i, fin, 0, sca->dndt_react, &app->local); 
  
    for (int i=0; i<ns->num_moments; ++i)
      gkyl_dg_mul_op_range(&app->basis, i, fin,
        i, fin, 0, sca->Jm0_init, &app->local);  

    gkyl_array_scale(fin, neut_scaling_fac*ns->info.mass);
  }
}

static void
gk_neut_species_scaling_apply_disabled(gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
  struct gk_scaling *sca, struct gkyl_array *fin, struct gkyl_array **bflux[])
{
  // Do nothing.
}

static void
gk_neut_species_scaling_write_enabled(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_scaling *sca, int ridx, double tm, int frame)
{
  // React diagnostics usually written from gk_species.
  // In the case of static gk_species, write_diagnostics flag
  // can be used to check reaction rates from gk_neut_species.
  struct timespec wst = gkyl_wall_clock();
//  // Compute reaction rate
//  const struct gkyl_array *fin[app->num_species];
//  const struct gkyl_array *fin_neut[app->num_neut_species];
//  for (int i=0; i<app->num_species; ++i) {
//    fin[i] = app->species[i].f;
//  }
//  for (int i=0; i<app->num_neut_species; ++i) {
//    fin_neut[i] = app->neut_species[i].f;
//  }
//  gk_neut_species_scaling_cross_moms(app, gkns, gkr, fin, fin_neut);
  app->stat.neut_species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
  
  struct timespec wtm = gkyl_wall_clock();
  app->stat.neut_species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
}

static void
gk_neut_species_scaling_write_disabled(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_scaling *sca, int ridx, double tm, int frame)
{
  // Do nothing
}

void 
gk_neut_species_scaling_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns, 
  struct gk_scaling *sca)
{
  struct gkyl_gyrokinetic_scaling_inp *sca_inp = &ns->info.scaling;

  sca->type = 0;
  sca->num_boundaries = 0;
  sca->cross_moms_func_neut = gk_neut_species_scaling_cross_moms_disabled;
  sca->rhs_func_neut = gk_neut_species_scaling_rhs_disabled;
  sca->apply_func_neut = gk_neut_species_scaling_apply_disabled;
  sca->write_func_neut = gk_neut_species_scaling_write_disabled;

  if (sca_inp->type == GKYL_GK_SPECIES_SCALING_RECYCLING_IZ_BALANCE) {
    assert(ns->is_fluid);

    sca->type = sca_inp->type;
    sca->num_boundaries = sca_inp->num_boundaries;
    sca->recycling_coeff = sca_inp->recycling_coeff;
    sca->write_diagnostics = sca_inp->write_diagnostics;

    // Initial number density times conf-space Jacobian.
    sca->Jm0_init = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

    // Create an updater that integrates an array.
    sca->integrate_op = gkyl_array_integrate_new(&app->grid, &app->basis, 1, GKYL_ARRAY_INTEGRATE_OP_NONE, app->use_gpu);
    int num_mom = 1;
    if (app->use_gpu){
      sca->react_vol_integ_local  = gkyl_cu_malloc(sizeof(double[num_mom]));
      sca->react_vol_integ_global = gkyl_cu_malloc(sizeof(double[num_mom]));
      sca->bflux_m0_vol_integ_local = gkyl_cu_malloc(sizeof(double));
    }
    else {
      sca->react_vol_integ_local  = gkyl_malloc(sizeof(double[num_mom]));
      sca->react_vol_integ_global = gkyl_malloc(sizeof(double[num_mom]));
      sca->bflux_m0_vol_integ_local = gkyl_malloc(sizeof(double));
    }

    sca->cross_moms_func_neut = gk_neut_species_scaling_cross_moms_enabled;
    sca->rhs_func_neut = gk_neut_species_scaling_rhs_enabled;
    sca->apply_func_neut = gk_neut_species_scaling_apply_enabled;
    if (sca->write_diagnostics)
      sca->write_func_neut = gk_neut_species_scaling_write_enabled;
    else
      sca->write_func_neut = gk_neut_species_scaling_write_disabled;
  }
}

void 
gk_neut_species_scaling_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
  struct gk_scaling *sca)
{
  if (sca->type == GKYL_GK_SPECIES_SCALING_RECYCLING_IZ_BALANCE) {
    struct gkyl_gyrokinetic_scaling_inp *sca_inp = &ns->info.scaling;

    // Fetch index of species for indexing arrays.
    sca->elc_idx = gk_find_species_idx(app, sca_inp->electron_name);
    sca->ion_idx = gk_find_species_idx(app, sca_inp->impacting_ion_name);

    struct gk_species *gks_ion = &app->species[sca->ion_idx];

    for (int j=0; j<sca->num_boundaries; ++j) {
      int dir  = sca_inp->boundaries_dir[j];
      int edge = sca_inp->boundaries_edge[j];
  
      // Source adaptation on periodic, zero flux, or reflect boundary is not allowed.
      assert(gks_ion->bc_is_np[dir]);
      if (edge == GKYL_LOWER_EDGE) {
        assert(gks_ion->lower_bc[dir].type != GKYL_BC_GK_SPECIES_ZERO_FLUX);
        assert(gks_ion->lower_bc[dir].type != GKYL_BC_GK_SPECIES_REFLECT);
      } else {
        assert(gks_ion->upper_bc[dir].type != GKYL_BC_GK_SPECIES_ZERO_FLUX);
        assert(gks_ion->upper_bc[dir].type != GKYL_BC_GK_SPECIES_REFLECT);
      }
  
      // Default scenario: we set the ranges to the full range of the ghost cells.
      sca->boundaries_conf_ghost[j] = edge == GKYL_LOWER_EDGE ? app->local_lower_ghost[dir] : app->local_upper_ghost[dir];
      sca->boundaries_dir[j]  = dir;
      sca->boundaries_edge[j] = edge;
  
      // Specific scenario if we are in a inner wall limited case. We select only SOL range in parallel direction.
      if (dir == app->cdim-1 && app->gk_geom->has_LCFS) {
        sca->boundaries_conf_ghost[j] = edge == GKYL_LOWER_EDGE? app->local_lower_ghost_par_sol
                                                               : app->local_upper_ghost_par_sol;
      }
    }

    struct gkyl_dg_iz_inp iz_inp = {
      .cbasis = &app->basis, 
      .conf_rng = &app->local, 
      .type_ion = sca_inp->impacting_ion_id, 
      .charge_state = 0,
      .type_self = GKYL_SELF_ION, // Could be GKYL_SELF_DONOR. It just can't be
                                  // GKYL_SELF_ELC because we don't need to
                                  // compute the ionization temperatures.
    };
    sca->iz_react_calc = gkyl_dg_iz_new(&iz_inp, app->use_gpu);

    // Reaction contribution.
    sca->reactivity = mkarr(app->use_gpu, ns->basis.num_basis, ns->local_ext.volume);
    sca->dndt_react = mkarr(app->use_gpu, ns->basis.num_basis, ns->local_ext.volume);
  }
}

void 
gk_neut_species_scaling_apply_ic_cross(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
  struct gk_scaling *sca)
{
  if (sca->type == GKYL_GK_SPECIES_SCALING_RECYCLING_IZ_BALANCE) {
    // Project and store the initial state. Need to re-project so that restarts use this as the state at t=0.
    if (ns->info.init_from_file.type == 0)
      gk_neut_species_projection_calc(app, ns, &ns->proj_init, ns->f1, 0.0);

    gkyl_array_set_offset(sca->Jm0_init, 1.0/ns->info.mass, ns->f1, 0);

    if (ns->info.init_from_file.type == 0)
      gkyl_array_clear(ns->f1, 0.0);
  }
}

void
gk_neut_species_scaling_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_neut_species *ns,
  struct gk_scaling *sca, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  sca->cross_moms_func_neut(app, ns, sca, fin, fin_neut);
}

void
gk_neut_species_scaling_rhs(gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
  struct gk_scaling *sca, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  sca->rhs_func_neut(app, ns, sca, fin, rhs);
}

void
gk_neut_species_scaling_apply(gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
  struct gk_scaling *sca, struct gkyl_array *fin, struct gkyl_array **bflux[])
{
  sca->apply_func_neut(app, ns, sca, fin, bflux);
}

void
gk_neut_species_scaling_write(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_scaling *sca, int ridx, double tm, int frame)
{
  sca->write_func_neut(app, gkns, sca, ridx, tm, frame);
}

void 
gk_neut_species_scaling_release(const struct gkyl_gyrokinetic_app *app,
  const struct gk_scaling *sca)
{
  if (sca->type == GKYL_GK_SPECIES_SCALING_RECYCLING_IZ_BALANCE) {
    gkyl_array_release(sca->Jm0_init);
    gkyl_array_integrate_release(sca->integrate_op);

    if (app->use_gpu){
      gkyl_cu_free(sca->react_vol_integ_local );
      gkyl_cu_free(sca->react_vol_integ_global);
      gkyl_cu_free(sca->bflux_m0_vol_integ_local);
    }
    else {
      gkyl_free(sca->react_vol_integ_local);
      gkyl_free(sca->react_vol_integ_global);
      gkyl_free(sca->bflux_m0_vol_integ_local);
    }
    gkyl_array_release(sca->reactivity);
    gkyl_array_release(sca->dndt_react);

    gkyl_dg_iz_release(sca->iz_react_calc);
  } 
}
