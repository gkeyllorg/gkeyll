#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gkns_react_cross_moms_disabled(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_react *react, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  // Do nothing.
}

static void
gkns_react_cross_moms_enabled(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_react *react, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  struct timespec wst = gkyl_wall_clock();    
  for (int i=0; i<react->num_react; ++i) {
    struct gk_species *gks_elc = &app->species[react->elc_idx[i]]; 
    struct gk_species *gks_ion = &app->species[react->ion_idx[i]]; 

    if (react->react_id[i] == GKYL_REACT_IZ) {
      // compute needed electron Maxwellian moments (J*n, u_par, T/m) 
      gk_species_moment_calc(&gks_elc->lte.moms, 
        gks_elc->local, app->local, fin[react->elc_idx[i]]);

      // divide out the Jacobian from the electron density for computing reaction rates
      gkyl_dg_div_op_range(gks_elc->lte.moms.mem_geo, &app->basis, 
        0, gks_elc->lte.moms.marr, 0, gks_elc->lte.moms.marr, 0, 
        app->gk_geom->geo_int.jacobgeo, &app->local); 

      // compute ionization reaction rate from input electron primitive moments
      gkyl_dg_iz_coll(react->iz[i], gks_elc->lte.moms.marr, 
        react->vt_sq_iz1[i], react->vt_sq_iz2[i], react->coeff_react[i], 0);
    }
    else if (react->react_id[i] == GKYL_REACT_RECOMB) {
      // compute needed electron Maxwellian moments (J*n, u_par, T/m) 
      gk_species_moment_calc(&gks_elc->lte.moms, 
        gks_elc->local, app->local, fin[react->elc_idx[i]]);

      // Copy J*n for use in final update
      gkyl_array_set_range(react->Jm0_elc[i], 1.0, gks_elc->lte.moms.marr, &app->local);

      // divide out the Jacobian from the electron density for computing reaction rates
      gkyl_dg_div_op_range(gks_elc->lte.moms.mem_geo, &app->basis, 
        0, gks_elc->lte.moms.marr, 0, gks_elc->lte.moms.marr, 0, 
        app->gk_geom->geo_int.jacobgeo, &app->local); 

      // compute needed ion Maxwellian moments (J*n, u_par, T/m) 
      gk_species_moment_calc(&gks_ion->lte.moms, 
        gks_ion->local, app->local, fin[react->ion_idx[i]]);

      // divide out the Jacobian from the ion density for use in Maxwellian projection
      gkyl_dg_div_op_range(gks_ion->lte.moms.mem_geo, &app->basis,
        0, gks_ion->lte.moms.marr, 0, gks_ion->lte.moms.marr, 0,
        app->gk_geom->geo_int.jacobgeo, &app->local);

      // Construct ion vector velocity upar b_i
      // if cdim = 1, u0 = upar, if cdim = 2, u1 = upar, if cdim = 3, u2 = upar
      gkyl_array_clear(react->upar_ion[i], 0.0);
      gkyl_array_set_offset(react->u_i_dot_b_i[i], 1.0, gks_ion->lte.moms.marr, 1*app->basis.num_basis);
      gkyl_array_set_offset(react->upar_ion[i], 1.0, react->u_i_dot_b_i[i], (app->cdim-1)*app->basis.num_basis);
      gkyl_array_clear(react->u_i_dot_b_i[i], 0.0);

      // Copy vt^2 = T/m of the ions
      gkyl_array_set_offset(react->vt_sq_ion[i], 1.0, gks_ion->lte.moms.marr, 2*app->basis.num_basis);
      
      // compute recombination reaction rate
      gkyl_dg_recomb_coll(react->recomb[i], gks_elc->lte.moms.marr, 
        react->coeff_react[i], 0);
    }
    else if (react->react_id[i] == GKYL_REACT_CX) {
      // compute needed ion Maxwellian moments (J*n, u_par, T/m) 
      gk_species_moment_calc(&gks_ion->lte.moms, 
        gks_ion->local, app->local, fin[react->ion_idx[i]]);

      // divide out the Jacobian from the ion density
      gkyl_dg_div_op_range(gks_ion->lte.moms.mem_geo, &app->basis,
        0, gks_ion->lte.moms.marr, 0, gks_ion->lte.moms.marr, 0,
        app->gk_geom->geo_int.jacobgeo, &app->local); 

      // Construct ion vector velocity upar b_i
      // if cdim = 1, u0 = upar, if cdim = 2, u1 = upar, if cdim = 3, u2 = upar
      gkyl_array_clear(react->upar_ion[i], 0.0);
      gkyl_array_set_offset(react->u_i_dot_b_i[i], 1.0, gks_ion->lte.moms.marr, 1*app->basis.num_basis);
      gkyl_array_set_offset(react->upar_ion[i], 1.0, react->u_i_dot_b_i[i], (app->cdim-1)*app->basis.num_basis);
      gkyl_array_clear(react->u_i_dot_b_i[i], 0.0);

      // Copy vt^2 = T/m of the ions (partner of the neutrals)
      gkyl_array_set_offset(react->vt_sq_ion[i], 1.0, gks_ion->lte.moms.marr, 2*app->basis.num_basis);

      // compute needed partner (neutral) Maxwellian moments (J*n, ux, uy, uz, T/m) 
      struct gk_neut_species *gkns_partner = &app->neut_species[react->partner_idx[i]];
      gk_neut_species_moment_calc(&gkns_partner->lte.moms, 
        gkns_partner->local, app->local, fin_neut[react->partner_idx[i]]);  

      // Copy J*n for use in final update
      gkyl_array_set_range(react->Jm0_partner[i], 1.0, gkns_partner->lte.moms.marr, &app->local);

      // Divide out the Jacobian from the partner density.
      gkyl_dg_div_op_range(gkns_partner->lte.moms.mem_geo, &app->basis,
        0, gkns_partner->lte.moms.marr, 0, gkns_partner->lte.moms.marr, 0,
        app->gk_geom->geo_int.jacobgeo, &app->local); 

      // prim_vars_neut_gk is returned to prim_vars[i] here.
      gkyl_dg_cx_coll(react->cx[i], gks_ion->lte.moms.marr, gkns_partner->lte.moms.marr, 
        react->upar_ion[i], react->coeff_react[i], 0);
    }
  }
  app->stat.neut_species_react_mom_tm += gkyl_time_diff_now_sec(wst);
}

static void
gkns_react_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_react *react, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  // Do nothing.
}

static void
gkns_react_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_react *react, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();  
  for (int i=0; i<react->num_react; ++i) {
    gkyl_array_clear(react->f_react, 0.0);
    gkyl_array_clear(react->react_lte_moms[i], 0.0);

    struct gk_species *gks_elc = &app->species[react->elc_idx[i]];
    struct gk_species *gks_ion = &app->species[react->ion_idx[i]];

    // neutral species can only be the donor or the receiver in reaction
    if (react->react_id[i] == GKYL_REACT_IZ) {
      // donor update is -n_elc*coeff_react*f_donor

      // Accumulate -n_elc*(J*f_donor) (*note* Jacobian factor already included in fin)
      gkyl_dg_mul_conf_phase_op_accumulate_range(&app->basis, &s->basis, react->f_react,
        -1.0, gks_elc->lte.moms.marr, fin, &app->local, &s->local);  
    }
    else if (react->react_id[i] == GKYL_REACT_RECOMB) {
      // Copy components of partner ions into reaction moments 
      gkyl_array_set_offset(react->react_lte_moms[i], 1.0, gks_ion->lte.moms.marr, 0*app->basis.num_basis); 
      // Overwrite flow velocity and vt^2 to be upar b_i (vector) and vt^2 of the ions
      gkyl_array_set_offset(react->react_lte_moms[i], 1.0, react->upar_ion[i], 1*app->basis.num_basis); 
      gkyl_array_set_offset(react->react_lte_moms[i], 1.0, react->vt_sq_ion[i], 4*app->basis.num_basis); 
      gk_neut_species_lte_from_moms(app, s, &s->lte, react->react_lte_moms[i]);

      // Accumulate J*n_elc*fmax(n_ion, upar bx, upar by, upar bz, vt_ion^2) onto f_react
      gkyl_dg_mul_conf_phase_op_accumulate_range(&app->basis, &s->basis, react->f_react,
        1.0, react->Jm0_elc[i], s->lte.f_lte, &app->local, &s->local);
    }
    else if (react->react_id[i] == GKYL_REACT_CX) {
      // neutral update is coeff_react*(n_partner*fmax(n_ion, upar bx, upar by, upar bz, vt_ion^2) - n_ion*f_partner)
      // Copy components of partner ions into reaction moments 
      gkyl_array_set_offset(react->react_lte_moms[i], 1.0, gks_ion->lte.moms.marr, 0*app->basis.num_basis); 
      // Overwrite flow velocity and vt^2 to be upar b_i (vector) and vt^2 of the ions
      gkyl_array_set_offset(react->react_lte_moms[i], 1.0, react->upar_ion[i], 1*app->basis.num_basis); 
      gkyl_array_set_offset(react->react_lte_moms[i], 1.0, react->vt_sq_ion[i], 4*app->basis.num_basis); 
      gk_neut_species_lte_from_moms(app, s, &s->lte, react->react_lte_moms[i]);

      // Accumulate J*n_partner*fmax(n_ion, upar bx, upar by, upar bz, vt_ion^2) onto f_react
      gkyl_dg_mul_conf_phase_op_accumulate_range(&app->basis, &s->basis, react->f_react,
        1.0, react->Jm0_partner[i], s->lte.f_lte, &app->local, &s->local);

      // Accumulate -n_ion*(J*f_partner) (*note* Jacobian factor already included in fin)
      gkyl_dg_mul_conf_phase_op_accumulate_range(&app->basis, &s->basis, react->f_react,
        -1.0, gks_ion->lte.moms.marr, fin, &app->local, &s->local); 
    }

    // Accumulate reaction update to rhs 
    gkyl_dg_mul_conf_phase_op_accumulate_range(&app->basis, &s->basis, rhs,
      1.0, react->coeff_react[i], react->f_react, &app->local, &s->local);  
  }
  app->stat.neut_species_react_tm += gkyl_time_diff_now_sec(wst);
}

static void
gkns_react_write_disabled(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, struct gk_react *gkr,
  int ridx, double tm, int frame)
{
  // Do nothing.
}

void
gkns_react_write_enabled(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, struct gk_react *gkr,
  int ridx, double tm, int frame)
{
  // React diagnostics usually written from gk_species.
  // In the case of static gk_species, write_diagnostics flag
  // can be used to check reaction rates from gk_neut_species.
  struct timespec wst = gkyl_wall_clock();
  // Compute reaction rate
  const struct gkyl_array *fin[app->num_species];
  const struct gkyl_array *fin_neut[app->num_neut_species];
  for (int i=0; i<app->num_species; ++i) {
    fin[i] = app->species[i].f;
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    fin_neut[i] = app->neut_species[i].f;
  }
  gk_neut_species_react_cross_moms(app, gkns, gkr, fin, fin_neut);
  app->stat.neut_species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
  
  struct timespec wtm = gkyl_wall_clock();

  // Package metadata.
  gkyl_msgpack_map_elem_set_double(gkns->io_meta_conf_len, gkns->io_meta_conf, "time", tm);
  gkyl_msgpack_map_elem_set_uint(gkns->io_meta_conf_len, gkns->io_meta_conf, "frame", frame);
  int io_meta_len[] = {gkns->io_meta_conf_len, app->gk_geom->io_meta_basic_len, 1};
  struct gkyl_msgpack_data *mt;

  if (app->use_gpu)
    gkyl_array_copy(gkr->coeff_react_host[ridx], gkr->coeff_react[ridx]);

  if (gkr->react_id[ridx] == GKYL_REACT_IZ) {
    const char *fmt = "%s-%s_%s_%s_iz_react_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gkns->info.name,
      gkr->react_type[ridx].elc_nm, gkr->react_type[ridx].ion_nm, frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gkns->info.name,
     gkr->react_type[ridx].elc_nm, gkr->react_type[ridx].ion_nm, frame);

    struct gkyl_msgpack_map_elem desc[] = {
      { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Ionization reaction rate." }
    };
    const struct gkyl_msgpack_map_elem* io_meta[] = {gkns->io_meta_conf, app->gk_geom->io_meta_basic, desc};
    mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);
  
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, 
      gkr->coeff_react_host[ridx], fileNm);
  }
  if (gkr->react_id[ridx] == GKYL_REACT_RECOMB) {
    const char *fmt = "%s-%s_%s_%s_recomb_react_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gkns->info.name,
      gkr->react_type[ridx].elc_nm, gkr->react_type[ridx].ion_nm, frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gkns->info.name,
      gkr->react_type[ridx].elc_nm, gkr->react_type[ridx].ion_nm, frame);
  
    struct gkyl_msgpack_map_elem desc[] = {
      { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Recombination reaction rate." }
    };
    const struct gkyl_msgpack_map_elem* io_meta[] = {gkns->io_meta_conf, app->gk_geom->io_meta_basic, desc};
    mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);
  
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, 
      gkr->coeff_react_host[ridx], fileNm);
  }
  if (gkr->react_id[ridx] == GKYL_REACT_CX) {
    const char *fmt = "%s-%s_%s_cx_react_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gkns->info.name,
      gkr->react_type[ridx].ion_nm, frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gkns->info.name,
      gkr->react_type[ridx].ion_nm, frame);
    
    struct gkyl_msgpack_map_elem desc[] = {
      { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Charge exchange reaction rate." }
    };
    const struct gkyl_msgpack_map_elem* io_meta[] = {gkns->io_meta_conf, app->gk_geom->io_meta_basic, desc};
    mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);
  
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, 
      gkr->coeff_react_host[ridx], fileNm);
  }
  app->stat.n_neut_diag_io += 1;
  
  gkyl_msgpack_data_release(mt); 
  app->stat.neut_species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
}

void 
gk_neut_species_react_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, 
  struct gkyl_gyrokinetic_react inp, struct gk_react *react)
{
  react->num_react = inp.num_react; 
  react->write_diagnostics = inp.write_diagnostics;

  react->cross_moms_func_neut = gkns_react_cross_moms_disabled;
  react->rhs_func_neut = gkns_react_rhs_disabled;
  react->write_func_neut = gkns_react_write_disabled;

  if (react->num_react) {
    // Initialize information about reactions from input struct.
    for (int i=0; i<react->num_react; ++i) 
      react->react_type[i] = inp.react_type[i];

    // Methods chosen at runtime:
    react->cross_moms_func_neut = gkns_react_cross_moms_enabled;
    react->rhs_func_neut = gkns_react_rhs_enabled;
    if (react->write_diagnostics)
      react->write_func_neut = gkns_react_write_enabled;
  }
}

static double
gk_species_react_get_vt_sq_min(struct gkyl_gyrokinetic_app *app, struct gk_species *s)
{
  double bmag_mid = app->bmag_ref;

  int vdim = s->info.vdim;
  double dv_min[vdim];
  gkyl_velocity_map_reduce_dv_range(s->vel_map, GKYL_MIN, dv_min, s->vel_map->local_vel);

  double tpar_min = (s->info.mass/6.0)*pow(dv_min[0],2);
  double tperp_min = vdim>1 ? (bmag_mid/3.0)*dv_min[1] : tpar_min;
  return (tpar_min + 2.0*tperp_min)/(3.0*s->info.mass);
}

static double
gk_neut_species_react_get_vt_sq_min(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s)
{
  double bmag_mid = app->bmag_ref;

  int vdim = s->info.vdim;
  double dv_min[vdim];
  gkyl_velocity_map_reduce_dv_range(s->vel_map, GKYL_MIN, dv_min, s->vel_map->local_vel);

  double t_min = 0.0;
  for (int i=0; i<vdim; i++)
    t_min += (s->info.mass/6.0)*pow(dv_min[0],2);
  return t_min/(3.0*s->info.mass);
}

void 
gk_neut_species_react_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, struct gk_react *react)
{
  if (react->num_react) {
    // Distribution function which holds update for each reaction
    // form depend on react->type_self, e.g., for recombination and react->type_self == GKYL_SELF_RECVR
    // react->f_react = n_elc*coeff_react*fmax(n_ion, upar_ion b_i, vt_ion^2)
    react->f_react = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);

    for (int i=0; i<react->num_react; ++i) {
      react->react_id[i] = react->react_type[i].react_id;
      react->type_self[i] = react->react_type[i].type_self;
      // Fetch pointers to species objects
      react->species_elc[i] = gk_find_species(app, react->react_type[i].elc_nm);
      react->species_ion[i] = gk_find_species(app, react->react_type[i].ion_nm);
      // Fetch index of species for indexing arrays
      react->elc_idx[i] = gk_find_species_idx(app, react->react_type[i].elc_nm);
      react->ion_idx[i] = gk_find_species_idx(app, react->react_type[i].ion_nm);

      // Compute a minimum representable temperature based on the smallest dv in the grid.
      double ion_vt_sq_min = gk_species_react_get_vt_sq_min(app,  &app->species[react->ion_idx[i]]);
      double neut_vt_sq_min; 
      if (gk_find_neut_species(app, react->react_type[i].donor_nm)) {
        react->donor_idx[i] = gk_find_neut_species_idx(app, react->react_type[i].donor_nm);
      }
      else if (gk_find_neut_species(app, react->react_type[i].partner_nm)) {
        react->partner_idx[i] = gk_find_neut_species_idx(app, react->react_type[i].partner_nm);
        neut_vt_sq_min = gk_neut_species_react_get_vt_sq_min(app, &app->neut_species[react->partner_idx[i]]);
      }

      react->coeff_react[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

      // Reaction LTE moments needed for projecting LTE distribution functions
      react->react_lte_moms[i] = mkarr(app->use_gpu, 5*app->basis.num_basis, app->local_ext.volume);

      react->vt_sq_iz1[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      react->vt_sq_iz2[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

      react->u_i_dot_b_i[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

      // J*n for use in final update formulae of reactions
      react->Jm0_elc[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      react->Jm0_partner[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

      // Ion flow velocity (upar b_i), ion vt^2 for projecting LTE distribution functions
      react->upar_ion[i] = mkarr(app->use_gpu, 3*app->basis.num_basis, app->local_ext.volume);
      react->vt_sq_ion[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

      if (react->react_id[i] == GKYL_REACT_IZ) {
        struct gkyl_dg_iz_inp iz_inp = {
          .cbasis = &app->basis, 
          .conf_rng = &app->local, 
          .type_ion = react->react_type[i].ion_id, 
          .charge_state = react->react_type[i].charge_state, 
          .type_self = react->type_self[i], 
        };
        react->iz[i] = gkyl_dg_iz_new(&iz_inp, app->use_gpu); 
      }
      else if (react->react_id[i] == GKYL_REACT_RECOMB) {
        struct gkyl_dg_recomb_inp recomb_inp = {
          .grid = &s->grid,
          .cbasis = &app->basis, 
          .pbasis = &s->basis, 
          .conf_rng = &app->local, 
          .conf_rng_ext = &app->local_ext,
          .phase_rng = &s->local, 
          .mass_self = s->info.mass, 
          .type_ion = react->react_type[i].ion_id, 
          .charge_state = react->react_type[i].charge_state, 
          .type_self = react->type_self[i], 
        };
        react->recomb[i] = gkyl_dg_recomb_new(&recomb_inp, app->use_gpu);       
      }
      else if (react->react_id[i] == GKYL_REACT_CX) {
        struct gk_species *gks = &app->species[react->ion_idx[i]];
        struct gkyl_dg_cx_inp cx_inp = {
          .cbasis = &app->basis,
          .conf_rng = &app->local,
          .vt_sq_ion_min = ion_vt_sq_min, 
          .vt_sq_neut_min = neut_vt_sq_min, 
          .type_ion = react->react_type[i].ion_id,
        };
        react->cx[i] = gkyl_dg_cx_new(&cx_inp, app->use_gpu);
      }

      if (react->write_diagnostics) {
        react->coeff_react_host[i] = app->use_gpu? mkarr(false, react->coeff_react[i]->ncomp, react->coeff_react[i]->size)
                                                 : gkyl_array_acquire(react->coeff_react[i]);
      }
    }
  }
}

void
gk_neut_species_react_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_react *react, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  react->cross_moms_func_neut(app, species, react, fin, fin_neut);
}

void
gk_neut_species_react_rhs(gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_react *react, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  react->rhs_func_neut(app, s, react, fin, rhs);
}

void
gk_neut_species_react_write(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, struct gk_react *gkr,
  int ridx, double tm, int frame)
{
  gkr->write_func_neut(app, gkns, gkr, ridx, tm, frame);
}

void 
gk_neut_species_react_release(const struct gkyl_gyrokinetic_app *app, const struct gk_react *react)
{
  if (react->num_react) {
    gkyl_array_release(react->f_react);

    for (int i=0; i<react->num_react; ++i) {
      gkyl_array_release(react->coeff_react[i]);
      if (react->write_diagnostics) {
        gkyl_array_release(react->coeff_react_host[i]);
      }
      gkyl_array_release(react->react_lte_moms[i]);
      gkyl_array_release(react->vt_sq_iz1[i]);
      gkyl_array_release(react->vt_sq_iz2[i]);
      gkyl_array_release(react->Jm0_elc[i]);
      gkyl_array_release(react->Jm0_partner[i]);
      gkyl_array_release(react->upar_ion[i]); 
      gkyl_array_release(react->vt_sq_ion[i]);
      gkyl_array_release(react->u_i_dot_b_i[i]);

      if (react->react_id[i] == GKYL_REACT_IZ) {
        gkyl_dg_iz_release(react->iz[i]);
      }
      else if (react->react_id[i] == GKYL_REACT_RECOMB) {
        gkyl_dg_recomb_release(react->recomb[i]);
      }
      else if (react->react_id[i] == GKYL_REACT_CX) {
        gkyl_dg_cx_release(react->cx[i]);
      }
    }
  }
}
