#include <assert.h>
#include <gkyl_vlasov_priv.h>

// Evaluate a p=1 configuration-space expansion at the computational coordinate xn (the
// leading cdim entries): the cell containing xn and the reference coordinate in it.
static double
eval_conf_field(
  const struct vm_pkpm_source_eval_ctx *c, const struct gkyl_array *arr, int comp, const double *xn
)
{
  int idx[GKYL_MAX_CDIM];
  double xc[GKYL_MAX_CDIM], eta[GKYL_MAX_CDIM];
  gkyl_rect_grid_coord_idx(c->conf_grid, xn, idx);
  for (int d = 0; d < c->cdim; ++d) {
    idx[d] = GKYL_MAX2(GKYL_MIN2(idx[d], c->local_ext->upper[d]), c->local_ext->lower[d]);
  }
  gkyl_rect_grid_cell_center(c->conf_grid, idx, xc);
  for (int d = 0; d < c->cdim; ++d) {
    eta[d] = 2.0 * (xn[d] - xc[d]) / c->conf_grid->dx[d];
  }
  const double *coeff = gkyl_array_cfetch(arr, gkyl_range_idx(c->local_ext, idx));
  return c->conf_basis->eval_expand(eta, &coeff[comp * c->conf_basis->num_basis]);
}

// The source distributions in the normalized variables of the species at a phase-space
// quadrature point (computational coordinates: the profiles were projected on the
// computational grid at the physical coordinates of the position map). With f_s the
// parallel Maxwellian of the source, n_s M(v_par; u_s.b, T_s/m), integrated over v_perp,
// the normalized F = int f dv_perp / n in w = (v_par - u.b)/xi receives
//   S_F = (xi/n) f_s(u.b + xi w),       S_G = (T_s/(m xi^2)) S_F
// (G is the perpendicular energy per xi^2: the source carries T_s per particle; a
// perpendicular drift of the source relative to the flow is neglected), and the
// normalization change -(m n_s/rho) [F, G] is added separately as a product with f.
static void
eval_source_dist(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct vm_pkpm_source_eval_ctx *c = ctx;
  int cdim = c->cdim;
  double w = xn[cdim];
  double n_s = eval_conf_field(c, c->prof, 0, xn), T_s = eval_conf_field(c, c->prof, 1, xn);
  double us[3], u[3], b[3];
  for (int i = 0; i < 3; ++i) {
    us[i] = eval_conf_field(c, c->prof, 2 + i, xn);
    u[i] = eval_conf_field(c, c->u, i, xn);
    b[i] = eval_conf_field(c, c->bvar, i, xn);
  }
  double rho = eval_conf_field(c, c->rho, 0, xn), xi = eval_conf_field(c, c->xi, 0, xn);
  double n = rho / c->mass;
  if (!(n > 0.0) || !(xi > 0.0) || !(T_s > 0.0) || n_s <= 0.0) {
    fout[0] = fout[1] = 0.0;
    return;
  }
  double u_par = b[0] * u[0] + b[1] * u[1] + b[2] * u[2];
  double us_par = b[0] * us[0] + b[1] * us[1] + b[2] * us[2];
  double vtsq = T_s / c->mass;
  double v = u_par + xi * w - us_par;
  double f_s = n_s * exp(-0.5 * v * v / vtsq) / sqrt(2.0 * M_PI * vtsq);
  fout[0] = (xi / n) * f_s;
  fout[1] = (T_s / (c->mass * xi * xi)) * fout[0];
}

// The fluid image of the source profile, [rho, M, E] = [m n_s, m n_s u_s, m n_s |u_s|^2/2 + 3 n_s T_s/2],
// and the profile itself [n_s, T_s, u_s] from the user functions.
static void
eval_source_profile(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct vm_pkpm_source_eval_ctx *c = ctx;
  double n_s = 0.0, T_s = 0.0, us[3] = {0.0};
  c->info->density(t, xn, &n_s, c->info->ctx_density);
  c->info->temp(t, xn, &T_s, c->info->ctx_temp);
  if (c->info->V_drift) {
    c->info->V_drift(t, xn, us, c->info->ctx_V_drift);
  }
  fout[0] = n_s;
  fout[1] = T_s;
  for (int i = 0; i < 3; ++i) {
    fout[2 + i] = us[i];
  }
  double m = c->mass;
  fout[5] = m * n_s;
  for (int i = 0; i < 3; ++i) {
    fout[6 + i] = m * n_s * us[i];
  }
  fout[9] = 0.5 * m * n_s * (us[0] * us[0] + us[1] * us[1] + us[2] * us[2]) + 1.5 * n_s * T_s;
}

void
vm_pkpm_species_source_init(
  gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_source *src
)
{
  int nc = app->basis.num_basis, np = s->kinetic->basis.num_basis;
  const struct gkyl_vlasov_source *info = &s->kinetic->info.source;
  src->source_id = info->source_id;
  src->calc_bflux = false;
  src->scale = 1.0;
  if (src->source_id == GKYL_NO_SOURCE) {
    return;
  }
  // One LTE profile, added as given or scaled by the wall losses; time-dependent and adaptive
  // sources are not available.
  assert(src->source_id == GKYL_PROJ_SOURCE || src->source_id == GKYL_BFLUX_SOURCE);
  assert(info->num_sources == 1 && info->projection[0].proj_id == GKYL_PROJ_VLASOV_LTE);
  assert(info->projection[0].density && info->projection[0].temp);
  assert(!info->evolve_source);
  if (src->source_id == GKYL_BFLUX_SOURCE) {
    assert(info->source_length > 0.0);
    src->calc_bflux = true;
    src->source_length = info->source_length;
    src->source_species = vlasov_find_species(app, info->source_species);
    assert(src->source_species);
    src->scale = 0.0;
  }
  src->eval_ctx = (struct vm_pkpm_source_eval_ctx){
    .cdim = app->cdim,
    .mass = s->mass,
    .conf_grid = &app->grid,
    .conf_basis = &app->basis,
    .local_ext = &app->local_ext,
    .info = &info->projection[0],
  };

  // Profile [n_s, T_s, u_s] and its fluid image, sampled at the physical coordinates.
  struct gkyl_array *prof10 = mkarr(false, 10 * nc, app->local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &app->grid,
    .basis = &app->basis,
    .num_quad = app->poly_order + 1,
    .num_ret_vals = 10,
    .eval = eval_source_profile,
    .ctx = &src->eval_ctx,
    .c2p_func = vm_pos_c2p,
    .c2p_func_ctx = app->pos_map,
  });
  gkyl_proj_on_basis_advance(proj, 0.0, &app->local_ext, prof10);
  gkyl_proj_on_basis_release(proj);
  src->prof_host = mkarr(false, 5 * nc, app->local_ext.volume);
  src->fluid_src_host = mkarr(false, 5 * nc, app->local_ext.volume);
  gkyl_array_set_offset(src->prof_host, 1.0, prof10, 0);
  gkyl_array_set_offset(src->fluid_src_host, 1.0, prof10, 5 * nc);
  gkyl_array_release(prof10);
  // The evolved fluid carries the position-map Jacobian.
  src->fluid_src = mkarr(app->use_gpu, 5 * nc, app->local_ext.volume);
  gkyl_array_copy(src->fluid_src, src->fluid_src_host);
  gkyl_vlasov_position_map_rescale_jacobpos_conf(
    app->pos_map, &app->local_ext, src->fluid_src, src->fluid_src
  );
  // The fluid aspect adds the fluid image of the profile, scaled by the rate of the stage.
  s->fluid->pkpm.src_fluid = src->fluid_src;
  s->fluid->pkpm.src_scale = src->scale;

  // Kinetic source of the stage: projected on the host from the species' nodal state.
  src->f_src = mkarr(app->use_gpu, 2 * np, s->kinetic->local_ext.volume);
  src->f_src_host = app->use_gpu ? mkarr(false, 2 * np, s->kinetic->local_ext.volume) :
                                   gkyl_array_acquire(src->f_src);
  gkyl_array_clear(src->f_src, 0.0);
  src->n_ratio = mkarr(app->use_gpu, nc, app->local_ext.volume);
  src->n_src = mkarr(app->use_gpu, nc, app->local_ext.volume);
  struct gkyl_array *n_src_host = mkarr(false, nc, app->local_ext.volume);
  gkyl_array_set_offset(n_src_host, s->mass, src->prof_host, 0); // m n_s (physical)
  gkyl_array_copy(src->n_src, n_src_host);
  gkyl_array_release(n_src_host);
  src->rho_host = app->use_gpu ? mkarr(false, nc, app->local_ext.volume) : 0;
  src->u_host = app->use_gpu ? mkarr(false, 3 * nc, app->local_ext.volume) : 0;
  src->xi_host = app->use_gpu ? mkarr(false, 3 * nc, app->local_ext.volume) : 0;
  src->bvar_host = app->use_gpu ? mkarr(false, 3 * nc, app->local_ext.volume) : 0;
  src->div_mem = app->use_gpu ? gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, nc) :
                                gkyl_dg_bin_op_mem_new(app->local.volume, nc);
  src->proj_f = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &s->kinetic->grid,
    .basis = &s->kinetic->basis,
    .num_quad = 4,
    .num_ret_vals = 2,
    .eval = eval_source_dist,
    .ctx = &src->eval_ctx,
  });
}

// Particle loss rate of the species through its absorbing walls from the stored
// absorbing-wall fluid fluxes of the stage (the mass flux F_rho, integrated over the
// faces; the lower wall's flux is in the skin cell's slot and negative when outgoing, the
// upper wall's in the ghost cell's slot and positive). Summed over ranks.
void
vm_pkpm_species_calc_bflux(gkyl_vlasov_app *app, struct vlasov_species *s)
{
  int cdim = app->cdim, nk = app->basis.num_basis / 2;
  double loss[2 * GKYL_MAX_CDIM] = {0.0};
  int is_np_bc[3] = {1, 1, 1};
  for (int d = 0; d < app->num_periodic_dir; ++d) {
    is_np_bc[app->periodic_dirs[d]] = 0;
  }
  const struct gkyl_array *flux = s->fluid->pkpm.fluid_flux_surf;
  if (app->use_gpu) {
    gkyl_array_copy(s->fluid->pkpm.fluid_flux_surf_host, s->fluid->pkpm.fluid_flux_surf);
    flux = s->fluid->pkpm.fluid_flux_surf_host;
  }
  for (int d = 0; d < cdim; ++d) {
    if (!is_np_bc[d]) {
      continue;
    }
    double face_fac = 1.0;
    for (int dd = 0; dd < cdim; ++dd) {
      if (dd != d) {
        face_fac *= app->grid.dx[dd] / sqrt(2.0);
      }
    }
    bool on_lower = app->local.lower[d] == app->global.lower[d];
    bool on_upper = app->local.upper[d] == app->global.upper[d];
    struct gkyl_range_iter iter;
    if (on_lower && s->kinetic->lower_bc[d].type == GKYL_SPECIES_ABSORB) {
      gkyl_range_iter_init(&iter, &app->lower_skin[d]);
      while (gkyl_range_iter_next(&iter)) {
        const double *fl = gkyl_array_cfetch(flux, gkyl_range_idx(&app->local_ext, iter.idx));
        loss[2 * d] -= fl[d * 5 * nk] * face_fac;
      }
    }
    if (on_upper && s->kinetic->upper_bc[d].type == GKYL_SPECIES_ABSORB) {
      gkyl_range_iter_init(&iter, &app->upper_ghost[d]);
      while (gkyl_range_iter_next(&iter)) {
        const double *fl = gkyl_array_cfetch(flux, gkyl_range_idx(&app->local_ext, iter.idx));
        loss[2 * d + 1] += fl[d * 5 * nk] * face_fac;
      }
    }
  }
  for (int i = 0; i < 2 * cdim; ++i) {
    loss[i] /= s->mass;
  }
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 2 * cdim, loss, s->kinetic->pkpm.bflux);
}

void
vm_pkpm_species_source_calc_rate(
  gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_source *src
)
{
  if (!src->calc_bflux) {
    return;
  }
  src->scale = vlasov_species_wall_loss_rate(app, src->source_species) / src->source_length;
  s->fluid->pkpm.src_scale = src->scale; // the fluid aspect scales the fluid image of the profile
}

// Source distributions of the stage [J S_F, J S_G] and the density ratio m n_s/rho, from the
// nodal state of the stage.
void
vm_pkpm_species_source_calc_dist(
  gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_source *src
)
{
  if (src->source_id == GKYL_NO_SOURCE || src->scale == 0.0) {
    return;
  }
  struct timespec tm = gkyl_wall_clock();
  // Host copies of the nodal state for the projection callback.
  if (app->use_gpu) {
    gkyl_array_copy(src->rho_host, s->fluid->pkpm.rho);
    gkyl_array_copy(src->u_host, s->fluid->pkpm.u);
    gkyl_array_copy(src->xi_host, s->fluid->pkpm.xi);
    gkyl_array_copy(src->bvar_host, app->field->bvar);
    src->eval_ctx.rho = src->rho_host;
    src->eval_ctx.u = src->u_host;
    src->eval_ctx.xi = src->xi_host;
    src->eval_ctx.bvar = src->bvar_host;
  } else {
    src->eval_ctx.rho = s->fluid->pkpm.rho;
    src->eval_ctx.u = s->fluid->pkpm.u;
    src->eval_ctx.xi = s->fluid->pkpm.xi;
    src->eval_ctx.bvar = app->field->bvar;
  }
  src->eval_ctx.prof = src->prof_host;
  gkyl_proj_on_basis_advance(src->proj_f, app->tcurr, &s->kinetic->local, src->f_src_host);
  if (app->use_gpu) {
    gkyl_array_copy(src->f_src, src->f_src_host);
  }
  gkyl_vlasov_position_map_rescale_jacobpos(
    app->pos_map, &s->kinetic->basis, &s->kinetic->local, src->f_src, src->f_src
  );
  gkyl_dg_div_op_range(
    src->div_mem, &app->basis, 0, src->n_ratio, 0, src->n_src, 0, s->fluid->pkpm.rho, &app->local
  );
  app->stat.species_rhs_tm += gkyl_time_diff_now_sec(tm);
}

// Kinetic source of the stage into the configuration-pass right-hand side (before the
// closures, which keep the constraints exact): scale [J S_F, J S_G] - scale (m n_s/rho) f.
void
vm_pkpm_species_source_kinetic_rhs(
  gkyl_vlasov_app *app, struct vm_species *kin, struct vm_pkpm_source *src,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  if (src->source_id == GKYL_NO_SOURCE || src->scale == 0.0) {
    return;
  }
  struct timespec tm = gkyl_wall_clock();
  // the source distribution is gyrotropic: it adds to the [F, G] blocks only
  gkyl_array_accumulate_offset_range(
    rhs, src->scale, src->f_src, PKPM_DIST_F * kin->basis.num_basis, &kin->local
  );
  // Normalization change: the density source removes (m n_s / rho) of every normalized
  // distribution [F, G, M0_perp, N_perp].
  vm_pkpm_species_accumulate_conf_phase(app, kin, rhs, -src->scale, src->n_ratio, fin);
  app->stat.species_rhs_tm += gkyl_time_diff_now_sec(tm);
}

void
vm_pkpm_species_source_write(
  gkyl_vlasov_app *app, struct vlasov_species *s, struct vm_pkpm_source *src, double tm, int frame
)
{
  if (src->source_id == GKYL_NO_SOURCE || !s->kinetic->info.source.write_source) {
    return;
  }
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new((struct vlasov_output_meta){
    .frame = frame,
    .stime = tm,
    .poly_order = app->poly_order,
    .basis_type = app->basis.id,
  });
  const char *fmt = "%s-%s_source_fluid_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, s->name, frame);
  char fileNm[sz + 1];
  snprintf(fileNm, sizeof fileNm, fmt, app->name, s->name, frame);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, src->fluid_src_host, fileNm);
  vlasov_array_meta_release(mt);
}

void
vm_pkpm_species_source_release(const gkyl_vlasov_app *app, const struct vm_pkpm_source *src)
{
  if (src->source_id == GKYL_NO_SOURCE) {
    return;
  }
  gkyl_proj_on_basis_release(src->proj_f);
  gkyl_dg_bin_op_mem_release(src->div_mem);
  gkyl_array_release(src->prof_host);
  gkyl_array_release(src->fluid_src_host);
  gkyl_array_release(src->fluid_src);
  gkyl_array_release(src->f_src);
  gkyl_array_release(src->f_src_host);
  gkyl_array_release(src->n_ratio);
  gkyl_array_release(src->n_src);
  if (app->use_gpu) {
    gkyl_array_release(src->rho_host);
    gkyl_array_release(src->u_host);
    gkyl_array_release(src->xi_host);
    gkyl_array_release(src->bvar_host);
  }
}
