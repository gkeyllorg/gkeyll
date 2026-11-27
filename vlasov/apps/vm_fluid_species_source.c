#include <assert.h>
#include <gkyl_vlasov_priv.h>
#include <gkyl_wv_euler.h>

// Context of the projection of the fluid image of the LTE source profiles.
struct vm_fluid_source_eval_ctx {
  const struct gkyl_vlasov_fluid_source *info; // source input
  double mass; // mass of the species
  double inv_gamma_m1; // 1/(gas_gamma - 1) (Euler)
  int num_equations; // components of the state: 5 (Euler) or 4 (isothermal Euler)
};

// The fluid image of the sum of the LTE profiles [n_s, T_s, u_s]:
// [m n_s, m n_s u_s, m n_s |u_s|^2/2 + n_s T_s/(gas_gamma - 1)] (no energy row for isothermal Euler).
static void
eval_source_image(double t, const double *xn, double *restrict fout, void *ctx)
{
  const struct vm_fluid_source_eval_ctx *c = ctx;
  for (int k = 0; k < c->num_equations; ++k) {
    fout[k] = 0.0;
  }
  for (int i = 0; i < c->info->num_sources; ++i) {
    const struct gkyl_vlasov_projection *proj = &c->info->projection[i];
    double n_s = 0.0, T_s = 0.0, us[3] = {0.0};
    proj->density(t, xn, &n_s, proj->ctx_density);
    proj->temp(t, xn, &T_s, proj->ctx_temp);
    if (proj->V_drift) {
      proj->V_drift(t, xn, us, proj->ctx_V_drift);
    }
    double rho_s = c->mass * n_s;
    fout[0] += rho_s;
    for (int d = 0; d < 3; ++d) {
      fout[1 + d] += rho_s * us[d];
    }
    if (c->num_equations == 5) {
      fout[4] +=
        0.5 * rho_s * (us[0] * us[0] + us[1] * us[1] + us[2] * us[2]) + n_s * T_s * c->inv_gamma_m1;
    }
  }
}

void
vm_fluid_species_source_init(
  struct gkyl_vlasov_app *app, struct vm_fluid_species *fluid_species, struct vm_fluid_source *src
)
{
  const struct gkyl_vlasov_fluid_source *info = &fluid_species->info.source;
  int ncomp = fluid_species->num_equations * app->basis.num_basis;
  src->source = mkarr(app->use_gpu, ncomp, app->local_ext.volume);
  src->source_host = mkarr(false, ncomp, app->local_ext.volume);
  src->calc_bflux = false;
  src->source_length = 1.0;
  src->source_species = 0;
  src->scale = 1.0;

  bool is_euler = fluid_species->eqn_type == GKYL_EQN_EULER ||
                  fluid_species->eqn_type == GKYL_EQN_ISO_EULER;
  struct vm_fluid_source_eval_ctx eval_ctx = {
    .info = info,
    .mass = fluid_species->mass,
    .num_equations = fluid_species->num_equations,
  };
  evalf_t eval = info->profile;
  void *eval_ctx_ptr = info->ctx;
  if (info->source_id == GKYL_PROJ_SOURCE || info->source_id == GKYL_BFLUX_SOURCE) {
    if (!is_euler || info->num_sources < 1) {
      gkyl_exit(
        "vm_fluid_species: a projected or boundary-flux source needs an Euler or isothermal Euler fluid and at least one LTE projection."
      );
    }
    for (int i = 0; i < info->num_sources; ++i) {
      if (info->projection[i].proj_id != GKYL_PROJ_VLASOV_LTE || !info->projection[i].density ||
          !info->projection[i].temp) {
        gkyl_exit(
          "vm_fluid_species: the projections of a fluid source are LTE projections with a density and a temperature."
        );
      }
    }
    if (fluid_species->eqn_type == GKYL_EQN_EULER) {
      eval_ctx.inv_gamma_m1 = 1.0 / (gkyl_wv_euler_gas_gamma(fluid_species->equation) - 1.0);
    }
    eval = eval_source_image;
    eval_ctx_ptr = &eval_ctx;
    if (info->source_id == GKYL_BFLUX_SOURCE) {
      if (!(info->source_length > 0.0)) {
        gkyl_exit("vm_fluid_species: a boundary-flux source needs a positive source_length.");
      }
      src->calc_bflux = true;
      src->source_length = info->source_length;
      src->source_species = vlasov_find_species(app, info->source_species);
      assert(src->source_species); // checked by the app
      src->scale = 0.0;
    }
  } else if (info->source_id != GKYL_FUNC_SOURCE) {
    gkyl_exit("vm_fluid_species: a fluid source is a function, projected or boundary-flux source.");
  }

  // The rate of every component, sampled at the physical coordinates of the position map.
  struct gkyl_proj_on_basis *proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &app->grid,
    .basis = &app->basis,
    .qtype = GKYL_GAUSS_QUAD,
    .num_quad = app->basis.poly_order + 1,
    .num_ret_vals = fluid_species->num_equations,
    .eval = eval,
    .ctx = eval_ctx_ptr,
    .c2p_func = vm_pos_c2p,
    .c2p_func_ctx = app->pos_map,
  });
  gkyl_proj_on_basis_advance(proj, 0.0, &app->local_ext, src->source_host);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(src->source, src->source_host);
  // The rate of the evolved state carries its position-map Jacobian.
  fluid_species->weight_func(app, src->source, &app->local_ext);
}

void
vm_fluid_species_source_calc_rate(gkyl_vlasov_app *app, struct vm_fluid_source *src)
{
  if (src->calc_bflux) {
    src->scale = vlasov_species_wall_loss_rate(app, src->source_species) / src->source_length;
  }
}

void
vm_fluid_species_source_rhs(
  gkyl_vlasov_app *app, const struct vm_fluid_species *species, struct vm_fluid_source *src,
  struct gkyl_array *rhs
)
{
  vm_fluid_species_source_calc_rate(app, src);
  gkyl_array_accumulate_range(rhs, src->scale, src->source, &app->local);
}

void
vm_fluid_species_source_write(
  gkyl_vlasov_app *app, const struct vm_fluid_species *fluid_species,
  const struct vm_fluid_source *src, double tm, int frame
)
{
  if (!fluid_species->info.source.write_source) {
    return;
  }
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new((struct vlasov_output_meta){
    .frame = frame,
    .stime = tm,
    .poly_order = app->poly_order,
    .basis_type = app->basis.id,
  });
  // The physical rate of the profile (unscaled by the wall losses).
  const char *fmt = "%s-%s_source_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, fluid_species->name, frame);
  char fileNm[sz + 1];
  snprintf(fileNm, sizeof fileNm, fmt, app->name, fluid_species->name, frame);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, src->source_host, fileNm);
  vlasov_array_meta_release(mt);
}

void
vm_fluid_species_source_release(const struct gkyl_vlasov_app *app, const struct vm_fluid_source *src)
{
  gkyl_array_release(src->source);
  gkyl_array_release(src->source_host);
}
