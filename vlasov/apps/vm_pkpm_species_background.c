#include <gkyl_vlasov_priv.h>

// Held background of a PKPM species (First_Fourier_Harmonic_Notes/background_density_gradient.tex):
// the eight normalized distributions (one function projection, the background of the kinetic
// block) and the fluid state (the background of the fluid block), both or neither. The force acts on
// the fluid only, since the kinetic equations do not see an applied acceleration; both aspects keep a
// remainder. The force, the remainders, the stage methods and the release are those of the fluid
// and the kinetic aspect (vm_fluid_species_background.c, vm_species_background.c). The capture is in
// gkyl_vlasov_app_apply_ic (vlasov.c).

void
vm_pkpm_species_background_init(struct gkyl_vlasov_app *app, struct vlasov_species *s)
{
  struct vm_species *kin = s->kinetic;
  struct vm_fluid_species *fl = s->fluid;
  const struct gkyl_vlasov_projection *inp = &kin->info.background;

  vm_species_background_init_none(&kin->bg);
  vm_fluid_species_background_init_none(&fl->bg);
  bool has_dist = inp->proj_id == GKYL_PROJ_FUNC && inp->func != 0;
  bool has_fluid = fl->info.background != 0;
  if (has_dist != has_fluid) {
    gkyl_exit(
      "vm_pkpm_species_background: a held background needs both the distributions (a GKYL_PROJ_FUNC background of the kinetic block) and the fluid state (the background of the fluid block)."
    );
  }
  if (has_dist) {
    vm_species_background_alloc_rhs(app, kin, &kin->bg);
    vm_fluid_species_background_alloc(app, fl, &fl->bg);
  }
}

void
vm_pkpm_species_background_project(gkyl_vlasov_app *app, struct vlasov_species *s, double t0)
{
  struct vm_species *kin = s->kinetic;
  struct vm_fluid_species *fl = s->fluid;
  vm_pkpm_kinetic_species_project(
    app, kin, kin->info.background.func, kin->info.background.ctx_func, t0
  );
  vm_pkpm_fluid_species_project(app, fl, fl->info.background, fl->info.background_ctx, t0);
  // The background is what the limiter leaves of it, as the initial condition is (the capture
  // discards the record).
  vm_pkpm_species_apply_positivity(app, fl, fl->fluid);
}
