#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>

// Input packaged as a struct
struct gkyl_mom_pkpm_inp {
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions (p=1 serendipity).
  const struct gkyl_basis *phase_basis; // Phase-space basis functions (p=1 x p=2 hybrid).
  double mass; // Mass of species (the moments are mass weighted).
  bool is_diag; // bool to determine if computing the diagnostic moments.
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create new Vlasov parallel-kinetic-perpendicular-moment (pkpm) moment type object.
 * Can either compute moments needed for update (rho, p_par, p_perp) or diagnostic moments
 * (rho, p_par, p_perp, q_par, q_perp, r_parpar, r_parperp)
 *
 * @param inp Input parameters defined in gkyl_mom_pkpm_inp struct.
 */
struct gkyl_mom_type *gkyl_mom_pkpm_inew(const struct gkyl_mom_pkpm_inp *inp);

/**
 * Create new Vlasov moment type object on NV-GPU: see new() method
 * above for documentation.
 */
struct gkyl_mom_type *gkyl_mom_pkpm_cu_dev_inew(const struct gkyl_mom_pkpm_inp *inp);
