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
 * Create new Vlasov parallel-kinetic-perpendicular-moment (pkpm) moment type object: the
 * mass-weighted w moments of the distributions [F, G, M0_perp, N_perp], with Mk the k-th w
 * moment of F and Mk_G that of G. The update set has 18 moments, [M0, M2, M0_G, M1, M3, M1_G,
 * int (1, w, w^2) M0_perp dw (x, y, z each), int N_perp dw (x, y, z)]; the diagnostic set
 * (is_diag) has 8, [M0, M1, M2, M0_G, M3, M1_G, M4, M2_G].
 *
 * @param inp Input parameters defined in gkyl_mom_pkpm_inp struct.
 * @return New moment type object.
 */
struct gkyl_mom_type *gkyl_mom_pkpm_inew(const struct gkyl_mom_pkpm_inp *inp);

/**
 * Create new Vlasov moment type object on NV-GPU: see new() method
 * above for documentation.
 */
struct gkyl_mom_type *gkyl_mom_pkpm_cu_dev_inew(const struct gkyl_mom_pkpm_inp *inp);
