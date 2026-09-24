#include <gkyl_prim_cross_m0deltas.h>
#include <gkyl_dg_bin_ops_priv.h>

struct gkyl_prim_cross_m0deltas {
  bool normNu; // Whether nu=nu(x,t).
  const struct gkyl_basis *basis; // DG basis representing m0 and nu.
  const struct gkyl_range *range; // Range to perform operation in.
  struct gkyl_dg_bin_op_mem *mem; // Memory needed for weak division.
  double betap1T2; // Greene's beta plus 1 times 2: (beta+1)*2.
  bool use_gpu; // Whether to run on GPU.
};

#ifdef GKYL_HAVE_CUDA
/**
 * On the NVIDIA GPU, compute 
 *   n_s*delta_s*(beta+1) = 2*(beta+1) * n_s * m_r * n_r * nu_rs / (m_s * n_s * nu_sr + m_r * n_r * nu_rs)
 * if collision frequency is constant, or
 *   nu_sr*n_s*delta_s*(beta+1) = 2*(beta+1) * n_s * nu_sr * m_r * n_r * nu_rs / (m_s * n_s * nu_sr + m_r * n_r * nu_rs)
 * if the collision frequency varies in space and time.
 *
 * @param up Struct defining this updater..
 * @param massself Mass of this species, m_s.
 * @param m0self Number density of this species, m0_s.
 * @param nuself Cross collision frequency of this species, nu_sr.
 * @param massother Mass of the other species, m_r.
 * @param m0other Number density of the other species, m0_r.
 * @param nuother Cross collision frequency of the other species, nu_rs.
 * @param out Output array.
 * @return New updater pointer.
 */
void gkyl_prim_cross_m0deltas_advance_cu(gkyl_prim_cross_m0deltas *up,
  double massself, const struct gkyl_array* m0self, const struct gkyl_array* nuself,
  double massother, const struct gkyl_array* m0other, const struct gkyl_array* nuother,
  struct gkyl_array* out);
#endif
