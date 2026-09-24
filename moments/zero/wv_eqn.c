#include <gkyl_wv_eqn.h>
#include <gkyl_alloc_flags_priv.h>

// these ensure inline functions are defined only once

extern inline double gkyl_wv_eqn_waves(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type, const double *delta, const double *ql, const double *qr, const double phil, const double phir, double *waves, double *speeds);

extern inline void gkyl_wv_eqn_qfluct(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type, const double *ql, const double *qr, const double phil, const double phir, const double *waves, const double *speeds, double *amdq, double *apdq);

extern inline void gkyl_wv_eqn_ffluct(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type, const double *ql, const double *qr, const double phil, const double phir, const double *waves, const double *speeds, double *amdq, double *apdq);

extern inline double gkyl_wv_eqn_flux_jump(const struct gkyl_wv_eqn *eqn,
  const double *ql, const double *qr, double *flux_jump);

extern inline bool gkyl_wv_eqn_check_inv(const struct gkyl_wv_eqn *eqn, const double *q);

extern inline double gkyl_wv_eqn_max_speed(const struct gkyl_wv_eqn *eqn, const double *q);

extern inline void gkyl_wv_eqn_rotate_to_local(const struct gkyl_wv_eqn *eqn,
  const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qglobal, double *GKYL_RESTRICT qlocal);

extern inline void gkyl_wv_eqn_rotate_to_global(const struct gkyl_wv_eqn *eqn,
  const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qlocal, double *GKYL_RESTRICT qglobal);

extern inline void gkyl_wv_eqn_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout);

bool
gkyl_wv_eqn_is_cu_dev(const struct gkyl_wv_eqn *eqn)
{
  return GKYL_IS_CU_ALLOC(eqn->flags);
}

struct gkyl_wv_eqn*
gkyl_wv_eqn_acquire(const struct gkyl_wv_eqn* eqn)
{
  gkyl_ref_count_inc(&eqn->ref_count);
  return (struct gkyl_wv_eqn*) eqn;
}

void
gkyl_wv_eqn_release(const struct gkyl_wv_eqn* eqn)
{
  gkyl_ref_count_dec(&eqn->ref_count);
}
