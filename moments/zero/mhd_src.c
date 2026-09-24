#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_mhd_src.h>


// Makes indexing cleaner
#define DN (0)
#define MX (1)
#define MY (2)
#define MZ (3)
#define ER (4)
#define BX (5)
#define BY (6)
#define BZ (7)
#define PSI_GLM (8)

gkyl_mhd_src *gkyl_mhd_src_new(struct gkyl_mhd_src_inp inp,
                               const struct gkyl_range *local_ext) {
  gkyl_mhd_src *up = gkyl_calloc(1, sizeof(gkyl_mhd_src));

  up->grid = *(inp.grid);
  up->ndim = up->grid.ndim;
  up->divergence_constraint = inp.divergence_constraint;
  up->glm_ch = inp.glm_ch;
  up->glm_alpha = inp.glm_alpha;
  up->dxyz_min = inp.dxyz_min;

  if (up->divergence_constraint >= GKYL_MHD_DIVB_EIGHT_WAVES
      && up->divergence_constraint <= GKYL_MHD_DIVB_GLM) {
    up->divB_array = gkyl_array_new(GKYL_DOUBLE, 1, local_ext->volume);
  }

  if (up->divergence_constraint == GKYL_MHD_DIVB_GLM) {
    up->B_dot_gradPsi_array = gkyl_array_new(GKYL_DOUBLE, 1, local_ext->volume);
  }

  return up;
}

///////////////////////////
// SOURCE UPDATE HELPERS //
///////////////////////////

static void calc_divB(const gkyl_mhd_src *up,
                      const struct gkyl_range *update_range,
                      struct gkyl_array *q_array) {
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(update_range, iter.idx);
    double *q = gkyl_array_fetch(q_array, lidx);
    double *divB = gkyl_array_fetch(up->divB_array, lidx);

    int idxl[3], idxr[3];
    
    double my_divB = 0.0;
    for(int d=0; d<up->grid.ndim; ++d) {
      gkyl_copy_int_arr(up->grid.ndim, iter.idx, idxl);
      gkyl_copy_int_arr(up->grid.ndim, iter.idx, idxr);
      idxl[d]--;
      idxr[d]++;

      lidx = gkyl_range_idx(update_range, idxl);
      double *ql = gkyl_array_fetch(q_array, lidx);
      lidx = gkyl_range_idx(update_range, idxr);
      double *qr = gkyl_array_fetch(q_array, lidx);

      double dx = 2.0 * up->grid.dx[d];
      int Bn = BX + d;
      my_divB += (qr[Bn] - ql[Bn]) / dx;
    }
    divB[0] = my_divB;
  }
}

static void calc_B_dot_gradPsi(const gkyl_mhd_src *up,
                               const struct gkyl_range *update_range,
                               struct gkyl_array *q_array) {
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(update_range, iter.idx);
    double *q = gkyl_array_fetch(q_array, lidx);
    double *B_dot_gradPsi = gkyl_array_fetch(up->B_dot_gradPsi_array, lidx);

    int idxl[3], idxr[3];
    
    double my_B_dot_gradPsi = 0.0;
    for(int d=0; d<up->grid.ndim; ++d) {
      gkyl_copy_int_arr(up->grid.ndim, iter.idx, idxl);
      gkyl_copy_int_arr(up->grid.ndim, iter.idx, idxr);
      idxl[d]--;
      idxr[d]++;

      lidx = gkyl_range_idx(update_range, idxl);
      double *ql = gkyl_array_fetch(q_array, lidx);
      lidx = gkyl_range_idx(update_range, idxr);
      double *qr = gkyl_array_fetch(q_array, lidx);

      double dx = 2.0 * up->grid.dx[d];
      int Bn = BX + d;
      my_B_dot_gradPsi += q[Bn] * (qr[PSI_GLM] - ql[PSI_GLM]) / dx;
    }
    B_dot_gradPsi[0] = my_B_dot_gradPsi;
  }
}

////////////////////////////
// VARIOUS SOURCE UPDATES //
////////////////////////////

static void gkyl_mhd_src_eight_wave(const gkyl_mhd_src *up, double dt,
                                    const struct gkyl_range *update_range,
                                    struct gkyl_array *q_array,
                                    const struct gkyl_array *acc_array) {
  // Powell et al., JCP (1999), 10.1006/jcph.1999.6299
  calc_divB(up, update_range, q_array);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(update_range, iter.idx);
    double *q = gkyl_array_fetch(q_array, lidx);
    double *divB =gkyl_array_fetch(up->divB_array, lidx);

    double ux = q[MX] / q[DN], uy = q[MY] / q[DN], uz = q[MZ] / q[DN];
    double Bx = q[BX], By = q[BY], Bz = q[BZ];

    // Equation (25)
    q[MX] -= dt * divB[0] * Bx;
    q[MY] -= dt * divB[0] * By;
    q[MZ] -= dt * divB[0] * Bz;
    q[BX] -= dt * divB[0] * ux;
    q[BY] -= dt * divB[0] * uy;
    q[BZ] -= dt * divB[0] * uz;
    q[ER] -= dt * divB[0] * (ux*Bx + uy*By + uz*Bz);
  }
}

static void gkyl_mhd_src_glm(const gkyl_mhd_src *up, double dt,
                             const struct gkyl_range *update_range,
                             struct gkyl_array *q_array,
                             const struct gkyl_array *acc_array) {
  // Dedner et al., JCP (2002), 10.1006/jcph.2001.6961
  // Mignone & Tzeferacos, JCP (2010), 10.1016/j.jcp.2009.11.026
  calc_divB(up, update_range, q_array);
  calc_B_dot_gradPsi(up, update_range, q_array);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(update_range, iter.idx);
    double *q = gkyl_array_fetch(q_array, lidx);
    double *divB =gkyl_array_fetch(up->divB_array, lidx);
    double *B_dot_gradPsi = gkyl_array_fetch(up->B_dot_gradPsi_array, lidx);

    double ch = up->glm_ch;
    double alpha = up->glm_alpha;

    // Dedner's equation (24e) or Mignone's equation (27)
    double rate = alpha * ch / up->dxyz_min;
    q[PSI_GLM] *= exp(- rate * dt);

    double Bx = q[BX], By = q[BY], Bz = q[BZ];
    // Dedner's equations (24b) and (24d)
    q[MX] -= dt * divB[0] * Bx;
    q[MY] -= dt * divB[0] * By;
    q[MZ] -= dt * divB[0] * Bz;
    q[ER] -= dt * B_dot_gradPsi[0];
  }
}

/////////////////////////////
// SOURCE UPDATES MAIN API //
/////////////////////////////

void gkyl_mhd_src_advance(const gkyl_mhd_src *up, double dt,
                          const struct gkyl_range *update_range,
                          struct gkyl_array *q_array,
                          const struct gkyl_array *acc_array) {
  switch (up->divergence_constraint) {
    case GKYL_MHD_DIVB_NONE:
      break;

    case GKYL_MHD_DIVB_EIGHT_WAVES:
      gkyl_mhd_src_eight_wave(up, dt, update_range, q_array, acc_array);
      break;

    case GKYL_MHD_DIVB_GLM:
      gkyl_mhd_src_glm(up, dt, update_range, q_array, acc_array);
      break;
  }
}

double gkyl_mhd_src_calc_divB(const gkyl_mhd_src *up,
                              const struct gkyl_range *update_range,
                              struct gkyl_array *q_array) {
  calc_divB(up, update_range, q_array);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);

  double total = 0.0;
  int count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(update_range, iter.idx);
    double *divB = gkyl_array_fetch(up->divB_array, lidx);
    total += fabs(divB[0]);
    count++;
  }
  return total / count;
}

void gkyl_mhd_src_release(gkyl_mhd_src *up) {
  if (up->divB_array)
    gkyl_array_release(up->divB_array);

  if (up->B_dot_gradPsi_array)
    gkyl_array_release(up->B_dot_gradPsi_array);

  gkyl_free(up);
}

////////////////////////////////
// MEMBER GETTERS AND SETTERS //
////////////////////////////////

void
gkyl_mhd_src_set_glm_ch(struct gkyl_mhd_src* up, const double glm_ch)
{
  up->glm_ch = glm_ch;
}

