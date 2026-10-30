#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_mom_type.h>
#include <gkyl_ref_count.h>
#include <gkyl_mom_vlasov_kernels.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_vlasov_velocity_map.h>
#include <gkyl_vlasov_kernel_sets.h>

struct mom_type_vlasov {
  struct gkyl_mom_type momt;
  bool use_vmap; // If using nonuniform velocity space grids.
  int hamil_dim; // Dimensionality of Hamiltonian.
  int hamil_offset; // Offset for indexing Hamiltonian from phase-space index.
  struct gkyl_range
    hamil_range; // Range for indexing Hamiltonian (either velocity-space range or full phase-space range).
  const struct gkyl_array *hamil; // Hamiltonian utilized to compute certain moments.
  struct gkyl_range vel_range; // Range for indexing velocity-space Jacobian and velocity map.
  const struct gkyl_vlasov_velocity_map
    *vel_map; // Velocity-space mapping object (acquired host-side; 0 if not given).
  const struct gkyl_array
    *vmap; // Velocity-space mapping (borrowed from vel_map; device pointer on GPUs).
  const struct gkyl_array
    *jacob_vel; // Velocity-space Jacobian (borrowed from vel_map; device pointer on GPUs).
  double v_thresh; // Threshold velocity for integration of moments over a subset of the domain.
  double f_thresh; // Threshold for whether we accumulate moment over subset of the domain.
  struct gkyl_range conf_range; // Range for indexing the conf-space Jacobian of the position map.
  const struct gkyl_vlasov_position_map
    *pos_map; // Configuration-space mapping object (acquired host-side; 0 if not mapped).
  const struct gkyl_array *jacob_pos; // Conf-space Jacobian (borrowed from pos_map; device pointer on
  // GPUs; 0 if not mapped).
};

// Threshold on the stored f used by the M0_UPPER and M0_LOWER kernels in the cell with phase-space
// index idx. f_thresh is a threshold on the physical f: the kernels weight it by the velocity-space
// Jacobian, and on a mapped conf mesh the stored f also carries the (per-cell constant) conf Jacobian.
GKYL_CU_DH static inline double
mom_vlasov_f_thresh(const struct mom_type_vlasov *mom_vlasov, const int *idx)
{
  if (!mom_vlasov->jacob_pos) {
    return mom_vlasov->f_thresh;
  }
  long cidx = gkyl_range_idx(&mom_vlasov->conf_range, idx);
  return mom_vlasov->f_thresh * ((const double *)gkyl_array_cfetch(mom_vlasov->jacob_pos, cidx))[0];
}

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
GKYL_CU_D static struct {
  int vdim[4];
} cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1, 0, 1, 2}, // 1x kernel indices
  {-1, 3, 4, 5}, // 2x kernel indices
  {-1, -1, -1, 6} // 3x kernel indices
};

// for use in kernel tables
typedef struct {
  momf_t kernels[4];
} gkyl_vlasov_mom_kern_list;

GKYL_CU_DH static void
kernel_mom_vlasov_M0_1x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_1x1v_ser_p1(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_1x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_1x1v_ser_p2(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_1x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_1x1v_ser_p3(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_1x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_1x2v_ser_p1(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_1x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_1x2v_ser_p2(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_1x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_1x3v_ser_p1(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_1x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_1x3v_ser_p2(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_2x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_2x1v_ser_p1(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_2x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_2x1v_ser_p2(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_2x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_2x1v_ser_p3(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_2x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_2x2v_ser_p1(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_2x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_2x2v_ser_p2(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_2x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_2x3v_ser_p1(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_2x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_2x3v_ser_p2(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_3x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_3x3v_ser_p1(xc, dx, idx, f, out);
}

// M0 kernel list (Serendipity basis)
GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_m0_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_M0_1x1v_ser_p1, kernel_mom_vlasov_M0_1x1v_ser_p2,
   kernel_mom_vlasov_M0_1x1v_ser_p3}, // 0
  {NULL, kernel_mom_vlasov_M0_1x2v_ser_p1, kernel_mom_vlasov_M0_1x2v_ser_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_M0_1x3v_ser_p1, kernel_mom_vlasov_M0_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_M0_2x1v_ser_p1, kernel_mom_vlasov_M0_2x1v_ser_p2,
   kernel_mom_vlasov_M0_2x1v_ser_p3}, // 3
  {NULL, kernel_mom_vlasov_M0_2x2v_ser_p1, kernel_mom_vlasov_M0_2x2v_ser_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_M0_2x3v_ser_p1, kernel_mom_vlasov_M0_2x3v_ser_p2, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_M0_3x3v_ser_p1, NULL, NULL} // 6
};

GKYL_CU_DH static void
kernel_mom_vlasov_M0_1x1v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_1x1v_tensor_p2(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_1x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_1x1v_tensor_p1(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_1x1v_tensor_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_1x1v_tensor_p3(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_1x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_1x2v_tensor_p2(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_1x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_1x2v_tensor_p1(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_1x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_1x3v_tensor_p2(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_1x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_1x3v_tensor_p1(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_2x1v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_2x1v_tensor_p2(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_2x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_2x1v_tensor_p1(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_2x1v_tensor_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_2x1v_tensor_p3(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_2x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_2x2v_tensor_p2(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_2x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_2x2v_tensor_p1(xc, dx, idx, f, out);
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_2x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_2x3v_tensor_p1(xc, dx, idx, f, out);
}

#ifdef GKYL_BUILD_VLASOV_TEN_2X3V_P2
GKYL_CU_DH static void
kernel_mom_vlasov_M0_2x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_2x3v_tensor_p2(xc, dx, idx, f, out);
}
#endif

GKYL_CU_DH static void
kernel_mom_vlasov_M0_3x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  return mom_vlasov_M0_3x3v_tensor_p1(xc, dx, idx, f, out);
}

// M0 kernel list (Tensor basis)
GKYL_CU_D static const gkyl_vlasov_mom_kern_list tensor_m0_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_M0_1x1v_tensor_p1, kernel_mom_vlasov_M0_1x1v_tensor_p2,
   kernel_mom_vlasov_M0_1x1v_tensor_p3}, // 0
  {NULL, kernel_mom_vlasov_M0_1x2v_tensor_p1, kernel_mom_vlasov_M0_1x2v_tensor_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_M0_1x3v_tensor_p1, kernel_mom_vlasov_M0_1x3v_tensor_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_M0_2x1v_tensor_p1, kernel_mom_vlasov_M0_2x1v_tensor_p2,
   kernel_mom_vlasov_M0_2x1v_tensor_p3}, // 3
  {NULL, kernel_mom_vlasov_M0_2x2v_tensor_p1, kernel_mom_vlasov_M0_2x2v_tensor_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_M0_2x3v_tensor_p1, GKYL_TEN_2X3V_P2(kernel_mom_vlasov_M0_2x3v_tensor_p2), NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_M0_3x3v_tensor_p1, NULL, NULL} // 6
};

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_1x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_1x1v_ser_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_1x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_1x1v_ser_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_1x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_1x1v_ser_p3(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_1x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_1x2v_ser_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_1x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_1x2v_ser_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_1x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_1x3v_ser_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_1x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_1x3v_ser_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_2x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_2x1v_ser_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_2x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_2x1v_ser_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_2x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_2x1v_ser_p3(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_2x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_2x2v_ser_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_2x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_2x2v_ser_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_2x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_2x3v_ser_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_2x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_2x3v_ser_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_3x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_3x3v_ser_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

// M2ij kernel list (Serendipity basis)
GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_m2ij_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_M2ij_1x1v_ser_p1, kernel_mom_vlasov_M2ij_1x1v_ser_p2,
   kernel_mom_vlasov_M2ij_1x1v_ser_p3}, // 0
  {NULL, kernel_mom_vlasov_M2ij_1x2v_ser_p1, kernel_mom_vlasov_M2ij_1x2v_ser_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_M2ij_1x3v_ser_p1, kernel_mom_vlasov_M2ij_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_M2ij_2x1v_ser_p1, kernel_mom_vlasov_M2ij_2x1v_ser_p2,
   kernel_mom_vlasov_M2ij_2x1v_ser_p3}, // 3
  {NULL, kernel_mom_vlasov_M2ij_2x2v_ser_p1, kernel_mom_vlasov_M2ij_2x2v_ser_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_M2ij_2x3v_ser_p1, kernel_mom_vlasov_M2ij_2x3v_ser_p2, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_M2ij_3x3v_ser_p1, NULL, NULL} // 6
};

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_1x1v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_1x1v_tensor_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_1x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_1x1v_tensor_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_1x1v_tensor_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_1x1v_tensor_p3(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_1x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_1x2v_tensor_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_1x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_1x2v_tensor_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_1x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_1x3v_tensor_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_1x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_1x3v_tensor_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_2x1v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_2x1v_tensor_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_2x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_2x1v_tensor_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_2x1v_tensor_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_2x1v_tensor_p3(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_2x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_2x2v_tensor_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_2x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_2x2v_tensor_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_2x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_2x3v_tensor_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

#ifdef GKYL_BUILD_VLASOV_TEN_2X3V_P2
GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_2x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_2x3v_tensor_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}
#endif

GKYL_CU_DH static void
kernel_mom_vlasov_M2ij_3x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M2ij_3x3v_tensor_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

// M2ij kernel list (Tensor basis)
GKYL_CU_D static const gkyl_vlasov_mom_kern_list tensor_m2ij_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_M2ij_1x1v_tensor_p1, kernel_mom_vlasov_M2ij_1x1v_tensor_p2,
   kernel_mom_vlasov_M2ij_1x1v_tensor_p3}, // 0
  {NULL, kernel_mom_vlasov_M2ij_1x2v_tensor_p1, kernel_mom_vlasov_M2ij_1x2v_tensor_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_M2ij_1x3v_tensor_p1, kernel_mom_vlasov_M2ij_1x3v_tensor_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_M2ij_2x1v_tensor_p1, kernel_mom_vlasov_M2ij_2x1v_tensor_p2,
   kernel_mom_vlasov_M2ij_2x1v_tensor_p3}, // 3
  {NULL, kernel_mom_vlasov_M2ij_2x2v_tensor_p1, kernel_mom_vlasov_M2ij_2x2v_tensor_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_M2ij_2x3v_tensor_p1, GKYL_TEN_2X3V_P2(kernel_mom_vlasov_M2ij_2x3v_tensor_p2), NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_M2ij_3x3v_tensor_p1, NULL, NULL} // 6
};

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_1x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_1x1v_ser_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_1x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_1x1v_ser_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_1x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_1x1v_ser_p3(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_1x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_1x2v_ser_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_1x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_1x2v_ser_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_1x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_1x3v_ser_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_1x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_1x3v_ser_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_2x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_2x1v_ser_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_2x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_2x1v_ser_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_2x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_2x1v_ser_p3(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_2x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_2x2v_ser_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_2x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_2x2v_ser_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_2x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_2x3v_ser_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_2x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_2x3v_ser_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_3x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_3x3v_ser_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

// M3ijk kernel list (Serendipity basis)
GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_m3ijk_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_M3ijk_1x1v_ser_p1, kernel_mom_vlasov_M3ijk_1x1v_ser_p2,
   kernel_mom_vlasov_M3ijk_1x1v_ser_p3}, // 0
  {NULL, kernel_mom_vlasov_M3ijk_1x2v_ser_p1, kernel_mom_vlasov_M3ijk_1x2v_ser_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_M3ijk_1x3v_ser_p1, kernel_mom_vlasov_M3ijk_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_M3ijk_2x1v_ser_p1, kernel_mom_vlasov_M3ijk_2x1v_ser_p2,
   kernel_mom_vlasov_M3ijk_2x1v_ser_p3}, // 3
  {NULL, kernel_mom_vlasov_M3ijk_2x2v_ser_p1, kernel_mom_vlasov_M3ijk_2x2v_ser_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_M3ijk_2x3v_ser_p1, kernel_mom_vlasov_M3ijk_2x3v_ser_p2, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_M3ijk_3x3v_ser_p1, NULL, NULL} // 6
};

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_1x1v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_1x1v_tensor_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_1x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_1x1v_tensor_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_1x1v_tensor_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_1x1v_tensor_p3(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_1x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_1x2v_tensor_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_1x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_1x2v_tensor_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_1x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_1x3v_tensor_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_1x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_1x3v_tensor_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_2x1v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_2x1v_tensor_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_2x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_2x1v_tensor_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_2x1v_tensor_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_2x1v_tensor_p3(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_2x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_2x2v_tensor_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_2x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_2x2v_tensor_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_2x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_2x3v_tensor_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

#ifdef GKYL_BUILD_VLASOV_TEN_2X3V_P2
GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_2x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_2x3v_tensor_p2(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}
#endif

GKYL_CU_DH static void
kernel_mom_vlasov_M3ijk_3x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M3ijk_3x3v_tensor_p1(
    xc, dx, idx, mom_vlasov->vmap ? (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx) : 0,
    f, out
  );
}

// M3ijk kernel list (Tensor basis)
GKYL_CU_D static const gkyl_vlasov_mom_kern_list tensor_m3ijk_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_M3ijk_1x1v_tensor_p1, kernel_mom_vlasov_M3ijk_1x1v_tensor_p2,
   kernel_mom_vlasov_M3ijk_1x1v_tensor_p3}, // 0
  {NULL, kernel_mom_vlasov_M3ijk_1x2v_tensor_p1, kernel_mom_vlasov_M3ijk_1x2v_tensor_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_M3ijk_1x3v_tensor_p1, kernel_mom_vlasov_M3ijk_1x3v_tensor_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_M3ijk_2x1v_tensor_p1, kernel_mom_vlasov_M3ijk_2x1v_tensor_p2,
   kernel_mom_vlasov_M3ijk_2x1v_tensor_p3}, // 3
  {NULL, kernel_mom_vlasov_M3ijk_2x2v_tensor_p1, kernel_mom_vlasov_M3ijk_2x2v_tensor_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_M3ijk_2x3v_tensor_p1, GKYL_TEN_2X3V_P2(kernel_mom_vlasov_M3ijk_2x3v_tensor_p2), NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_M3ijk_3x3v_tensor_p1, NULL, NULL} // 6
};

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_1x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_1x1v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_1x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_1x1v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_1x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_1x1v_ser_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_1x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_1x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M1i_1x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M1i_1x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_1x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_1x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M1i_1x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M1i_1x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_1x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_1x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M1i_1x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M1i_1x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_1x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_1x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M1i_1x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M1i_1x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_2x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_2x1v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_2x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_2x1v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_2x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_2x1v_ser_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_2x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_2x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M1i_2x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M1i_2x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_2x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_2x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M1i_2x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M1i_2x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_2x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_2x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M1i_2x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M1i_2x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_2x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_2x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M1i_2x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M1i_2x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_3x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_3x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M1i_3x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M1i_3x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

// M1i (dH/dv moment, velocity-space Hamiltonian) kernel list (Serendipity basis)
GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_hamil_vel_dense_m1i_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M1i_1x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M1i_1x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_M1i_1x1v_ser_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M1i_1x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M1i_1x2v_ser_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M1i_1x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M1i_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M1i_2x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M1i_2x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_M1i_2x1v_ser_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M1i_2x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M1i_2x2v_ser_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M1i_2x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M1i_2x3v_ser_p2, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M1i_3x3v_ser_p1, NULL, NULL} // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_hamil_vel_sparse_m1i_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M1i_1x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M1i_1x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_M1i_1x1v_ser_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M1i_1x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_M1i_1x2v_ser_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M1i_1x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_M1i_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M1i_2x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M1i_2x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_M1i_2x1v_ser_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M1i_2x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_M1i_2x2v_ser_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M1i_2x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_M1i_2x3v_ser_p2, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M1i_3x3v_ser_p1, NULL, NULL} // 6
};

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_1x1v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_1x1v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_1x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_1x1v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_1x1v_tensor_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_1x1v_tensor_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_1x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_1x2v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_1x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_1x2v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M1i_1x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M1i_1x2v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M1i_1x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M1i_1x2v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_1x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_1x3v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_1x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_1x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M1i_1x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M1i_1x3v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M1i_1x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M1i_1x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_2x1v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_2x1v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_2x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_2x1v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_2x1v_tensor_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_2x1v_tensor_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_2x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_2x2v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_2x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_2x2v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_2x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_2x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

#ifdef GKYL_BUILD_VLASOV_TEN_2X3V_P2
GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_2x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_2x3v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}
#endif

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M1i_3x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M1i_3x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M1i_2x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M1i_2x2v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M1i_2x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M1i_2x2v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M1i_2x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M1i_2x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

#ifdef GKYL_BUILD_VLASOV_TEN_2X3V_P2
GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M1i_2x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M1i_2x3v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}
#endif

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M1i_3x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M1i_3x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

// M1i (dH/dv moment, velocity-space Hamiltonian) kernel list (Tensor basis)
GKYL_CU_D static const gkyl_vlasov_mom_kern_list tensor_hamil_vel_dense_m1i_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M1i_1x1v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_M1i_1x1v_tensor_p2,
   kernel_mom_vlasov_hamil_vel_dense_M1i_1x1v_tensor_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M1i_1x2v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_M1i_1x2v_tensor_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M1i_1x3v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_M1i_1x3v_tensor_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M1i_2x1v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_M1i_2x1v_tensor_p2,
   kernel_mom_vlasov_hamil_vel_dense_M1i_2x1v_tensor_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M1i_2x2v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_M1i_2x2v_tensor_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M1i_2x3v_tensor_p1, GKYL_TEN_2X3V_P2(kernel_mom_vlasov_hamil_vel_dense_M1i_2x3v_tensor_p2), NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M1i_3x3v_tensor_p1, NULL, NULL} // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D static const gkyl_vlasov_mom_kern_list tensor_hamil_vel_sparse_m1i_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M1i_1x1v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_M1i_1x1v_tensor_p2,
   kernel_mom_vlasov_hamil_vel_dense_M1i_1x1v_tensor_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M1i_1x2v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_sparse_M1i_1x2v_tensor_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M1i_1x3v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_sparse_M1i_1x3v_tensor_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M1i_2x1v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_M1i_2x1v_tensor_p2,
   kernel_mom_vlasov_hamil_vel_dense_M1i_2x1v_tensor_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M1i_2x2v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_sparse_M1i_2x2v_tensor_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M1i_2x3v_tensor_p1, GKYL_TEN_2X3V_P2(kernel_mom_vlasov_hamil_vel_sparse_M1i_2x3v_tensor_p2), NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M1i_3x3v_tensor_p1, NULL, NULL} // 6
};

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_1x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_1x1v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_1x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_1x1v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_1x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_1x1v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_1x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_1x1v_ser_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_1x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_1x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_1x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_1x2v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_1x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_1x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_1x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_1x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_1x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_1x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_1x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_1x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_2x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_2x1v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_2x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_2x1v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_2x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_2x1v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_2x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_2x1v_ser_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_2x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_2x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_2x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_2x2v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_2x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_2x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_2x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_2x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_2x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_2x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_3x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_3x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M1i_3x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M1i_3x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

// M1i (dH/dv moment, phase-space Hamiltonian) kernel list (Serendipity basis)
GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_hamil_phase_m1i_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_M1i_1x1v_ser_p1,
   kernel_mom_vlasov_hamil_phase_M1i_1x1v_ser_p2,
   kernel_mom_vlasov_hamil_phase_M1i_1x1v_ser_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_phase_M1i_1x2v_ser_p1,
   kernel_mom_vlasov_hamil_phase_M1i_1x2v_ser_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_phase_M1i_1x3v_ser_p1,
   kernel_mom_vlasov_hamil_phase_M1i_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_M1i_2x1v_ser_p1,
   kernel_mom_vlasov_hamil_phase_M1i_2x1v_ser_p2,
   kernel_mom_vlasov_hamil_phase_M1i_2x1v_ser_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_phase_M1i_2x2v_ser_p1,
   kernel_mom_vlasov_hamil_phase_M1i_2x2v_ser_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_phase_M1i_2x3v_ser_p1, NULL, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_M1i_3x3v_ser_p1, NULL, NULL} // 6
};

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_1x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_1x1v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_1x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_1x1v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_1x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_1x1v_ser_p3(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_1x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_1x2v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M2_1x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M2_1x2v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_1x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_1x2v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M2_1x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M2_1x2v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_1x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_1x3v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M2_1x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M2_1x3v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_1x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_1x3v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M2_1x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M2_1x3v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_2x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_2x1v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_2x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_2x1v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_2x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_2x1v_ser_p3(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_2x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_2x2v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M2_2x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M2_2x2v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_2x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_2x2v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M2_2x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M2_2x2v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_2x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_2x3v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M2_2x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M2_2x3v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_2x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_2x3v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M2_2x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M2_2x3v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_3x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_3x3v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M2_3x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M2_3x3v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

// M2 (H moment, velocity-space Hamiltonian) kernel list (Serendipity basis)
GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_hamil_vel_dense_m2_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M2_1x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M2_1x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_M2_1x1v_ser_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M2_1x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M2_1x2v_ser_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M2_1x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M2_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M2_2x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M2_2x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_M2_2x1v_ser_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M2_2x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M2_2x2v_ser_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M2_2x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M2_2x3v_ser_p2, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M2_3x3v_ser_p1, NULL, NULL} // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_hamil_vel_sparse_m2_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M2_1x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M2_1x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_M2_1x1v_ser_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M2_1x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_M2_1x2v_ser_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M2_1x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_M2_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M2_2x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M2_2x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_M2_2x1v_ser_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M2_2x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_M2_2x2v_ser_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M2_2x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_M2_2x3v_ser_p2, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M2_3x3v_ser_p1, NULL, NULL} // 6
};

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_1x1v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_1x1v_tensor_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_1x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_1x1v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_1x1v_tensor_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_1x1v_tensor_p3(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_1x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_1x2v_tensor_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_1x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_1x2v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M2_1x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M2_1x2v_tensor_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M2_1x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M2_1x2v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_1x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_1x3v_tensor_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_1x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_1x3v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M2_1x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M2_1x3v_tensor_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M2_1x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M2_1x3v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_2x1v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_2x1v_tensor_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_2x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_2x1v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_2x1v_tensor_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_2x1v_tensor_p3(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_2x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_2x2v_tensor_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_2x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_2x2v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_2x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_2x3v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

#ifdef GKYL_BUILD_VLASOV_TEN_2X3V_P2
GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_2x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_2x3v_tensor_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}
#endif

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M2_3x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M2_3x3v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M2_2x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M2_2x2v_tensor_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M2_2x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M2_2x2v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M2_2x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M2_2x3v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

#ifdef GKYL_BUILD_VLASOV_TEN_2X3V_P2
GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M2_2x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M2_2x3v_tensor_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}
#endif

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M2_3x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M2_3x3v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

// M2 (H moment, velocity-space Hamiltonian) kernel list (Tensor basis)
GKYL_CU_D static const gkyl_vlasov_mom_kern_list tensor_hamil_vel_dense_m2_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M2_1x1v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_M2_1x1v_tensor_p2,
   kernel_mom_vlasov_hamil_vel_dense_M2_1x1v_tensor_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M2_1x2v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_M2_1x2v_tensor_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M2_1x3v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_M2_1x3v_tensor_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M2_2x1v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_M2_2x1v_tensor_p2,
   kernel_mom_vlasov_hamil_vel_dense_M2_2x1v_tensor_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M2_2x2v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_M2_2x2v_tensor_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M2_2x3v_tensor_p1, GKYL_TEN_2X3V_P2(kernel_mom_vlasov_hamil_vel_dense_M2_2x3v_tensor_p2), NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M2_3x3v_tensor_p1, NULL, NULL} // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D static const gkyl_vlasov_mom_kern_list tensor_hamil_vel_sparse_m2_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M2_1x1v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_M2_1x1v_tensor_p2,
   kernel_mom_vlasov_hamil_vel_dense_M2_1x1v_tensor_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M2_1x2v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_sparse_M2_1x2v_tensor_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M2_1x3v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_sparse_M2_1x3v_tensor_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M2_2x1v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_M2_2x1v_tensor_p2,
   kernel_mom_vlasov_hamil_vel_dense_M2_2x1v_tensor_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M2_2x2v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_sparse_M2_2x2v_tensor_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M2_2x3v_tensor_p1, GKYL_TEN_2X3V_P2(kernel_mom_vlasov_hamil_vel_sparse_M2_2x3v_tensor_p2), NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M2_3x3v_tensor_p1, NULL, NULL} // 6
};

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_1x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_1x1v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_1x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_1x1v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_1x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_1x1v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_1x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_1x1v_ser_p3(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_1x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_1x2v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_1x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_1x2v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_1x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_1x2v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_1x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_1x3v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_1x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_1x3v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_1x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_1x3v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_2x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_2x1v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_2x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_2x1v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_2x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_2x1v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_2x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_2x1v_ser_p3(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_2x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_2x2v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_2x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_2x2v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_2x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_2x2v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_2x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_2x3v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_2x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_2x3v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_3x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_3x3v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_M2_3x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_M2_3x3v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

// M2 (H moment, phase-space Hamiltonian) kernel list (Serendipity basis)
GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_hamil_phase_m2_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_M2_1x1v_ser_p1, kernel_mom_vlasov_hamil_phase_M2_1x1v_ser_p2,
   kernel_mom_vlasov_hamil_phase_M2_1x1v_ser_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_phase_M2_1x2v_ser_p1, kernel_mom_vlasov_hamil_phase_M2_1x2v_ser_p2,
   NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_phase_M2_1x3v_ser_p1, kernel_mom_vlasov_hamil_phase_M2_1x3v_ser_p2,
   NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_M2_2x1v_ser_p1, kernel_mom_vlasov_hamil_phase_M2_2x1v_ser_p2,
   kernel_mom_vlasov_hamil_phase_M2_2x1v_ser_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_phase_M2_2x2v_ser_p1, kernel_mom_vlasov_hamil_phase_M2_2x2v_ser_p2,
   NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_phase_M2_2x3v_ser_p1, NULL, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_M2_3x3v_ser_p1, NULL, NULL} // 6
};

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M3i_1x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M3i_1x1v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M3i_1x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M3i_1x1v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M3i_1x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M3i_1x1v_ser_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M3i_1x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M3i_1x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M3i_1x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M3i_1x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M3i_1x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M3i_1x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M3i_1x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M3i_1x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M3i_1x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M3i_1x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M3i_1x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M3i_1x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M3i_1x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M3i_1x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M3i_1x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M3i_1x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M3i_2x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M3i_2x1v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M3i_2x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M3i_2x1v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M3i_2x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M3i_2x1v_ser_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M3i_2x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M3i_2x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M3i_2x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M3i_2x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M3i_2x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M3i_2x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M3i_2x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M3i_2x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M3i_2x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M3i_2x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M3i_2x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M3i_2x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M3i_2x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M3i_2x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M3i_2x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M3i_2x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_M3i_3x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_M3i_3x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_M3i_3x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_M3i_3x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

// M3i (H*dH/dv moment, velocity-space Hamiltonian) kernel list (Serendipity basis)
GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_hamil_vel_dense_m3i_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M3i_1x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M3i_1x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_M3i_1x1v_ser_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M3i_1x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M3i_1x2v_ser_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M3i_1x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M3i_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M3i_2x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M3i_2x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_M3i_2x1v_ser_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M3i_2x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M3i_2x2v_ser_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M3i_2x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M3i_2x3v_ser_p2, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M3i_3x3v_ser_p1, NULL, NULL} // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_hamil_vel_sparse_m3i_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M3i_1x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M3i_1x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_M3i_1x1v_ser_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M3i_1x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_M3i_1x2v_ser_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M3i_1x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_M3i_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_M3i_2x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_M3i_2x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_M3i_2x1v_ser_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M3i_2x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_M3i_2x2v_ser_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M3i_2x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_M3i_2x3v_ser_p2, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_M3i_3x3v_ser_p1, NULL, NULL} // 6
};

// M3i (H*dH/dv moment, velocity-space Hamiltonian) kernel list (Tensor basis)
// Not yet implemented
GKYL_CU_D static const gkyl_vlasov_mom_kern_list tensor_hamil_vel_dense_m3i_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, NULL, NULL, NULL} // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D static const gkyl_vlasov_mom_kern_list tensor_hamil_vel_sparse_m3i_kernels[] = {
  // 1x kernels
  {NULL, NULL, NULL, NULL}, // 0
  {NULL, NULL, NULL, NULL}, // 1
  {NULL, NULL, NULL, NULL}, // 2
  // 2x kernels
  {NULL, NULL, NULL, NULL}, // 3
  {NULL, NULL, NULL, NULL}, // 4
  {NULL, NULL, NULL, NULL}, // 5
  // 3x kernels
  {NULL, NULL, NULL, NULL} // 6
};

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_1x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_1x1v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_1x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_1x1v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_1x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_1x1v_ser_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_1x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_1x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_five_moments_1x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_five_moments_1x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_1x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_1x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_five_moments_1x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_five_moments_1x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_1x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_1x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_five_moments_1x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_five_moments_1x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_1x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_1x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_five_moments_1x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_five_moments_1x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_2x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_2x1v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_2x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_2x1v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_2x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_2x1v_ser_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_2x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_2x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_five_moments_2x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_five_moments_2x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_2x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_2x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_five_moments_2x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_five_moments_2x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_2x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_2x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_five_moments_2x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_five_moments_2x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_2x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_2x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_five_moments_2x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_five_moments_2x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_3x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_3x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_five_moments_3x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_five_moments_3x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

// Five moments (Zeroth, First, and Second moment together) kernel list (Serendipity basis)
// {1, dH/dv, H} moments for velocity-space Hamiltonian.
GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_hamil_vel_dense_five_moments_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_five_moments_1x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_1x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_1x1v_ser_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_vel_dense_five_moments_1x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_1x2v_ser_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_vel_dense_five_moments_1x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_five_moments_2x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_2x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_2x1v_ser_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_vel_dense_five_moments_2x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_2x2v_ser_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_vel_dense_five_moments_2x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_2x3v_ser_p2, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_five_moments_3x3v_ser_p1, NULL, NULL} // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_hamil_vel_sparse_five_moments_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_five_moments_1x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_1x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_1x1v_ser_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_five_moments_1x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_five_moments_1x2v_ser_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_five_moments_1x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_five_moments_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_five_moments_2x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_2x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_2x1v_ser_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_five_moments_2x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_five_moments_2x2v_ser_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_five_moments_2x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_five_moments_2x3v_ser_p2, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_five_moments_3x3v_ser_p1, NULL, NULL} // 6
};

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_1x1v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_1x1v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_1x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_1x1v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_1x1v_tensor_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_1x1v_tensor_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_1x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_1x2v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_1x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_1x2v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_five_moments_1x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_five_moments_1x2v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_five_moments_1x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_five_moments_1x2v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_1x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_1x3v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_1x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_1x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_five_moments_1x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_five_moments_1x3v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_five_moments_1x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_five_moments_1x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_2x1v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_2x1v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_2x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_2x1v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_2x1v_tensor_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_2x1v_tensor_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_2x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_2x2v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_2x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_2x2v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_2x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_2x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

#ifdef GKYL_BUILD_VLASOV_TEN_2X3V_P2
GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_2x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_2x3v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}
#endif

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_five_moments_3x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_five_moments_3x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_five_moments_2x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_five_moments_2x2v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_five_moments_2x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_five_moments_2x2v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_five_moments_2x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_five_moments_2x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

#ifdef GKYL_BUILD_VLASOV_TEN_2X3V_P2
GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_five_moments_2x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_five_moments_2x3v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}
#endif

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_five_moments_3x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_five_moments_3x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

// Five moments (Zeroth, First, and Second moment together) kernel list (Tensor basis)
// {1, dH/dv, H} moments for velocity-space Hamiltonian.
GKYL_CU_D static const gkyl_vlasov_mom_kern_list tensor_hamil_vel_dense_five_moments_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_five_moments_1x1v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_1x1v_tensor_p2,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_1x1v_tensor_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_vel_dense_five_moments_1x2v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_1x2v_tensor_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_vel_dense_five_moments_1x3v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_1x3v_tensor_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_five_moments_2x1v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_2x1v_tensor_p2,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_2x1v_tensor_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_vel_dense_five_moments_2x2v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_2x2v_tensor_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_vel_dense_five_moments_2x3v_tensor_p1, GKYL_TEN_2X3V_P2(kernel_mom_vlasov_hamil_vel_dense_five_moments_2x3v_tensor_p2), NULL
  }, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_five_moments_3x3v_tensor_p1, NULL, NULL
  } // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D static const gkyl_vlasov_mom_kern_list tensor_hamil_vel_sparse_five_moments_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_five_moments_1x1v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_1x1v_tensor_p2,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_1x1v_tensor_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_five_moments_1x2v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_sparse_five_moments_1x2v_tensor_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_five_moments_1x3v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_sparse_five_moments_1x3v_tensor_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_five_moments_2x1v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_2x1v_tensor_p2,
   kernel_mom_vlasov_hamil_vel_dense_five_moments_2x1v_tensor_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_five_moments_2x2v_tensor_p1,
   kernel_mom_vlasov_hamil_vel_sparse_five_moments_2x2v_tensor_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_five_moments_2x3v_tensor_p1, GKYL_TEN_2X3V_P2(kernel_mom_vlasov_hamil_vel_sparse_five_moments_2x3v_tensor_p2), NULL
  }, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_five_moments_3x3v_tensor_p1, NULL, NULL
  } // 6
};

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_1x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_1x1v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_1x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_1x1v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_1x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_1x1v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_1x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_1x1v_ser_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_1x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_1x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_1x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_1x2v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_1x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_1x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_1x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_1x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_1x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_1x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_1x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_1x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_2x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_2x1v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_2x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_2x1v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_2x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_2x1v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_2x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_2x1v_ser_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_2x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_2x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_2x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_2x2v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_2x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_2x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_2x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_2x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_2x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_2x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_3x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_3x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_five_moments_3x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_five_moments_3x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

// Five moments (Zeroth, First, and Second moment together) kernel list (Serendipity basis)
// {1, dH/dv, H} moments for phase-space Hamiltonian.
GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_hamil_phase_five_moments_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_five_moments_1x1v_ser_p1,
   kernel_mom_vlasov_hamil_phase_five_moments_1x1v_ser_p2,
   kernel_mom_vlasov_hamil_phase_five_moments_1x1v_ser_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_phase_five_moments_1x2v_ser_p1,
   kernel_mom_vlasov_hamil_phase_five_moments_1x2v_ser_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_phase_five_moments_1x3v_ser_p1,
   kernel_mom_vlasov_hamil_phase_five_moments_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_five_moments_2x1v_ser_p1,
   kernel_mom_vlasov_hamil_phase_five_moments_2x1v_ser_p2,
   kernel_mom_vlasov_hamil_phase_five_moments_2x1v_ser_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_phase_five_moments_2x2v_ser_p1,
   kernel_mom_vlasov_hamil_phase_five_moments_2x2v_ser_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_phase_five_moments_2x3v_ser_p1, NULL, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_five_moments_3x3v_ser_p1, NULL, NULL} // 6
};

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_ser_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_1x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_1x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_int_five_moments_1x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_1x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_1x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_int_five_moments_1x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_1x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_1x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_int_five_moments_1x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_1x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_1x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_int_five_moments_1x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_ser_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_2x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_2x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_int_five_moments_2x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_2x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_2x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_int_five_moments_2x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_2x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_2x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_int_five_moments_2x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_2x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_2x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_int_five_moments_2x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_3x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_3x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_3x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_int_five_moments_3x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

// Integrated five moments (Zeroth, First, and Second moment together) kernel list (Serendipity basis)
// {1, dH/dv, H} moments for velocity-space Hamiltonian.
GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_hamil_vel_dense_int_five_moments_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_ser_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x2v_ser_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_ser_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x2v_ser_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x3v_ser_p2, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_int_five_moments_3x3v_ser_p1, NULL, NULL} // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_hamil_vel_sparse_int_five_moments_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_ser_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_1x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_1x2v_ser_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_1x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_ser_p1,
   kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_ser_p2,
   kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_ser_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_2x2v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_2x2v_ser_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_2x3v_ser_p1,
   kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_2x3v_ser_p2, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_3x3v_ser_p1, NULL, NULL} // 6
};

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_tensor_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_tensor_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_1x2v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_1x2v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_1x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_int_five_moments_1x2v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_1x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_int_five_moments_1x2v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_1x3v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_1x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_1x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_int_five_moments_1x3v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_1x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_int_five_moments_1x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_tensor_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_tensor_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_2x2v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_2x2v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_2x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

#ifdef GKYL_BUILD_VLASOV_TEN_2X3V_P2
GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_2x3v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}
#endif

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_dense_int_five_moments_3x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_dense_int_five_moments_3x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_2x2v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_int_five_moments_2x2v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_2x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_int_five_moments_2x2v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_2x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_int_five_moments_2x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

#ifdef GKYL_BUILD_VLASOV_TEN_2X3V_P2
GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_2x3v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_int_five_moments_2x3v_tensor_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}
#endif

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_3x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_vel_sparse_int_five_moments_3x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

// Integrated five moments (Zeroth, First, and Second moment together) kernel list (Tensor basis)
// {1, dH/dv, H} moments for velocity-space Hamiltonian.
GKYL_CU_D static const gkyl_vlasov_mom_kern_list tensor_hamil_vel_dense_int_five_moments_kernels[] =
  {
    // 1x kernels
    {NULL, kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_tensor_p1,
     kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_tensor_p2,
     kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_tensor_p3}, // 0
    {NULL, kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x2v_tensor_p1,
     kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x2v_tensor_p2, NULL}, // 1
    {NULL, kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x3v_tensor_p1,
     kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x3v_tensor_p2, NULL}, // 2
    // 2x kernels
    {NULL, kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_tensor_p1,
     kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_tensor_p2,
     kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_tensor_p3}, // 3
    {NULL, kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x2v_tensor_p1,
     kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x2v_tensor_p2, NULL}, // 4
    {NULL, kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x3v_tensor_p1, GKYL_TEN_2X3V_P2(kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x3v_tensor_p2),
     NULL}, // 5
    // 3x kernels
    {NULL, kernel_mom_vlasov_hamil_vel_dense_int_five_moments_3x3v_tensor_p1, NULL,
     NULL} // 6
};
// Sparse-Hamiltonian variant; currently identical to the dense table
// (points at the same kernels) until the sparse kernels land.
GKYL_CU_D static const gkyl_vlasov_mom_kern_list
  tensor_hamil_vel_sparse_int_five_moments_kernels[] = {
    // 1x kernels
    {NULL, kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_tensor_p1,
     kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_tensor_p2,
     kernel_mom_vlasov_hamil_vel_dense_int_five_moments_1x1v_tensor_p3}, // 0
    {NULL, kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_1x2v_tensor_p1,
     kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_1x2v_tensor_p2, NULL}, // 1
    {NULL, kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_1x3v_tensor_p1,
     kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_1x3v_tensor_p2, NULL}, // 2
    // 2x kernels
    {NULL, kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_tensor_p1,
     kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_tensor_p2,
     kernel_mom_vlasov_hamil_vel_dense_int_five_moments_2x1v_tensor_p3}, // 3
    {NULL, kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_2x2v_tensor_p1,
     kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_2x2v_tensor_p2, NULL}, // 4
    {NULL, kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_2x3v_tensor_p1, GKYL_TEN_2X3V_P2(kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_2x3v_tensor_p2),
     NULL}, // 5
    // 3x kernels
    {NULL, kernel_mom_vlasov_hamil_vel_sparse_int_five_moments_3x3v_tensor_p1, NULL,
     NULL} // 6
};

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_1x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_1x1v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_1x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_1x1v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_1x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_1x1v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_1x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_1x1v_ser_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_1x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_1x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_1x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_1x2v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_1x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_1x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_1x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_1x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_1x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_1x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_1x3v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_1x3v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_2x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_2x1v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_2x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_2x1v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_2x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_2x1v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_2x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_2x1v_ser_p3(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_2x2v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_2x2v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_2x2v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_2x2v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_2x2v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_2x2v_ser_p2(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_2x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_2x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_2x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_2x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_3x3v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_3x3v_ser_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_hamil_phase_int_five_moments_3x3v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;

  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  int idx_hamil[GKYL_MAX_DIM];
  for (int i = 0; i < mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset + i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);

  return mom_vlasov_hamil_phase_int_five_moments_3x3v_tensor_p1(
    xc, dx, idx,
    mom_vlasov->jacob_vel ? (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx) : 0,
    (const double *)gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out
  );
}

// Integrated five moments (Zeroth, First, and Second moment together) kernel list (Serendipity basis)
// {1, dH/dv, H} moments for phase-space Hamiltonian.
GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_hamil_phase_int_five_moments_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_int_five_moments_1x1v_ser_p1,
   kernel_mom_vlasov_hamil_phase_int_five_moments_1x1v_ser_p2,
   kernel_mom_vlasov_hamil_phase_int_five_moments_1x1v_ser_p3}, // 0
  {NULL, kernel_mom_vlasov_hamil_phase_int_five_moments_1x2v_ser_p1,
   kernel_mom_vlasov_hamil_phase_int_five_moments_1x2v_ser_p2, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_phase_int_five_moments_1x3v_ser_p1,
   kernel_mom_vlasov_hamil_phase_int_five_moments_1x3v_ser_p2, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_int_five_moments_2x1v_ser_p1,
   kernel_mom_vlasov_hamil_phase_int_five_moments_2x1v_ser_p2,
   kernel_mom_vlasov_hamil_phase_int_five_moments_2x1v_ser_p3}, // 3
  {NULL, kernel_mom_vlasov_hamil_phase_int_five_moments_2x2v_ser_p1,
   kernel_mom_vlasov_hamil_phase_int_five_moments_2x2v_ser_p2, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_phase_int_five_moments_2x3v_ser_p1, NULL, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_int_five_moments_3x3v_ser_p1, NULL, NULL} // 6
};

// Phase-space Hamiltonian moment kernels (Tensor basis): only the p=1
// tensor hybrid has a phase-space Hamiltonian representation.
GKYL_CU_D static const gkyl_vlasov_mom_kern_list tensor_hamil_phase_m1i_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_M1i_1x1v_tensor_p1, NULL, NULL}, // 0
  {NULL, kernel_mom_vlasov_hamil_phase_M1i_1x2v_tensor_p1, NULL, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_phase_M1i_1x3v_tensor_p1, NULL, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_M1i_2x1v_tensor_p1, NULL, NULL}, // 3
  {NULL, kernel_mom_vlasov_hamil_phase_M1i_2x2v_tensor_p1, NULL, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_phase_M1i_2x3v_tensor_p1, NULL, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_M1i_3x3v_tensor_p1, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_vlasov_mom_kern_list tensor_hamil_phase_m2_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_M2_1x1v_tensor_p1, NULL, NULL}, // 0
  {NULL, kernel_mom_vlasov_hamil_phase_M2_1x2v_tensor_p1, NULL, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_phase_M2_1x3v_tensor_p1, NULL, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_M2_2x1v_tensor_p1, NULL, NULL}, // 3
  {NULL, kernel_mom_vlasov_hamil_phase_M2_2x2v_tensor_p1, NULL, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_phase_M2_2x3v_tensor_p1, NULL, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_M2_3x3v_tensor_p1, NULL, NULL} // 6
};

GKYL_CU_D static const gkyl_vlasov_mom_kern_list tensor_hamil_phase_five_moments_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_five_moments_1x1v_tensor_p1, NULL, NULL}, // 0
  {NULL, kernel_mom_vlasov_hamil_phase_five_moments_1x2v_tensor_p1, NULL, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_phase_five_moments_1x3v_tensor_p1, NULL, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_five_moments_2x1v_tensor_p1, NULL, NULL}, // 3
  {NULL, kernel_mom_vlasov_hamil_phase_five_moments_2x2v_tensor_p1, NULL, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_phase_five_moments_2x3v_tensor_p1, NULL, NULL}, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_five_moments_3x3v_tensor_p1, NULL, NULL
  } // 6
};

GKYL_CU_D static const gkyl_vlasov_mom_kern_list tensor_hamil_phase_int_five_moments_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_int_five_moments_1x1v_tensor_p1, NULL, NULL}, // 0
  {NULL, kernel_mom_vlasov_hamil_phase_int_five_moments_1x2v_tensor_p1, NULL, NULL}, // 1
  {NULL, kernel_mom_vlasov_hamil_phase_int_five_moments_1x3v_tensor_p1, NULL, NULL}, // 2
  // 2x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_int_five_moments_2x1v_tensor_p1, NULL, NULL}, // 3
  {NULL, kernel_mom_vlasov_hamil_phase_int_five_moments_2x2v_tensor_p1, NULL, NULL}, // 4
  {NULL, kernel_mom_vlasov_hamil_phase_int_five_moments_2x3v_tensor_p1, NULL, NULL
  }, // 5
  // 3x kernels
  {NULL, kernel_mom_vlasov_hamil_phase_int_five_moments_3x3v_tensor_p1, NULL,
   NULL} // 6
};

// Special moment kernels for integrating the density over a subset of the domain.
// Only work in 1V right now for specific relativistic Vlasov applications JJ: 08/01/25
GKYL_CU_DH static void
kernel_mom_vlasov_M0_upper_1x1v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_upper_1x1v_tensor_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_upper_1x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_upper_1x1v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_upper_1x1v_tensor_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_upper_1x1v_tensor_p3(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_upper_2x1v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_upper_2x1v_tensor_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_upper_2x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_upper_2x1v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_upper_2x1v_tensor_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_upper_2x1v_tensor_p3(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_lower_1x1v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_lower_1x1v_tensor_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_lower_1x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_lower_1x1v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_lower_1x1v_tensor_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_lower_1x1v_tensor_p3(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_lower_2x1v_tensor_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_lower_2x1v_tensor_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_lower_2x1v_tensor_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_lower_2x1v_tensor_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_lower_2x1v_tensor_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_lower_2x1v_tensor_p3(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

// M0 upper half-plane kernel list

GKYL_CU_DH static void
kernel_mom_vlasov_M0_upper_1x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_upper_1x1v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_upper_1x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_upper_1x1v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_upper_1x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_upper_1x1v_ser_p3(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_upper_2x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_upper_2x1v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_upper_2x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_upper_2x1v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_upper_2x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_upper_2x1v_ser_p3(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_lower_1x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_lower_1x1v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_lower_1x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_lower_1x1v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_lower_1x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_lower_1x1v_ser_p3(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_lower_2x1v_ser_p2(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_lower_2x1v_ser_p2(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_lower_2x1v_ser_p1(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_lower_2x1v_ser_p1(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

GKYL_CU_DH static void
kernel_mom_vlasov_M0_lower_2x1v_ser_p3(
  const struct gkyl_mom_type *momt, const double *xc, const double *dx, const int *idx,
  const double *f, double *out, enum gkyl_vel_edge edge
)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);

  int cdim = mom_vlasov->momt.cdim;
  int pdim = mom_vlasov->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i = 0; i < pdim - cdim; ++i) {
    idx_vel[i] = idx[cdim + i];
  }
  long vidx = gkyl_range_idx(&mom_vlasov->vel_range, idx_vel);

  return mom_vlasov_M0_lower_2x1v_ser_p3(
    xc, dx, idx, (const double *)gkyl_array_cfetch(mom_vlasov->vmap, vidx),
    (const double *)gkyl_array_cfetch(mom_vlasov->jacob_vel, vidx), mom_vlasov->v_thresh,
    mom_vlasov_f_thresh(mom_vlasov, idx), f, out
  );
}

// Serendipity-basis variants of the threshold-density tables (M0 over a
// velocity sub-domain; 1V only).
GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_m0_upper_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_M0_upper_1x1v_ser_p1, kernel_mom_vlasov_M0_upper_1x1v_ser_p2,
   kernel_mom_vlasov_M0_upper_1x1v_ser_p3}, // 0
  // 2x kernels
  {NULL, kernel_mom_vlasov_M0_upper_2x1v_ser_p1, kernel_mom_vlasov_M0_upper_2x1v_ser_p2,
   kernel_mom_vlasov_M0_upper_2x1v_ser_p3}, // 1
  // 3x kernels
  {NULL, NULL, NULL, NULL} // 2
};

GKYL_CU_D static const gkyl_vlasov_mom_kern_list ser_m0_lower_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_M0_lower_1x1v_ser_p1, kernel_mom_vlasov_M0_lower_1x1v_ser_p2,
   kernel_mom_vlasov_M0_lower_1x1v_ser_p3}, // 0
  // 2x kernels
  {NULL, kernel_mom_vlasov_M0_lower_2x1v_ser_p1, kernel_mom_vlasov_M0_lower_2x1v_ser_p2,
   kernel_mom_vlasov_M0_lower_2x1v_ser_p3}, // 1
  // 3x kernels
  {NULL, NULL, NULL, NULL} // 2
};

GKYL_CU_D static const gkyl_vlasov_mom_kern_list tensor_m0_upper_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_M0_upper_1x1v_tensor_p1, kernel_mom_vlasov_M0_upper_1x1v_tensor_p2,
   kernel_mom_vlasov_M0_upper_1x1v_tensor_p3}, // 0
  // 2x kernels
  {NULL, kernel_mom_vlasov_M0_upper_2x1v_tensor_p1, kernel_mom_vlasov_M0_upper_2x1v_tensor_p2,
   kernel_mom_vlasov_M0_upper_2x1v_tensor_p3}, // 1
  // 3x kernels
  {NULL, NULL, NULL, NULL} // 2
};

// M0 lower half-plane kernel list
GKYL_CU_D static const gkyl_vlasov_mom_kern_list tensor_m0_lower_kernels[] = {
  // 1x kernels
  {NULL, kernel_mom_vlasov_M0_lower_1x1v_tensor_p1, kernel_mom_vlasov_M0_lower_1x1v_tensor_p2,
   kernel_mom_vlasov_M0_lower_1x1v_tensor_p3}, // 0
  // 2x kernels
  {NULL, kernel_mom_vlasov_M0_lower_2x1v_tensor_p1, kernel_mom_vlasov_M0_lower_2x1v_tensor_p2,
   kernel_mom_vlasov_M0_lower_2x1v_tensor_p3}, // 1
  // 3x kernels
  {NULL, NULL, NULL, NULL} // 2
};

/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_mom_vlasov_free(const struct gkyl_ref_count *ref);