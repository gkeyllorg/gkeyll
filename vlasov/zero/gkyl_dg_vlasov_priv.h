#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_vlasov_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

// Types for various kernels
typedef double (*vlasov_stream_surf_t)(const double *w, const double *dxv, 
  const double *alpha_surf_l, const double *alpha_surf_r, 
  const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
  const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef double (*vlasov_stream_boundary_surf_t)(const double *w, const double *dxv, 
  const double *alpha_surf_edge, const double *alpha_surf_skin, 
  const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
  const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
  const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out);

typedef double (*vlasov_accel_surf_t)(const double *w, const double *dxv,
  const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef double (*vlasov_accel_boundary_surf_t)(const double *w, const double *dxv,
  const double *field, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_vlasov_stream_vol_kern_list;
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_vlasov_vol_kern_list;

typedef struct { vlasov_stream_surf_t kernels[3]; } gkyl_dg_vlasov_stream_surf_kern_list;

typedef struct { vlasov_stream_boundary_surf_t kernels[3]; } gkyl_dg_vlasov_stream_boundary_surf_kern_list;

typedef struct { vlasov_accel_surf_t kernels[3]; } gkyl_dg_vlasov_accel_surf_kern_list;

typedef struct { vlasov_accel_boundary_surf_t kernels[3]; } gkyl_dg_vlasov_accel_boundary_surf_kern_list;

struct dg_vlasov {
  struct gkyl_dg_eqn eqn; // Base object.
  int cdim; // Config-space dimensions.
  int pdim; // Phase-space dimensions.
  vlasov_stream_surf_t stream_surf[3]; // Surface terms for streaming.
  vlasov_stream_boundary_surf_t stream_boundary_surf[3]; // Boundary surface terms for streaming
  vlasov_accel_surf_t accel_surf[3]; // Surface terms for acceleration.
  vlasov_accel_boundary_surf_t accel_boundary_surf[3]; // Surface terms for acceleration
  struct gkyl_range conf_range; // Configuration space range (for indexing fields)
  struct gkyl_range phase_range; // Phase space range (for indexing alpha_geo in geometry)
  struct gkyl_dg_vlasov_auxfields auxfields; // Auxiliary fields.
};

//
// Serendipity volume kernels (streaming only, no geometry)
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_1x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_1x1v_ser_p1(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_1x1v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_1x1v_ser_p2(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_1x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_1x2v_ser_p1(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_1x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_1x2v_ser_p2(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_1x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_1x3v_ser_p1(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_1x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_1x3v_ser_p2(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_2x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_2x2v_ser_p1(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_2x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_2x2v_ser_p2(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_2x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_2x3v_ser_p1(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_2x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_2x3v_ser_p2(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_3x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_3x3v_ser_p1(xc, dx, qIn, qRhsOut);
}

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_vlasov_stream_vol_kern_list ser_stream_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_stream_vol_1x1v_ser_p1, kernel_vlasov_stream_vol_1x1v_ser_p2 }, // 0
  { NULL, kernel_vlasov_stream_vol_1x2v_ser_p1, kernel_vlasov_stream_vol_1x2v_ser_p2 }, // 1
  { NULL, kernel_vlasov_stream_vol_1x3v_ser_p1, kernel_vlasov_stream_vol_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_stream_vol_2x2v_ser_p1, kernel_vlasov_stream_vol_2x2v_ser_p2 }, // 3
  { NULL, kernel_vlasov_stream_vol_2x3v_ser_p1, kernel_vlasov_stream_vol_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, kernel_vlasov_stream_vol_3x3v_ser_p1, NULL               }, // 5
};

//
// Tensor volume kernels (streaming only, no geometry)
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_1x1v_tensor_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_1x1v_tensor_p1(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_1x1v_tensor_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_1x1v_tensor_p2(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_1x2v_tensor_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_1x2v_tensor_p1(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_1x2v_tensor_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_1x2v_tensor_p2(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_1x3v_tensor_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_1x3v_tensor_p1(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_1x3v_tensor_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_1x3v_tensor_p2(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_2x2v_tensor_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_2x2v_tensor_p1(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_2x2v_tensor_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_2x2v_tensor_p2(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_2x3v_tensor_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_2x3v_tensor_p1(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_2x3v_tensor_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_2x3v_tensor_p2(xc, dx, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_stream_vol_3x3v_tensor_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  return vlasov_stream_vol_3x3v_tensor_p1(xc, dx, qIn, qRhsOut);
}

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_vlasov_stream_vol_kern_list tensor_stream_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_stream_vol_1x1v_tensor_p1, kernel_vlasov_stream_vol_1x1v_tensor_p2 }, // 0
  { NULL, kernel_vlasov_stream_vol_1x2v_tensor_p1, kernel_vlasov_stream_vol_1x2v_tensor_p2 }, // 1
  { NULL, kernel_vlasov_stream_vol_1x3v_tensor_p1, kernel_vlasov_stream_vol_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_stream_vol_2x2v_tensor_p1, kernel_vlasov_stream_vol_2x2v_tensor_p2 }, // 3
  { NULL, kernel_vlasov_stream_vol_2x3v_tensor_p1, kernel_vlasov_stream_vol_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, kernel_vlasov_stream_vol_3x3v_tensor_p1, NULL               }, // 5
};

//
// Serendipity volume kernels (full Vlasov-Maxwell)
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_vlasov_vol_1x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_1x1v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_1x1v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_1x1v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_1x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_1x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_1x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_1x2v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_1x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_1x3v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_1x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_1x3v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_2x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_2x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_2x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_2x2v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_2x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_2x3v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_2x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_2x3v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_3x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_3x3v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_D
static const gkyl_dg_vlasov_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_vol_1x1v_ser_p1, kernel_vlasov_vol_1x1v_ser_p2 }, // 0
  { NULL, kernel_vlasov_vol_1x2v_ser_p1, kernel_vlasov_vol_1x2v_ser_p2 }, // 1
  { NULL, kernel_vlasov_vol_1x3v_ser_p1, kernel_vlasov_vol_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_vol_2x2v_ser_p1, kernel_vlasov_vol_2x2v_ser_p2 }, // 3
  { NULL, kernel_vlasov_vol_2x3v_ser_p1, kernel_vlasov_vol_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, kernel_vlasov_vol_3x3v_ser_p1, NULL               }, // 5
};

//
// Tensor volume kernels (full Vlasov-Maxwell)
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_vlasov_vol_1x1v_tensor_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_1x1v_tensor_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_1x1v_tensor_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_1x1v_tensor_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_1x2v_tensor_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_1x2v_tensor_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_1x2v_tensor_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_1x2v_tensor_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_1x3v_tensor_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_1x3v_tensor_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_1x3v_tensor_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_1x3v_tensor_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_2x2v_tensor_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_2x2v_tensor_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_2x2v_tensor_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_2x2v_tensor_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_2x3v_tensor_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_2x3v_tensor_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_2x3v_tensor_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_2x3v_tensor_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_vol_3x3v_tensor_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_vol_3x3v_tensor_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx),
    qIn, qRhsOut);
}

GKYL_CU_D
static const gkyl_dg_vlasov_vol_kern_list tensor_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_vol_1x1v_tensor_p1, kernel_vlasov_vol_1x1v_tensor_p2 }, // 0
  { NULL, kernel_vlasov_vol_1x2v_tensor_p1, kernel_vlasov_vol_1x2v_tensor_p2 }, // 1
  { NULL, kernel_vlasov_vol_1x3v_tensor_p1, kernel_vlasov_vol_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_vol_2x2v_tensor_p1, kernel_vlasov_vol_2x2v_tensor_p2 }, // 3
  { NULL, kernel_vlasov_vol_2x3v_tensor_p1, kernel_vlasov_vol_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, kernel_vlasov_vol_3x3v_tensor_p1, NULL               }, // 5
};

// Streaming surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list ser_stream_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_surfx_1x1v_ser_p1, vlasov_surfx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_surfx_1x2v_ser_p1, vlasov_surfx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_surfx_1x3v_ser_p1, vlasov_surfx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_surfx_2x2v_ser_p1, vlasov_surfx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_surfx_2x3v_ser_p1, vlasov_surfx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surfx_3x3v_ser_p1, NULL                  }, // 5
};

// Streaming surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list tensor_stream_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_surfx_1x1v_tensor_p1, vlasov_surfx_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_surfx_1x2v_tensor_p1, vlasov_surfx_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_surfx_1x3v_tensor_p1, vlasov_surfx_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_surfx_2x2v_tensor_p1, vlasov_surfx_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_surfx_2x3v_tensor_p1, vlasov_surfx_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surfx_3x3v_tensor_p1, NULL                  }, // 5
};


// Streaming surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list ser_stream_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, vlasov_surfy_2x2v_ser_p1, vlasov_surfy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_surfy_2x3v_ser_p1, vlasov_surfy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surfy_3x3v_ser_p1, NULL                  }, // 5
};

// Streaming surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list tensor_stream_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, vlasov_surfy_2x2v_tensor_p1, vlasov_surfy_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_surfy_2x3v_tensor_p1, vlasov_surfy_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surfy_3x3v_tensor_p1, NULL                  }, // 5
};

// Streaming surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list ser_stream_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, vlasov_surfz_3x3v_ser_p1, NULL }, // 5
};

// Streaming surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list tensor_stream_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, vlasov_surfz_3x3v_tensor_p1, NULL }, // 5
};

// Acceleration (full Vlasov-Maxwell) surface kernel list: vx-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list ser_accel_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_surfvx_1x1v_ser_p1, vlasov_surfvx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_surfvx_1x2v_ser_p1, vlasov_surfvx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_surfvx_1x3v_ser_p1, vlasov_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_surfvx_2x2v_ser_p1, vlasov_surfvx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_surfvx_2x3v_ser_p1, vlasov_surfvx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surfvx_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration (full Vlasov-Maxwell) surface kernel list: vx-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list tensor_accel_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_surfvx_1x1v_tensor_p1, vlasov_surfvx_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_surfvx_1x2v_tensor_p1, vlasov_surfvx_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_surfvx_1x3v_tensor_p1, vlasov_surfvx_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_surfvx_2x2v_tensor_p1, vlasov_surfvx_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_surfvx_2x3v_tensor_p1, vlasov_surfvx_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surfvx_3x3v_tensor_p1, NULL                   }, // 5
};

// Acceleration (full Vlasov-Maxwell) surface kernel list: vy-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list ser_accel_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, vlasov_surfvy_1x2v_ser_p1, vlasov_surfvy_1x2v_ser_p2 }, // 1
  { NULL, vlasov_surfvy_1x3v_ser_p1, vlasov_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_surfvy_2x2v_ser_p1, vlasov_surfvy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_surfvy_2x3v_ser_p1, vlasov_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surfvy_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration (full Vlasov-Maxwell) surface kernel list: vy-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list tensor_accel_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, vlasov_surfvy_1x2v_tensor_p1, vlasov_surfvy_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_surfvy_1x3v_tensor_p1, vlasov_surfvy_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_surfvy_2x2v_tensor_p1, vlasov_surfvy_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_surfvy_2x3v_tensor_p1, vlasov_surfvy_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surfvy_3x3v_tensor_p1, NULL                   }, // 5
};

// Acceleration (full Vlasov-Maxwell) surface kernel list: vz-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list ser_accel_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_surfvz_1x3v_ser_p1, vlasov_surfvz_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, vlasov_surfvz_2x3v_ser_p1, vlasov_surfvz_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surfvz_3x3v_ser_p1, NULL }, // 5
};

// Acceleration (full Vlasov-Maxwell) surface kernel list: vz-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list tensor_accel_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_surfvz_1x3v_tensor_p1, vlasov_surfvz_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, vlasov_surfvz_2x3v_tensor_p1, vlasov_surfvz_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surfvz_3x3v_tensor_p1, NULL }, // 5
};

// Streaming boundary surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_kern_list ser_stream_boundary_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_boundary_surfx_1x1v_ser_p1, vlasov_boundary_surfx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_boundary_surfx_1x2v_ser_p1, vlasov_boundary_surfx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_boundary_surfx_1x3v_ser_p1, vlasov_boundary_surfx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_boundary_surfx_2x2v_ser_p1, vlasov_boundary_surfx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_boundary_surfx_2x3v_ser_p1, vlasov_boundary_surfx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_surfx_3x3v_ser_p1, NULL                   }, // 5
};

// Streaming boundary surface kernel list: x-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_kern_list tensor_stream_boundary_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_boundary_surfx_1x1v_tensor_p1, vlasov_boundary_surfx_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_boundary_surfx_1x2v_tensor_p1, vlasov_boundary_surfx_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_boundary_surfx_1x3v_tensor_p1, vlasov_boundary_surfx_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_boundary_surfx_2x2v_tensor_p1, vlasov_boundary_surfx_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_boundary_surfx_2x3v_tensor_p1, vlasov_boundary_surfx_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_surfx_3x3v_tensor_p1, NULL                   }, // 5
};

// Streaming boundary surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_kern_list ser_stream_boundary_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_boundary_surfy_2x2v_ser_p1, vlasov_boundary_surfy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_boundary_surfy_2x3v_ser_p1, vlasov_boundary_surfy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_surfy_3x3v_ser_p1, NULL                   }, // 5
};

// Streaming boundary surface kernel list: y-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_kern_list tensor_stream_boundary_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_boundary_surfy_2x2v_tensor_p1, vlasov_boundary_surfy_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_boundary_surfy_2x3v_tensor_p1, vlasov_boundary_surfy_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_surfy_3x3v_tensor_p1, NULL                   }, // 5
};

// Streaming boundary surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_kern_list ser_stream_boundary_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_surfz_3x3v_ser_p1, NULL }, // 5
};

// Streaming boundary surface kernel list: z-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_stream_boundary_surf_kern_list tensor_stream_boundary_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_surfz_3x3v_tensor_p1, NULL }, // 5
};

// Acceleration (full Vlasov-Maxwell) boundary surface kernel (zero-flux BCs) list: vx-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list ser_accel_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_boundary_surfvx_1x1v_ser_p1, vlasov_boundary_surfvx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_boundary_surfvx_1x2v_ser_p1, vlasov_boundary_surfvx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_boundary_surfvx_1x3v_ser_p1, vlasov_boundary_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_boundary_surfvx_2x2v_ser_p1, vlasov_boundary_surfvx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_boundary_surfvx_2x3v_ser_p1, vlasov_boundary_surfvx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_surfvx_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration (full Vlasov-Maxwell) boundary surface kernel (zero-flux BCs) list: vx-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list tensor_accel_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_boundary_surfvx_1x1v_tensor_p1, vlasov_boundary_surfvx_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_boundary_surfvx_1x2v_tensor_p1, vlasov_boundary_surfvx_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_boundary_surfvx_1x3v_tensor_p1, vlasov_boundary_surfvx_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_boundary_surfvx_2x2v_tensor_p1, vlasov_boundary_surfvx_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_boundary_surfvx_2x3v_tensor_p1, vlasov_boundary_surfvx_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_surfvx_3x3v_tensor_p1, NULL                   }, // 5
};

// Acceleration (full Vlasov-Maxwell) boundary surface kernel (zero-flux BCs) list: vy-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list ser_accel_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, vlasov_boundary_surfvy_1x2v_ser_p1, vlasov_boundary_surfvy_1x2v_ser_p2 }, // 1
  { NULL, vlasov_boundary_surfvy_1x3v_ser_p1, vlasov_boundary_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_boundary_surfvy_2x2v_ser_p1, vlasov_boundary_surfvy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_boundary_surfvy_2x3v_ser_p1, vlasov_boundary_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_surfvy_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration (full Vlasov-Maxwell) boundary surface kernel (zero-flux BCs) list: vy-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list tensor_accel_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, vlasov_boundary_surfvy_1x2v_tensor_p1, vlasov_boundary_surfvy_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_boundary_surfvy_1x3v_tensor_p1, vlasov_boundary_surfvy_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_boundary_surfvy_2x2v_tensor_p1, vlasov_boundary_surfvy_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_boundary_surfvy_2x3v_tensor_p1, vlasov_boundary_surfvy_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_surfvy_3x3v_tensor_p1, NULL                   }, // 5
};

// Acceleration (full Vlasov-Maxwell) boundary surface kernel (zero-flux BCs) list: vz-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list ser_accel_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_boundary_surfvz_1x3v_ser_p1, vlasov_boundary_surfvz_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, vlasov_boundary_surfvz_2x3v_ser_p1, vlasov_boundary_surfvz_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_surfvz_3x3v_ser_p1, NULL }, // 5
};

// Acceleration (full Vlasov-Maxwell) boundary surface kernel (zero-flux BCs) list: vz-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list tensor_accel_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_boundary_surfvz_1x3v_tensor_p1, vlasov_boundary_surfvz_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, vlasov_boundary_surfvz_2x3v_tensor_p1, vlasov_boundary_surfvz_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_surfvz_3x3v_tensor_p1, NULL }, // 5
};

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst,cdim,vd,poly_order) lst[cv_index[cdim].vdim[vd]].kernels[poly_order]

/**
 * Free vlasov eqn object.
 *
 * @param ref Reference counter for vlasov eqn
 */
void gkyl_vlasov_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double* xcL, const double* xcC, const double* xcR, 
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  if (dir < vlasov->cdim) {
    // Each cell owns the *lower* edge surface alpha
    // Since alpha is continuous, fetch alpha_surf in center cell for lower edge
    // and fetch alpha_surf in right cell for upper edge
    long pidxC = gkyl_range_idx(&vlasov->phase_range, idxC);
    long pidxR = gkyl_range_idx(&vlasov->phase_range, idxR);
    return vlasov->stream_surf[dir]
      (xcC, dxC,
       vlasov->auxfields.alpha_surf ? (const double*) gkyl_array_cfetch(vlasov->auxfields.alpha_surf, pidxC) : 0, 
       vlasov->auxfields.alpha_surf ? (const double*) gkyl_array_cfetch(vlasov->auxfields.alpha_surf, pidxR) : 0,
       vlasov->auxfields.sgn_alpha_surf ? (const double*) gkyl_array_cfetch(vlasov->auxfields.sgn_alpha_surf, pidxC) : 0, 
       vlasov->auxfields.sgn_alpha_surf ? (const double*) gkyl_array_cfetch(vlasov->auxfields.sgn_alpha_surf, pidxR) : 0,
       vlasov->auxfields.const_sgn_alpha ? (const int*) gkyl_array_cfetch(vlasov->auxfields.const_sgn_alpha, pidxC) : 0,
       vlasov->auxfields.const_sgn_alpha ? (const int*) gkyl_array_cfetch(vlasov->auxfields.const_sgn_alpha, pidxR) : 0,
       qInL, qInC, qInR, qRhsOut);
  }
  else {
    long cidx = gkyl_range_idx(&vlasov->conf_range, idxC);
    return vlasov->accel_surf[dir-vlasov->cdim]
      (xcC, dxC,
        vlasov->auxfields.field ? (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx) : 0,
        qInL, qInC, qInR, qRhsOut);
  }
}

GKYL_CU_D
static double
boundary_surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double* xcEdge, const double* xcSkin,
  const double* dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  if (dir < vlasov->cdim) {
    // Each cell owns the *lower* edge surface alpha
    long pidxEdge = gkyl_range_idx(&vlasov->phase_range, idxEdge);
    long pidxSkin = gkyl_range_idx(&vlasov->phase_range, idxSkin);
    return vlasov->stream_boundary_surf[dir]
      (xcSkin, dxSkin,
       vlasov->auxfields.alpha_surf ? (const double*) gkyl_array_cfetch(vlasov->auxfields.alpha_surf, pidxEdge) : 0,
       vlasov->auxfields.alpha_surf ? (const double*) gkyl_array_cfetch(vlasov->auxfields.alpha_surf, pidxSkin) : 0,
       vlasov->auxfields.sgn_alpha_surf ? (const double*) gkyl_array_cfetch(vlasov->auxfields.sgn_alpha_surf, pidxEdge) : 0,
       vlasov->auxfields.sgn_alpha_surf ? (const double*) gkyl_array_cfetch(vlasov->auxfields.sgn_alpha_surf, pidxSkin) : 0,
       vlasov->auxfields.const_sgn_alpha ? (const int*) gkyl_array_cfetch(vlasov->auxfields.const_sgn_alpha, pidxEdge) : 0,
       vlasov->auxfields.const_sgn_alpha ? (const int*) gkyl_array_cfetch(vlasov->auxfields.const_sgn_alpha, pidxSkin) : 0,
       edge, qInEdge, qInSkin, qRhsOut);
  } else if (dir >= vlasov->cdim) {
    long cidx = gkyl_range_idx(&vlasov->conf_range, idxSkin);
    return vlasov->accel_boundary_surf[dir-vlasov->cdim]
      (xcSkin, dxSkin,
        vlasov->auxfields.field ? (const double*) gkyl_array_cfetch(vlasov->auxfields.field, cidx) : 0,
        edge, qInEdge, qInSkin, qRhsOut);
  }
  return 0.;
}

GKYL_CU_D
static double
boundary_diag(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double* xcEdge, const double* xcSkin,
  const double* dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  return boundary_surf(eqn, dir, xcEdge, xcSkin, dxEdge, dxSkin, idxEdge, idxSkin, edge, qInEdge, qInSkin, qRhsOut);
}
