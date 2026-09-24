#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_canonical_pb_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

// Types for various kernels
typedef double (*canonical_pb_stream_surf_t)(const double *w, const double *dxv,
  const double *hamil, 
  const double *alpha_surf_edge, const double *alpha_surf_skin, 
  const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
  const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef double (*canonical_pb_accel_surf_t)(const double *w, const double *dxv,
  const double *hamil,
  const double *alpha_surf_l, const double *alpha_surf_r, 
  const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
  const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef double (*canonical_pb_accel_boundary_surf_t)(const double *w, const double *dxv,
  const double *hamil, 
  const double *alpha_surf_edge, const double *alpha_surf_skin, 
  const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
  const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);

typedef double (*canonical_pb_stream_boundary_surf_t)(const double *w, const double *dxv,
  const double *hamil, 
  const double *alpha_surf_edge, const double *alpha_surf_skin, 
  const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
  const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);
  

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_canonical_pb_vol_kern_list;
typedef struct { canonical_pb_stream_surf_t kernels[3]; } gkyl_dg_canonical_pb_stream_surf_kern_list;
typedef struct { canonical_pb_accel_surf_t kernels[3]; } gkyl_dg_canonical_pb_accel_surf_kern_list;
typedef struct { canonical_pb_stream_boundary_surf_t kernels[3]; } gkyl_dg_canonical_pb_stream_boundary_surf_kern_list;
typedef struct { canonical_pb_accel_boundary_surf_t kernels[3]; } gkyl_dg_canonical_pb_accel_boundary_surf_kern_list;

struct dg_canonical_pb {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Config-space dimensions
  int pdim; // Phase-space dimensions
  canonical_pb_stream_surf_t stream_surf[3]; // Surface terms for streaming
  canonical_pb_accel_surf_t accel_surf[3]; // Surface terms for acceleration
  canonical_pb_stream_boundary_surf_t stream_boundary_surf[3]; // Surface terms for streaming
  canonical_pb_accel_boundary_surf_t accel_boundary_surf[3]; // Surface terms for acceleration
  struct gkyl_range conf_range; // configuration space range
  struct gkyl_range vel_range; // velocity space range
  struct gkyl_range phase_range; // velocity space range
  struct gkyl_dg_canonical_pb_auxfields auxfields; // Auxiliary fields.
};

// The cv_index[cd].vdim[cd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};


//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_canonical_pb_vol_1x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_1x1v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_1x1v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_1x1v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_1x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_1x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_1x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_1x2v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_1x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_1x3v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_1x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_1x3v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_2x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_2x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_2x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_2x2v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_2x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_2x3v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_2x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_2x3v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_canonical_pb_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_canonical_pb_vol_1x1v_ser_p1, kernel_canonical_pb_vol_1x1v_ser_p2 }, // 0 
  { NULL, kernel_canonical_pb_vol_1x2v_ser_p1, kernel_canonical_pb_vol_1x2v_ser_p2 }, // 1 
  { NULL, kernel_canonical_pb_vol_1x3v_ser_p1, kernel_canonical_pb_vol_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_canonical_pb_vol_2x2v_ser_p1, kernel_canonical_pb_vol_2x2v_ser_p2 }, // 3
  { NULL, kernel_canonical_pb_vol_2x3v_ser_p1, kernel_canonical_pb_vol_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL               }, // 5
};

//
// Serendipity surface kernels
//

// Streaming surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_stream_surf_kern_list ser_stream_surf_x_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_surfx_1x1v_ser_p1, canonical_pb_surfx_1x1v_ser_p2 }, // 0
  { NULL, canonical_pb_surfx_1x2v_ser_p1, canonical_pb_surfx_1x2v_ser_p2 }, // 1
  { NULL, canonical_pb_surfx_1x3v_ser_p1, canonical_pb_surfx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_surfx_2x2v_ser_p1, canonical_pb_surfx_2x2v_ser_p2 }, // 3
  { NULL, canonical_pb_surfx_2x3v_ser_p1, canonical_pb_surfx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                  }, // 5
};

// Streaming surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_stream_surf_kern_list ser_stream_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, canonical_pb_surfy_2x2v_ser_p1, canonical_pb_surfy_2x2v_ser_p2 }, // 3
  { NULL, canonical_pb_surfy_2x3v_ser_p1, canonical_pb_surfy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                  }, // 5
};

// Streaming surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_stream_surf_kern_list ser_stream_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// Acceleration surface kernel list: vx-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_accel_surf_kern_list ser_accel_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_surfvx_1x1v_ser_p1, canonical_pb_surfvx_1x1v_ser_p2 }, // 0
  { NULL, canonical_pb_surfvx_1x2v_ser_p1, canonical_pb_surfvx_1x2v_ser_p2 }, // 1
  { NULL, canonical_pb_surfvx_1x3v_ser_p1, canonical_pb_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_surfvx_2x2v_ser_p1, canonical_pb_surfvx_2x2v_ser_p2 }, // 3
  { NULL, canonical_pb_surfvx_2x3v_ser_p1, canonical_pb_surfvx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};

// Acceleration surface kernel list: vy-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_accel_surf_kern_list ser_accel_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, canonical_pb_surfvy_1x2v_ser_p1, canonical_pb_surfvy_1x2v_ser_p2 }, // 1
  { NULL, canonical_pb_surfvy_1x3v_ser_p1, canonical_pb_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_surfvy_2x2v_ser_p1, canonical_pb_surfvy_2x2v_ser_p2 }, // 3
  { NULL, canonical_pb_surfvy_2x3v_ser_p1, canonical_pb_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};

// Acceleration surface kernel list: vz-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_accel_surf_kern_list ser_accel_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, canonical_pb_surfvz_1x3v_ser_p1, canonical_pb_surfvz_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, canonical_pb_surfvz_2x3v_ser_p1, canonical_pb_surfvz_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vx-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_accel_boundary_surf_kern_list ser_accel_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_boundary_surfvx_1x1v_ser_p1, canonical_pb_boundary_surfvx_1x1v_ser_p2 }, // 0
  { NULL, canonical_pb_boundary_surfvx_1x2v_ser_p1, canonical_pb_boundary_surfvx_1x2v_ser_p2 }, // 1
  { NULL, canonical_pb_boundary_surfvx_1x3v_ser_p1, canonical_pb_boundary_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_boundary_surfvx_2x2v_ser_p1, canonical_pb_boundary_surfvx_2x2v_ser_p2 }, // 3
  { NULL, canonical_pb_boundary_surfvx_2x3v_ser_p1, canonical_pb_boundary_surfvx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vy-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_accel_boundary_surf_kern_list ser_accel_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, canonical_pb_boundary_surfvy_1x2v_ser_p1, canonical_pb_boundary_surfvy_1x2v_ser_p2 }, // 1
  { NULL, canonical_pb_boundary_surfvy_1x3v_ser_p1, canonical_pb_boundary_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_boundary_surfvy_2x2v_ser_p1, canonical_pb_boundary_surfvy_2x2v_ser_p2 }, // 3
  { NULL, canonical_pb_boundary_surfvy_2x3v_ser_p1, canonical_pb_boundary_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vz-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_accel_boundary_surf_kern_list ser_accel_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, canonical_pb_boundary_surfvz_1x3v_ser_p1, canonical_pb_boundary_surfvz_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, canonical_pb_boundary_surfvz_2x3v_ser_p1, canonical_pb_boundary_surfvz_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};


// Stream boundary surface kernel (zero-flux BCs) list: x-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_stream_boundary_surf_kern_list ser_stream_boundary_surf_x_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_boundary_surfx_1x1v_ser_p1, canonical_pb_boundary_surfx_1x1v_ser_p2 }, // 0
  { NULL, canonical_pb_boundary_surfx_1x2v_ser_p1, canonical_pb_boundary_surfx_1x2v_ser_p2 }, // 1
  { NULL, canonical_pb_boundary_surfx_1x3v_ser_p1, canonical_pb_boundary_surfx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_boundary_surfx_2x2v_ser_p1, canonical_pb_boundary_surfx_2x2v_ser_p2 }, // 3
  { NULL, canonical_pb_boundary_surfx_2x3v_ser_p1, canonical_pb_boundary_surfx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};

// Stream boundary surface kernel (zero-flux BCs) list: y-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_stream_boundary_surf_kern_list ser_stream_boundary_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, canonical_pb_boundary_surfy_2x2v_ser_p1, canonical_pb_boundary_surfy_2x2v_ser_p2 }, // 3
  { NULL, canonical_pb_boundary_surfy_2x3v_ser_p1, canonical_pb_boundary_surfy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};

// Stream boundary surface kernel (zero-flux BCs) list: z-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_stream_boundary_surf_kern_list ser_stream_boundary_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};


//
// Tensor volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_canonical_pb_vol_1x1v_tensor_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_1x1v_tensor_p1(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_1x1v_tensor_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_1x1v_tensor_p2(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_1x2v_tensor_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_1x2v_tensor_p1(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_1x2v_tensor_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_1x2v_tensor_p2(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_1x3v_tensor_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_1x3v_tensor_p1(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_1x3v_tensor_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_1x3v_tensor_p2(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_2x2v_tensor_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_2x2v_tensor_p1(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_2x2v_tensor_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_2x2v_tensor_p2(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}


GKYL_CU_DH
static double
kernel_canonical_pb_vol_2x3v_tensor_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_2x3v_tensor_p1(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_3x3v_tensor_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_3x3v_tensor_p1(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

// Volume kernel list for relativistic streaming + EM
GKYL_CU_D
static const gkyl_dg_canonical_pb_vol_kern_list tensor_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_canonical_pb_vol_1x1v_tensor_p1, kernel_canonical_pb_vol_1x1v_tensor_p2 }, // 0
  { NULL, kernel_canonical_pb_vol_1x2v_tensor_p1, kernel_canonical_pb_vol_1x2v_tensor_p2 }, // 1
  { NULL, kernel_canonical_pb_vol_1x3v_tensor_p1, kernel_canonical_pb_vol_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, kernel_canonical_pb_vol_2x2v_tensor_p1, kernel_canonical_pb_vol_2x2v_tensor_p2 }, // 3
  { NULL, kernel_canonical_pb_vol_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, kernel_canonical_pb_vol_3x3v_tensor_p1, NULL               }, // 2
};

//
// tensor surface kernels
//

// Streaming surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_stream_surf_kern_list tensor_stream_surf_x_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_surfx_1x1v_tensor_p1, canonical_pb_surfx_1x1v_tensor_p2 }, // 0
  { NULL, canonical_pb_surfx_1x2v_tensor_p1, canonical_pb_surfx_1x2v_tensor_p2 }, // 1
  { NULL, canonical_pb_surfx_1x3v_tensor_p1, canonical_pb_surfx_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_surfx_2x2v_tensor_p1, canonical_pb_surfx_2x2v_tensor_p2 }, // 3
  { NULL, canonical_pb_surfx_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_surfx_3x3v_tensor_p1, NULL }, // 5
};

// Streaming surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_stream_surf_kern_list tensor_stream_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, canonical_pb_surfy_2x2v_tensor_p1, canonical_pb_surfy_2x2v_tensor_p2 }, // 3
  { NULL, canonical_pb_surfy_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_surfy_3x3v_tensor_p1, NULL }, // 2
};

// Streaming surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_stream_surf_kern_list tensor_stream_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_surfz_3x3v_tensor_p1, NULL }, // 5
};

// Acceleration surface kernel list: vx-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_accel_surf_kern_list tensor_accel_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_surfvx_1x1v_tensor_p1, canonical_pb_surfvx_1x1v_tensor_p2 }, // 0
  { NULL, canonical_pb_surfvx_1x2v_tensor_p1, canonical_pb_surfvx_1x2v_tensor_p2 }, // 1
  { NULL, canonical_pb_surfvx_1x3v_tensor_p1, canonical_pb_surfvx_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_surfvx_2x2v_tensor_p1, canonical_pb_surfvx_2x2v_tensor_p2 }, // 3
  { NULL, canonical_pb_surfvx_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_surfvx_3x3v_tensor_p1, NULL }, // 5
};

// Acceleration surface kernel list: vy-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_accel_surf_kern_list tensor_accel_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, canonical_pb_surfvy_1x2v_tensor_p1, canonical_pb_surfvy_1x2v_tensor_p2 }, // 1
  { NULL, canonical_pb_surfvy_1x3v_tensor_p1, canonical_pb_surfvy_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_surfvy_2x2v_tensor_p1, canonical_pb_surfvy_2x2v_tensor_p2 }, // 3
  { NULL, canonical_pb_surfvy_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_surfvy_3x3v_tensor_p1, NULL }, // 5
};

// Acceleration surface kernel list: vz-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_accel_surf_kern_list tensor_accel_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, canonical_pb_surfvz_1x3v_tensor_p1, canonical_pb_surfvz_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, canonical_pb_surfvz_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_surfvz_3x3v_tensor_p1, NULL }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vx-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_accel_boundary_surf_kern_list tensor_accel_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_boundary_surfvx_1x1v_tensor_p1, canonical_pb_boundary_surfvx_1x1v_tensor_p2 }, // 0
  { NULL, canonical_pb_boundary_surfvx_1x2v_tensor_p1, canonical_pb_boundary_surfvx_1x2v_tensor_p2 }, // 1
  { NULL, canonical_pb_boundary_surfvx_1x3v_tensor_p1, canonical_pb_boundary_surfvx_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_boundary_surfvx_2x2v_tensor_p1, canonical_pb_boundary_surfvx_2x2v_tensor_p2 }, // 3
  { NULL, canonical_pb_boundary_surfvx_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_boundary_surfvx_3x3v_tensor_p1, NULL }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vy-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_accel_boundary_surf_kern_list tensor_accel_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, canonical_pb_boundary_surfvy_1x2v_tensor_p1, canonical_pb_boundary_surfvy_1x2v_tensor_p2  }, // 1
  { NULL, canonical_pb_boundary_surfvy_1x3v_tensor_p1, canonical_pb_boundary_surfvy_1x3v_tensor_p2  }, // 2
  // 2x kernels
  { NULL, canonical_pb_boundary_surfvy_2x2v_tensor_p1, canonical_pb_boundary_surfvy_2x2v_tensor_p2 }, // 3
  { NULL, canonical_pb_boundary_surfvy_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_boundary_surfvy_3x3v_tensor_p1, NULL }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vz-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_accel_boundary_surf_kern_list tensor_accel_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, canonical_pb_boundary_surfvz_1x3v_tensor_p1, canonical_pb_boundary_surfvz_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, canonical_pb_boundary_surfvz_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_boundary_surfvz_3x3v_tensor_p1, NULL }, // 5
};


// Stream boundary surface kernel (zero-flux BCs) list: x-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_stream_boundary_surf_kern_list tensor_stream_boundary_surf_x_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_boundary_surfx_1x1v_tensor_p1, canonical_pb_boundary_surfx_1x1v_tensor_p2 }, // 0
  { NULL, canonical_pb_boundary_surfx_1x2v_tensor_p1, canonical_pb_boundary_surfx_1x2v_tensor_p2 }, // 1
  { NULL, canonical_pb_boundary_surfx_1x3v_tensor_p1, canonical_pb_boundary_surfx_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, canonical_pb_boundary_surfx_2x2v_tensor_p1, canonical_pb_boundary_surfx_2x2v_tensor_p2 }, // 3
  { NULL, canonical_pb_boundary_surfx_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_boundary_surfx_3x3v_tensor_p1, NULL }, // 5
};

// Stream boundary surface kernel (zero-flux BCs) list: y-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_stream_boundary_surf_kern_list tensor_stream_boundary_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL  }, // 1
  { NULL, NULL, NULL  }, // 2
  // 2x kernels
  { NULL, canonical_pb_boundary_surfy_2x2v_tensor_p1, canonical_pb_boundary_surfy_2x2v_tensor_p2 }, // 3
  { NULL, canonical_pb_boundary_surfy_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_boundary_surfy_3x3v_tensor_p1, NULL }, // 5
};

// Stream boundary surface kernel (zero-flux BCs) list: z-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_stream_boundary_surf_kern_list tensor_stream_boundary_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, canonical_pb_boundary_surfz_3x3v_tensor_p1, NULL }, // 5
};

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst,cv_index,poly_order) lst[cv_index].kernels[poly_order]

/**
 * Free vlasov eqn object.
 *
 * @param ref Reference counter for vlasov eqn
 */
void gkyl_canonical_pb_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{    
  // Each cell owns the *lower* edge surface alpha
  // Since alpha is continuous, fetch alpha_surf in center cell for lower edge
  // and fetch alpha_surf in right cell for upper edge
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidxC = gkyl_range_idx(&canonical_pb->phase_range, idxC);
  long pidxR = gkyl_range_idx(&canonical_pb->phase_range, idxR);
  if (dir < canonical_pb->cdim) {
    return canonical_pb->stream_surf[dir](xcC, dxC,
        (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidxC),
        (const double*) gkyl_array_cfetch(canonical_pb->auxfields.alpha_surf, pidxC), 
        (const double*) gkyl_array_cfetch(canonical_pb->auxfields.alpha_surf, pidxR), 
        (const double*) gkyl_array_cfetch(canonical_pb->auxfields.sgn_alpha_surf, pidxC), 
        (const double*) gkyl_array_cfetch(canonical_pb->auxfields.sgn_alpha_surf, pidxR), 
        (const int*) gkyl_array_cfetch(canonical_pb->auxfields.const_sgn_alpha, pidxC), 
        (const int*) gkyl_array_cfetch(canonical_pb->auxfields.const_sgn_alpha, pidxR), 
        qInL, qInC, qInR, qRhsOut);    
  }
  else {
    return canonical_pb->accel_surf[dir-canonical_pb->cdim](xcC, dxC,
        (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidxC),
        (const double*) gkyl_array_cfetch(canonical_pb->auxfields.alpha_surf, pidxC), 
        (const double*) gkyl_array_cfetch(canonical_pb->auxfields.alpha_surf, pidxR), 
        (const double*) gkyl_array_cfetch(canonical_pb->auxfields.sgn_alpha_surf, pidxC), 
        (const double*) gkyl_array_cfetch(canonical_pb->auxfields.sgn_alpha_surf, pidxR), 
        (const int*) gkyl_array_cfetch(canonical_pb->auxfields.const_sgn_alpha, pidxC), 
        (const int*) gkyl_array_cfetch(canonical_pb->auxfields.const_sgn_alpha, pidxR), 
        qInL, qInC, qInR, qRhsOut);
  }
  return 0.;
}

GKYL_CU_D
static double
boundary_surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcEdge, const double*  xcSkin,
  const double*  dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  if (dir < canonical_pb->cdim) {
    // Each cell owns the *lower* edge surface alpha
    long pidxEdge = gkyl_range_idx(&canonical_pb->phase_range, idxEdge);
    long pidxSkin = gkyl_range_idx(&canonical_pb->phase_range, idxSkin);
    return canonical_pb->stream_boundary_surf[dir](xcSkin, dxSkin,
      (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidxSkin),
      (const double*) gkyl_array_cfetch(canonical_pb->auxfields.alpha_surf, pidxEdge), 
      (const double*) gkyl_array_cfetch(canonical_pb->auxfields.alpha_surf, pidxSkin), 
      (const double*) gkyl_array_cfetch(canonical_pb->auxfields.sgn_alpha_surf, pidxEdge), 
      (const double*) gkyl_array_cfetch(canonical_pb->auxfields.sgn_alpha_surf, pidxSkin), 
      (const int*) gkyl_array_cfetch(canonical_pb->auxfields.const_sgn_alpha, pidxEdge), 
      (const int*) gkyl_array_cfetch(canonical_pb->auxfields.const_sgn_alpha, pidxSkin), 
      edge, qInEdge, qInSkin, qRhsOut);
  }
  if (dir >= canonical_pb->cdim) {
    // Each cell owns the *lower* edge surface alpha
    long pidxEdge = gkyl_range_idx(&canonical_pb->phase_range, idxEdge);
    long pidxSkin = gkyl_range_idx(&canonical_pb->phase_range, idxSkin);
    return canonical_pb->accel_boundary_surf[dir-canonical_pb->cdim](xcSkin, dxSkin,
      (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidxSkin),
      (const double*) gkyl_array_cfetch(canonical_pb->auxfields.alpha_surf, pidxEdge), 
      (const double*) gkyl_array_cfetch(canonical_pb->auxfields.alpha_surf, pidxSkin), 
      (const double*) gkyl_array_cfetch(canonical_pb->auxfields.sgn_alpha_surf, pidxEdge), 
      (const double*) gkyl_array_cfetch(canonical_pb->auxfields.sgn_alpha_surf, pidxSkin), 
      (const int*) gkyl_array_cfetch(canonical_pb->auxfields.const_sgn_alpha, pidxEdge), 
      (const int*) gkyl_array_cfetch(canonical_pb->auxfields.const_sgn_alpha, pidxSkin), 
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

