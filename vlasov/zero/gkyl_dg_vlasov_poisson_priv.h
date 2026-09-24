#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_vlasov_poisson_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

// Types for various kernels
typedef double (*vlasov_poisson_stream_surf_t)(const double *w, const double *dxv, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef double (*vlasov_poisson_stream_boundary_surf_t)(const double *w, const double *dxv, 
  const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out);

typedef double (*vlasov_poisson_accel_surf_t)(const double *w, const double *dxv,
  const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef double (*vlasov_poisson_accel_boundary_surf_t)(const double *w, const double *dxv,
  const double *pots, const double *EBext, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out);

// Null kernels used for vxB term in case without external fields.
GKYL_CU_DH static double kernel_vlasov_poisson_zero_accel_surf(const double *w, const double *dxv,
  const double *pots, const double *EBext, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out)
{
  return 0.0;
}
GKYL_CU_DH static double kernel_vlasov_poisson_zero_accel_boundary_surf(const double *w, const double *dxv,
  const double *pots, const double *EBext, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out)
{
  return 0.0;
}

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_vlasov_poisson_vol_kern_list;

typedef struct { vlasov_poisson_stream_surf_t kernels[3]; } gkyl_dg_vlasov_poisson_stream_surf_kern_list;
typedef struct { vlasov_poisson_accel_surf_t kernels[3]; } gkyl_dg_vlasov_poisson_accel_surf_kern_list;

typedef struct { vlasov_poisson_stream_boundary_surf_t kernels[3]; } gkyl_dg_vlasov_poisson_stream_boundary_surf_kern_list;
typedef struct { vlasov_poisson_accel_boundary_surf_t kernels[3]; } gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list;

struct dg_vlasov_poisson {
  struct gkyl_dg_eqn eqn; // Base object.
  int cdim; // Config-space dimensions.
  int pdim; // Phase-space dimensions.
  vlasov_poisson_stream_surf_t stream_surf[3]; // Surface terms for streaming.
  vlasov_poisson_stream_boundary_surf_t stream_boundary_surf[3]; // Boundary surface terms for streaming
  vlasov_poisson_accel_surf_t accel_surf[3]; // Surface terms for acceleration.
  vlasov_poisson_accel_boundary_surf_t accel_boundary_surf[3]; // Surface terms for acceleration
  struct gkyl_range conf_range; // Configuration space range (for indexing fields)
  struct gkyl_range phase_range; // Phase space range (for indexing alpha_geo in geometry)
  struct gkyl_dg_vlasov_poisson_auxfields auxfields; // Auxiliary fields.
};

//
// Serendipity streaming surface and boundary surface kernels.
//

// Streaming surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_stream_surf_kern_list ser_poisson_stream_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_poisson_surfx_1x1v_ser_p1, vlasov_poisson_surfx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_poisson_surfx_1x2v_ser_p1, vlasov_poisson_surfx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_poisson_surfx_1x3v_ser_p1, vlasov_poisson_surfx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_surfx_2x2v_ser_p1, vlasov_poisson_surfx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_surfx_2x3v_ser_p1, vlasov_poisson_surfx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                  }, // 5
};

// Streaming surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_stream_surf_kern_list ser_poisson_stream_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, vlasov_poisson_surfy_2x2v_ser_p1, vlasov_poisson_surfy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_surfy_2x3v_ser_p1, vlasov_poisson_surfy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                  }, // 5
};

// Streaming surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_stream_surf_kern_list ser_poisson_stream_surf_z_kernels[] = {
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

// Streaming boundary surface kernel list: x-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_stream_boundary_surf_kern_list ser_poisson_stream_boundary_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_poisson_boundary_surfx_1x1v_ser_p1, vlasov_poisson_boundary_surfx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_poisson_boundary_surfx_1x2v_ser_p1, vlasov_poisson_boundary_surfx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_poisson_boundary_surfx_1x3v_ser_p1, vlasov_poisson_boundary_surfx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_boundary_surfx_2x2v_ser_p1, vlasov_poisson_boundary_surfx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_boundary_surfx_2x3v_ser_p1, vlasov_poisson_boundary_surfx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};

// Streaming boundary surface kernel list: y-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_stream_boundary_surf_kern_list ser_poisson_stream_boundary_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_boundary_surfy_2x2v_ser_p1, vlasov_poisson_boundary_surfy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_boundary_surfy_2x3v_ser_p1, vlasov_poisson_boundary_surfy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};

// Streaming boundary surface kernel list: z-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_stream_boundary_surf_kern_list ser_poisson_stream_boundary_surf_z_kernels[] = {
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
// Serendipity volume kernels (Vlasov-Poisson)
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_vlasov_poisson_vol_1x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_vol_1x1v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_vol_1x1v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_vol_1x1v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_vol_1x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_vol_1x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_vol_1x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_vol_1x2v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_vol_1x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_vol_1x3v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_vol_1x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_vol_1x3v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_vol_2x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_vol_2x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_vol_2x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_vol_2x2v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_vol_2x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_vol_2x3v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_vol_2x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_vol_2x3v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

// Volume kernel list, phi only
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_vol_kern_list ser_poisson_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_poisson_vol_1x1v_ser_p1, kernel_vlasov_poisson_vol_1x1v_ser_p2 }, // 0
  { NULL, kernel_vlasov_poisson_vol_1x2v_ser_p1, kernel_vlasov_poisson_vol_1x2v_ser_p2 }, // 1
  { NULL, kernel_vlasov_poisson_vol_1x3v_ser_p1, kernel_vlasov_poisson_vol_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_poisson_vol_2x2v_ser_p1, kernel_vlasov_poisson_vol_2x2v_ser_p2 }, // 3
  { NULL, kernel_vlasov_poisson_vol_2x3v_ser_p1, kernel_vlasov_poisson_vol_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL               }, // 5
};

// Acceleration (phi only) surface kernel list: vx-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_surf_kern_list ser_poisson_accel_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_poisson_surfvx_1x1v_ser_p1, vlasov_poisson_surfvx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_poisson_surfvx_1x2v_ser_p1, vlasov_poisson_surfvx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_poisson_surfvx_1x3v_ser_p1, vlasov_poisson_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_surfvx_2x2v_ser_p1, vlasov_poisson_surfvx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_surfvx_2x3v_ser_p1, vlasov_poisson_surfvx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};

// Acceleration (phi only) surface kernel list: vy-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_surf_kern_list ser_poisson_accel_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, kernel_vlasov_poisson_zero_accel_surf, kernel_vlasov_poisson_zero_accel_surf }, // 1
  { NULL, kernel_vlasov_poisson_zero_accel_surf, kernel_vlasov_poisson_zero_accel_surf }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_surfvy_2x2v_ser_p1, vlasov_poisson_surfvy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_surfvy_2x3v_ser_p1, vlasov_poisson_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};

// Acceleration (phi only) surface kernel list: vz-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_surf_kern_list ser_poisson_accel_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, kernel_vlasov_poisson_zero_accel_surf, kernel_vlasov_poisson_zero_accel_surf }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, kernel_vlasov_poisson_zero_accel_surf, kernel_vlasov_poisson_zero_accel_surf }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// Acceleration (phi only) boundary surface kernel (zero-flux BCs) list: vx-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list ser_poisson_accel_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_poisson_boundary_surfvx_1x1v_ser_p1, vlasov_poisson_boundary_surfvx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_poisson_boundary_surfvx_1x2v_ser_p1, vlasov_poisson_boundary_surfvx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_poisson_boundary_surfvx_1x3v_ser_p1, vlasov_poisson_boundary_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_boundary_surfvx_2x2v_ser_p1, vlasov_poisson_boundary_surfvx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_boundary_surfvx_2x3v_ser_p1, vlasov_poisson_boundary_surfvx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};

// Acceleration (phi only) boundary surface kernel (zero-flux BCs) list: vy-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list ser_poisson_accel_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, kernel_vlasov_poisson_zero_accel_boundary_surf, kernel_vlasov_poisson_zero_accel_boundary_surf }, // 1
  { NULL, kernel_vlasov_poisson_zero_accel_boundary_surf, kernel_vlasov_poisson_zero_accel_boundary_surf }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_boundary_surfvy_2x2v_ser_p1, vlasov_poisson_boundary_surfvy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_boundary_surfvy_2x3v_ser_p1, vlasov_poisson_boundary_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};

// Acceleration (phi only) boundary surface kernel (zero-flux BCs) list: vz-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list ser_poisson_accel_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, kernel_vlasov_poisson_zero_accel_boundary_surf, kernel_vlasov_poisson_zero_accel_boundary_surf }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, kernel_vlasov_poisson_zero_accel_boundary_surf, kernel_vlasov_poisson_zero_accel_boundary_surf }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

//
// Serendipity volume kernels (Vlasov-Poisson + external phi and A)
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_phiA_vol_1x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_phiA_vol_1x1v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_phiA_vol_1x1v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_phiA_vol_1x1v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_phiA_vol_1x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_phiA_vol_1x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_phiA_vol_1x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_phiA_vol_1x2v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_phiA_vol_1x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_phiA_vol_1x3v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_phiA_vol_1x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_phiA_vol_1x3v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_phiA_vol_2x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_phiA_vol_2x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_phiA_vol_2x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_phiA_vol_2x2v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_phiA_vol_2x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_phiA_vol_2x3v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_phiA_vol_2x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_phiA_vol_2x3v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    0, qIn, qRhsOut);
}

// Volume kernel list (with external phi and A).
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_vol_kern_list ser_poisson_ext_phiA_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_poisson_ext_phiA_vol_1x1v_ser_p1, kernel_vlasov_poisson_ext_phiA_vol_1x1v_ser_p2 }, // 0
  { NULL, kernel_vlasov_poisson_ext_phiA_vol_1x2v_ser_p1, kernel_vlasov_poisson_ext_phiA_vol_1x2v_ser_p2 }, // 1
  { NULL, kernel_vlasov_poisson_ext_phiA_vol_1x3v_ser_p1, kernel_vlasov_poisson_ext_phiA_vol_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_poisson_ext_phiA_vol_2x2v_ser_p1, kernel_vlasov_poisson_ext_phiA_vol_2x2v_ser_p2 }, // 3
  { NULL, kernel_vlasov_poisson_ext_phiA_vol_2x3v_ser_p1, kernel_vlasov_poisson_ext_phiA_vol_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL               }, // 5
};

// Acceleration (phi and A) surface kernel list: vx-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_surf_kern_list ser_poisson_ext_phiA_accel_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_poisson_ext_phiA_surfvx_1x1v_ser_p1, vlasov_poisson_ext_phiA_surfvx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_poisson_ext_phiA_surfvx_1x2v_ser_p1, vlasov_poisson_ext_phiA_surfvx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_poisson_ext_phiA_surfvx_1x3v_ser_p1, vlasov_poisson_ext_phiA_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_ext_phiA_surfvx_2x2v_ser_p1, vlasov_poisson_ext_phiA_surfvx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_ext_phiA_surfvx_2x3v_ser_p1, vlasov_poisson_ext_phiA_surfvx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};

// Acceleration (phi and A) surface kernel list: vy-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_surf_kern_list ser_poisson_ext_phiA_accel_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, vlasov_poisson_ext_phiA_surfvy_1x2v_ser_p1, vlasov_poisson_ext_phiA_surfvy_1x2v_ser_p2 }, // 1
  { NULL, vlasov_poisson_ext_phiA_surfvy_1x3v_ser_p1, vlasov_poisson_ext_phiA_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_ext_phiA_surfvy_2x2v_ser_p1, vlasov_poisson_ext_phiA_surfvy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_ext_phiA_surfvy_2x3v_ser_p1, vlasov_poisson_ext_phiA_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};

// Acceleration (phi and A) surface kernel list: vz-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_surf_kern_list ser_poisson_ext_phiA_accel_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_poisson_ext_phiA_surfvz_1x3v_ser_p1, vlasov_poisson_ext_phiA_surfvz_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, vlasov_poisson_ext_phiA_surfvz_2x3v_ser_p1, vlasov_poisson_ext_phiA_surfvz_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// Acceleration (phi and A) boundary surface kernel (zero-flux BCs) list: vx-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list ser_poisson_ext_phiA_accel_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_poisson_ext_phiA_boundary_surfvx_1x1v_ser_p1, vlasov_poisson_ext_phiA_boundary_surfvx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_poisson_ext_phiA_boundary_surfvx_1x2v_ser_p1, vlasov_poisson_ext_phiA_boundary_surfvx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_poisson_ext_phiA_boundary_surfvx_1x3v_ser_p1, vlasov_poisson_ext_phiA_boundary_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_ext_phiA_boundary_surfvx_2x2v_ser_p1, vlasov_poisson_ext_phiA_boundary_surfvx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_ext_phiA_boundary_surfvx_2x3v_ser_p1, vlasov_poisson_ext_phiA_boundary_surfvx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};

// Acceleration (phi and A) boundary surface kernel (zero-flux BCs) list: vy-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list ser_poisson_ext_phiA_accel_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, vlasov_poisson_ext_phiA_boundary_surfvy_1x2v_ser_p1, vlasov_poisson_ext_phiA_boundary_surfvy_1x2v_ser_p2 }, // 1
  { NULL, vlasov_poisson_ext_phiA_boundary_surfvy_1x3v_ser_p1, vlasov_poisson_ext_phiA_boundary_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_ext_phiA_boundary_surfvy_2x2v_ser_p1, vlasov_poisson_ext_phiA_boundary_surfvy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_ext_phiA_boundary_surfvy_2x3v_ser_p1, vlasov_poisson_ext_phiA_boundary_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};
// Acceleration (phi and A) boundary surface kernel (zero-flux BCs) list: vz-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list ser_poisson_ext_phiA_accel_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_poisson_ext_phiA_boundary_surfvz_1x3v_ser_p1, vlasov_poisson_ext_phiA_boundary_surfvz_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, vlasov_poisson_ext_phiA_boundary_surfvz_2x3v_ser_p1, vlasov_poisson_ext_phiA_boundary_surfvz_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

//
// Serendipity volume kernels (Vlasov-Poisson + external E and B)
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_EB_vol_1x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_EB_vol_1x1v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    (const double*) gkyl_array_cfetch(vlasov->auxfields.fields_ext, cidx), 
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_EB_vol_1x1v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_EB_vol_1x1v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    (const double*) gkyl_array_cfetch(vlasov->auxfields.fields_ext, cidx), 
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_EB_vol_1x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_EB_vol_1x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    (const double*) gkyl_array_cfetch(vlasov->auxfields.fields_ext, cidx), 
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_EB_vol_1x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_EB_vol_1x2v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    (const double*) gkyl_array_cfetch(vlasov->auxfields.fields_ext, cidx), 
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_EB_vol_1x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_EB_vol_1x3v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    (const double*) gkyl_array_cfetch(vlasov->auxfields.fields_ext, cidx), 
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_EB_vol_1x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_EB_vol_1x3v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    (const double*) gkyl_array_cfetch(vlasov->auxfields.fields_ext, cidx), 
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_EB_vol_2x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_EB_vol_2x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    (const double*) gkyl_array_cfetch(vlasov->auxfields.fields_ext, cidx), 
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_EB_vol_2x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_EB_vol_2x2v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    (const double*) gkyl_array_cfetch(vlasov->auxfields.fields_ext, cidx), 
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_EB_vol_2x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_EB_vol_2x3v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    (const double*) gkyl_array_cfetch(vlasov->auxfields.fields_ext, cidx), 
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_vlasov_poisson_ext_EB_vol_2x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov_poisson_ext_EB_vol_2x3v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx), 
    (const double*) gkyl_array_cfetch(vlasov->auxfields.fields_ext, cidx), 
    qIn, qRhsOut);
}

// Volume kernel list (with external E and B).
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_vol_kern_list ser_poisson_ext_EB_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_poisson_ext_EB_vol_1x1v_ser_p1, kernel_vlasov_poisson_ext_EB_vol_1x1v_ser_p2 }, // 0
  { NULL, kernel_vlasov_poisson_ext_EB_vol_1x2v_ser_p1, kernel_vlasov_poisson_ext_EB_vol_1x2v_ser_p2 }, // 1
  { NULL, kernel_vlasov_poisson_ext_EB_vol_1x3v_ser_p1, kernel_vlasov_poisson_ext_EB_vol_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_poisson_ext_EB_vol_2x2v_ser_p1, kernel_vlasov_poisson_ext_EB_vol_2x2v_ser_p2 }, // 3
  { NULL, kernel_vlasov_poisson_ext_EB_vol_2x3v_ser_p1, kernel_vlasov_poisson_ext_EB_vol_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL               }, // 5
};

// Acceleration (E and B) surface kernel list: vx-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_surf_kern_list ser_poisson_ext_EB_accel_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_poisson_ext_EB_surfvx_1x1v_ser_p1, vlasov_poisson_ext_EB_surfvx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_poisson_ext_EB_surfvx_1x2v_ser_p1, vlasov_poisson_ext_EB_surfvx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_poisson_ext_EB_surfvx_1x3v_ser_p1, vlasov_poisson_ext_EB_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_ext_EB_surfvx_2x2v_ser_p1, vlasov_poisson_ext_EB_surfvx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_ext_EB_surfvx_2x3v_ser_p1, vlasov_poisson_ext_EB_surfvx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};

// Acceleration (E and B) surface kernel list: vy-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_surf_kern_list ser_poisson_ext_EB_accel_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, vlasov_poisson_ext_EB_surfvy_1x2v_ser_p1, vlasov_poisson_ext_EB_surfvy_1x2v_ser_p2 }, // 1
  { NULL, vlasov_poisson_ext_EB_surfvy_1x3v_ser_p1, vlasov_poisson_ext_EB_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_ext_EB_surfvy_2x2v_ser_p1, vlasov_poisson_ext_EB_surfvy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_ext_EB_surfvy_2x3v_ser_p1, vlasov_poisson_ext_EB_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};

// Acceleration (E and B) surface kernel list: vz-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_surf_kern_list ser_poisson_ext_EB_accel_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_poisson_ext_EB_surfvz_1x3v_ser_p1, vlasov_poisson_ext_EB_surfvz_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, vlasov_poisson_ext_EB_surfvz_2x3v_ser_p1, vlasov_poisson_ext_EB_surfvz_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// Acceleration (E and B) boundary surface kernel (zero-flux BCs) list: vx-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list ser_poisson_ext_EB_accel_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_poisson_ext_EB_boundary_surfvx_1x1v_ser_p1, vlasov_poisson_ext_EB_boundary_surfvx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_poisson_ext_EB_boundary_surfvx_1x2v_ser_p1, vlasov_poisson_ext_EB_boundary_surfvx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_poisson_ext_EB_boundary_surfvx_1x3v_ser_p1, vlasov_poisson_ext_EB_boundary_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_ext_EB_boundary_surfvx_2x2v_ser_p1, vlasov_poisson_ext_EB_boundary_surfvx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_ext_EB_boundary_surfvx_2x3v_ser_p1, vlasov_poisson_ext_EB_boundary_surfvx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};

// Acceleration (E and B) boundary surface kernel (zero-flux BCs) list: vy-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list ser_poisson_ext_EB_accel_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, vlasov_poisson_ext_EB_boundary_surfvy_1x2v_ser_p1, vlasov_poisson_ext_EB_boundary_surfvy_1x2v_ser_p2 }, // 1
  { NULL, vlasov_poisson_ext_EB_boundary_surfvy_1x3v_ser_p1, vlasov_poisson_ext_EB_boundary_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_ext_EB_boundary_surfvy_2x2v_ser_p1, vlasov_poisson_ext_EB_boundary_surfvy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_ext_EB_boundary_surfvy_2x3v_ser_p1, vlasov_poisson_ext_EB_boundary_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                   }, // 5
};
// Acceleration (E and B) boundary surface kernel (zero-flux BCs) list: vz-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list ser_poisson_ext_EB_accel_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_poisson_ext_EB_boundary_surfvz_1x3v_ser_p1, vlasov_poisson_ext_EB_boundary_surfvz_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, vlasov_poisson_ext_EB_boundary_surfvz_2x3v_ser_p1, vlasov_poisson_ext_EB_boundary_surfvz_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst,cdim,vd,poly_order) lst[cv_index[cdim].vdim[vd]].kernels[poly_order]

/**
 * Free vlasov eqn object.
 *
 * @param ref Reference counter for vlasov eqn
 */
void gkyl_vlasov_poisson_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double* xcL, const double* xcC, const double* xcR, 
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  if (dir < vlasov->cdim) {
    return vlasov->stream_surf[dir](xcC, dxC,
      qInL, qInC, qInR, qRhsOut);
  }
  else {
    long cidx = gkyl_range_idx(&vlasov->conf_range, idxC);
    return vlasov->accel_surf[dir-vlasov->cdim](xcC, dxC,
        (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx),
        (const double*) gkyl_array_cfetch(vlasov->auxfields.fields_ext, cidx), 
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
  struct dg_vlasov_poisson *vlasov = container_of(eqn, struct dg_vlasov_poisson, eqn);

  if (dir < vlasov->cdim) {
    return vlasov->stream_boundary_surf[dir](xcSkin, dxSkin,
       edge, qInEdge, qInSkin, qRhsOut);
  }
  else {
    long cidx = gkyl_range_idx(&vlasov->conf_range, idxSkin);
    return vlasov->accel_boundary_surf[dir-vlasov->cdim](xcSkin, dxSkin,
        (const double*) gkyl_array_cfetch(vlasov->auxfields.potentials, cidx),
        (const double*) gkyl_array_cfetch(vlasov->auxfields.fields_ext, cidx), 
        edge, qInEdge, qInSkin, qRhsOut);
  }
  return 0.;
}

#ifdef GKYL_HAVE_CUDA
/**
 * Create a new Vlasov equation object that lives on NV-GPU
 *
 * @param cbasis Configuration space basis functions.
 * @param pbasis Phase-space basis functions.
 * @param conf_range Configuration space range for use in indexing field.
 * @param phase_range Phase space range.
 * @param model_id enum to determine what type of Vlasov model.
 * @param field_id enum to determine what type of fields (e.g. phi or phi and A_ext).
 * @return Pointer to Vlasov equation object
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_poisson_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, const struct gkyl_range* phase_range,
  enum gkyl_model_id model_id, enum gkyl_field_id field_id);

/**
 * CUDA device function to set auxiliary fields (e.g. q/m*EM) needed in updating the force terms.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_vlasov_poisson_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_vlasov_poisson_auxfields auxin);
#endif
