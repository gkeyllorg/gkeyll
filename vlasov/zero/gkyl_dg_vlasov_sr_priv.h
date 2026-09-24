#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_vlasov_sr_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

// Types for various kernels
typedef double (*vlasov_sr_stream_surf_t)(const double *w, const double *dxv,
  const double *gamma, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef double (*vlasov_sr_accel_surf_t)(const double *w, const double *dxv,
  const double *gamma, const double *qmem, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef double (*vlasov_sr_accel_boundary_surf_t)(const double *w, const double *dxv,
  const double *gamma, const double *qmem, 
  const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_vlasov_sr_stream_vol_kern_list;
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_vlasov_sr_vol_kern_list;
typedef struct { vlasov_sr_stream_surf_t kernels[3]; } gkyl_dg_vlasov_sr_stream_surf_kern_list;
typedef struct { vlasov_sr_accel_surf_t kernels[3]; } gkyl_dg_vlasov_sr_accel_surf_kern_list;
typedef struct { vlasov_sr_accel_boundary_surf_t kernels[3]; } gkyl_dg_vlasov_sr_accel_boundary_surf_kern_list;

struct dg_vlasov_sr {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Config-space dimensions
  int pdim; // Phase-space dimensions
  vlasov_sr_stream_surf_t stream_surf[3]; // Surface terms for streaming
  vlasov_sr_accel_surf_t accel_surf[3]; // Surface terms for acceleration
  vlasov_sr_accel_boundary_surf_t accel_boundary_surf[3]; // Surface terms for acceleration
  struct gkyl_range conf_range; // configuration space range
  struct gkyl_range vel_range; // velocity space range
  struct gkyl_dg_vlasov_sr_auxfields auxfields; // Auxiliary fields.
};

//
// Serendipity volume kernels relativistic streaming only
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_vlasov_sr_stream_vol_1x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_stream_vol_1x1v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 0,
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_stream_vol_1x1v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_stream_vol_1x1v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 0,
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_stream_vol_1x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_stream_vol_1x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 0,
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_stream_vol_1x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_stream_vol_1x2v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 0,
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_stream_vol_1x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_stream_vol_1x3v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 0,
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_stream_vol_1x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_stream_vol_1x3v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 0,
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_stream_vol_2x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_stream_vol_2x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 0,
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_stream_vol_2x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_stream_vol_2x2v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 0,
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_stream_vol_2x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_stream_vol_2x3v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 0,
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_stream_vol_2x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_stream_vol_2x3v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 0,
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_stream_vol_3x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_stream_vol_3x3v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 0,
    qIn, qRhsOut);  
}

// Volume kernel list for relativistic streaming only
GKYL_CU_D
static const gkyl_dg_vlasov_sr_stream_vol_kern_list ser_stream_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_sr_stream_vol_1x1v_ser_p1, kernel_vlasov_sr_stream_vol_1x1v_ser_p2 }, // 0
  { NULL, kernel_vlasov_sr_stream_vol_1x2v_ser_p1, kernel_vlasov_sr_stream_vol_1x2v_ser_p2 }, // 1
  { NULL, kernel_vlasov_sr_stream_vol_1x3v_ser_p1, kernel_vlasov_sr_stream_vol_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_sr_stream_vol_2x2v_ser_p1, kernel_vlasov_sr_stream_vol_2x2v_ser_p2 }, // 3
  { NULL, kernel_vlasov_sr_stream_vol_2x3v_ser_p1, kernel_vlasov_sr_stream_vol_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, kernel_vlasov_sr_stream_vol_3x3v_ser_p1, NULL               }, // 5
};

//
// Serendipity volume kernels (streaming + EM) kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_vlasov_sr_vol_1x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_vol_1x1v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.qmem, cidx),
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_vol_1x1v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_vol_1x1v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.qmem, cidx),
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_vol_1x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_vol_1x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.qmem, cidx),
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_vol_1x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_vol_1x2v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.qmem, cidx),
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_vol_1x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_vol_1x3v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.qmem, cidx),
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_vol_1x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_vol_1x3v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.qmem, cidx),
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_vol_2x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_vol_2x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.qmem, cidx),
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_vol_2x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_vol_2x2v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.qmem, cidx),
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_vol_2x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_vol_2x3v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.qmem, cidx),
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_vol_2x3v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_vol_2x3v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.qmem, cidx),
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_vlasov_sr_vol_3x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idx[vlasov_sr->cdim+i];

  long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idx);
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  return vlasov_sr_vol_3x3v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 
    (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.qmem, cidx),
    qIn, qRhsOut);  
}

// Volume kernel list for relativistic streaming + EM
GKYL_CU_D
static const gkyl_dg_vlasov_sr_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_vlasov_sr_vol_1x1v_ser_p1, kernel_vlasov_sr_vol_1x1v_ser_p2 }, // 0
  { NULL, kernel_vlasov_sr_vol_1x2v_ser_p1, kernel_vlasov_sr_vol_1x2v_ser_p2 }, // 1
  { NULL, kernel_vlasov_sr_vol_1x3v_ser_p1, kernel_vlasov_sr_vol_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_vlasov_sr_vol_2x2v_ser_p1, kernel_vlasov_sr_vol_2x2v_ser_p2 }, // 3
  { NULL, kernel_vlasov_sr_vol_2x3v_ser_p1, kernel_vlasov_sr_vol_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, kernel_vlasov_sr_vol_3x3v_ser_p1, NULL               }, // 5
};

//
// Serendipity surface kernels
//

// Streaming surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_vlasov_sr_stream_surf_kern_list ser_stream_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_sr_surfx_1x1v_ser_p1, vlasov_sr_surfx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_sr_surfx_1x2v_ser_p1, vlasov_sr_surfx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_sr_surfx_1x3v_ser_p1, vlasov_sr_surfx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_sr_surfx_2x2v_ser_p1, vlasov_sr_surfx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_sr_surfx_2x3v_ser_p1, vlasov_sr_surfx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_sr_surfx_3x3v_ser_p1, NULL                  }, // 5
};

// Streaming surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_vlasov_sr_stream_surf_kern_list ser_stream_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, vlasov_sr_surfy_2x2v_ser_p1, vlasov_sr_surfy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_sr_surfy_2x3v_ser_p1, vlasov_sr_surfy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_sr_surfy_3x3v_ser_p1, NULL                  }, // 5
};

// Streaming surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_vlasov_sr_stream_surf_kern_list ser_stream_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, vlasov_sr_surfz_3x3v_ser_p1, NULL }, // 5
};

// Acceleration surface kernel list: vx-direction
GKYL_CU_D
static const gkyl_dg_vlasov_sr_accel_surf_kern_list ser_accel_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_sr_surfvx_1x1v_ser_p1, vlasov_sr_surfvx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_sr_surfvx_1x2v_ser_p1, vlasov_sr_surfvx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_sr_surfvx_1x3v_ser_p1, vlasov_sr_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_sr_surfvx_2x2v_ser_p1, vlasov_sr_surfvx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_sr_surfvx_2x3v_ser_p1, vlasov_sr_surfvx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_sr_surfvx_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration surface kernel list: vy-direction
GKYL_CU_D
static const gkyl_dg_vlasov_sr_accel_surf_kern_list ser_accel_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, vlasov_sr_surfvy_1x2v_ser_p1, vlasov_sr_surfvy_1x2v_ser_p2 }, // 1
  { NULL, vlasov_sr_surfvy_1x3v_ser_p1, vlasov_sr_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_sr_surfvy_2x2v_ser_p1, vlasov_sr_surfvy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_sr_surfvy_2x3v_ser_p1, vlasov_sr_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_sr_surfvy_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration surface kernel list: vz-direction
GKYL_CU_D
static const gkyl_dg_vlasov_sr_accel_surf_kern_list ser_accel_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_sr_surfvz_1x3v_ser_p1, vlasov_sr_surfvz_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, vlasov_sr_surfvz_2x3v_ser_p1, vlasov_sr_surfvz_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_sr_surfvz_3x3v_ser_p1, NULL }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vx-direction
GKYL_CU_D
static const gkyl_dg_vlasov_sr_accel_boundary_surf_kern_list ser_accel_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_sr_boundary_surfvx_1x1v_ser_p1, vlasov_sr_boundary_surfvx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_sr_boundary_surfvx_1x2v_ser_p1, vlasov_sr_boundary_surfvx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_sr_boundary_surfvx_1x3v_ser_p1, vlasov_sr_boundary_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_sr_boundary_surfvx_2x2v_ser_p1, vlasov_sr_boundary_surfvx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_sr_boundary_surfvx_2x3v_ser_p1, vlasov_sr_boundary_surfvx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_sr_boundary_surfvx_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vy-direction
GKYL_CU_D
static const gkyl_dg_vlasov_sr_accel_boundary_surf_kern_list ser_accel_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, vlasov_sr_boundary_surfvy_1x2v_ser_p1, vlasov_sr_boundary_surfvy_1x2v_ser_p2 }, // 1
  { NULL, vlasov_sr_boundary_surfvy_1x3v_ser_p1, vlasov_sr_boundary_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_sr_boundary_surfvy_2x2v_ser_p1, vlasov_sr_boundary_surfvy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_sr_boundary_surfvy_2x3v_ser_p1, vlasov_sr_boundary_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_sr_boundary_surfvy_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vz-direction
GKYL_CU_D
static const gkyl_dg_vlasov_sr_accel_boundary_surf_kern_list ser_accel_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_sr_boundary_surfvz_1x3v_ser_p1, vlasov_sr_boundary_surfvz_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, vlasov_sr_boundary_surfvz_2x3v_ser_p1, vlasov_sr_boundary_surfvz_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_sr_boundary_surfvz_3x3v_ser_p1, NULL }, // 5
};

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst,cdim,vd,poly_order) lst[cv_index[cdim].vdim[vd]].kernels[poly_order]

/**
 * Free vlasov eqn object.
 *
 * @param ref Reference counter for vlasov eqn
 */
void gkyl_vlasov_sr_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);

  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idxC[vlasov_sr->cdim+i];
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  if (dir < vlasov_sr->cdim) {
    return vlasov_sr->stream_surf[dir]
      (xcC, dxC, 
        (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx), 
        qInL, qInC, qInR, qRhsOut);
  }
  else {
    long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idxC);
    return vlasov_sr->accel_surf[dir-vlasov_sr->cdim]
      (xcC, dxC,
        (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx),
        vlasov_sr->auxfields.qmem ? (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.qmem, cidx) : 0,
        qInL, qInC, qInR, qRhsOut);
  }
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
  struct dg_vlasov_sr *vlasov_sr = container_of(eqn, struct dg_vlasov_sr, eqn);

  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<vlasov_sr->pdim-vlasov_sr->cdim; ++i)
    idx_vel[i] = idxSkin[vlasov_sr->cdim+i];
  long vidx = gkyl_range_idx(&vlasov_sr->vel_range, idx_vel);

  if (dir >= vlasov_sr->cdim) {
    long cidx = gkyl_range_idx(&vlasov_sr->conf_range, idxSkin);
    return vlasov_sr->accel_boundary_surf[dir-vlasov_sr->cdim]
      (xcSkin, dxSkin,
        (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.gamma, vidx),
        vlasov_sr->auxfields.qmem ? (const double*) gkyl_array_cfetch(vlasov_sr->auxfields.qmem, cidx) : 0,
        edge, qInEdge, qInSkin, qRhsOut);
  }
  return 0.;
}
