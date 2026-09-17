#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_ref_count.h>
#include <gkyl_mom_fpo_vlasov_kernels.h>

#include <assert.h>

struct mom_type_bcorr_fpo_vlasov {
  struct gkyl_mom_type momt;
  double vBoundary[6];
  struct gkyl_range phase_range; // velocity space range
  struct gkyl_mom_bcorr_fpo_vlasov_auxfields auxfields; // Auxiliary fields.
  bool use_gpu;
};

// for use in kernel tables
typedef struct {
  momf_t kernels[4];
} gkyl_mom_bcorr_fpo_vlasov_kern_list;

// Fill pdim phase index given edge information and idx array of remaining directions
GKYL_CU_DH
static void
set_phase_idx(enum gkyl_vel_edge edge, const struct gkyl_range *phase_range, int cdim, int pdim, const int *idx, int fidx[GKYL_MAX_DIM])
{
  int vdim = pdim - cdim;
  int vel_idx = cdim+edge%vdim;
  int cell = edge < vdim ? 1 : phase_range->upper[vel_idx];

  for (int i=0; i<(vel_idx); ++i) fidx[i] = idx[i];
  fidx[vel_idx] = cell;
  for (int i=vel_idx+1; i<pdim; ++i) fidx[i] = idx[i-1];
}

GKYL_CU_DH
static void
kernel_mom_bcorr_fpo_vlasov_1x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_bcorr_fpo_vlasov *mom_fpo_vlasov = container_of(momt, struct mom_type_bcorr_fpo_vlasov, momt);
  enum gkyl_vel_edge edge = *(enum gkyl_vel_edge *)param;

  int fidx[GKYL_MAX_DIM];
  double w[GKYL_MAX_DIM];

  int cdim = mom_fpo_vlasov->momt.cdim;
  int pdim = mom_fpo_vlasov->momt.pdim;

  gkyl_copy_int_arr(GKYL_MAX_DIM, idx, fidx);

  int vdim = 3;
  int vel_idx = edge%vdim;
  for (int i=0; i<cdim+vel_idx; ++i) w[i] = xc[i];
  for (int i=cdim+vel_idx; i<pdim-1; ++i) w[i] = xc[i+1];

  long linc = gkyl_range_idx(&mom_fpo_vlasov->phase_range, fidx);

  return mom_bcorr_fpo_vlasov_1x3v_ser_p1(w, fidx, edge, mom_fpo_vlasov->vBoundary, dx,
    (const double*) gkyl_array_cfetch(mom_fpo_vlasov->auxfields.D, linc), 
    f, out);
}

GKYL_CU_DH
static void
kernel_mom_bcorr_fpo_vlasov_1x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_bcorr_fpo_vlasov *mom_fpo_vlasov = container_of(momt, struct mom_type_bcorr_fpo_vlasov, momt);
  enum gkyl_vel_edge edge = *(enum gkyl_vel_edge *)param;

  int fidx[GKYL_MAX_DIM];
  double w[GKYL_MAX_DIM];

  int cdim = mom_fpo_vlasov->momt.cdim;
  int pdim = mom_fpo_vlasov->momt.pdim;

  gkyl_copy_int_arr(GKYL_MAX_DIM, idx, fidx);

  int vdim = 3;
  int vel_idx = edge%vdim;
  for (int i=0; i<cdim+vel_idx; ++i) w[i] = xc[i];
  for (int i=cdim+vel_idx; i<pdim-1; ++i) w[i] = xc[i+1];

  long linc = gkyl_range_idx(&mom_fpo_vlasov->phase_range, fidx);

  return mom_bcorr_fpo_vlasov_1x3v_ser_p2(w, fidx, edge, mom_fpo_vlasov->vBoundary, dx,
    (const double*) gkyl_array_cfetch(mom_fpo_vlasov->auxfields.D, linc), 
    f, out);
}

GKYL_CU_DH
static void
kernel_mom_bcorr_fpo_vlasov_2x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_bcorr_fpo_vlasov *mom_fpo_vlasov = container_of(momt, struct mom_type_bcorr_fpo_vlasov, momt);
  enum gkyl_vel_edge edge = *(enum gkyl_vel_edge *)param;

  int fidx[GKYL_MAX_DIM];
  double w[GKYL_MAX_DIM];

  int cdim = mom_fpo_vlasov->momt.cdim;
  int pdim = mom_fpo_vlasov->momt.pdim;

  gkyl_copy_int_arr(GKYL_MAX_DIM, idx, fidx);

  int vdim = 3;
  int vel_idx = edge%vdim;
  for (int i=0; i<cdim+vel_idx; ++i) w[i] = xc[i];
  for (int i=cdim+vel_idx; i<pdim-1; ++i) w[i] = xc[i+1];

  long linc = gkyl_range_idx(&mom_fpo_vlasov->phase_range, fidx);

  return mom_bcorr_fpo_vlasov_2x3v_ser_p1(w, fidx, edge, mom_fpo_vlasov->vBoundary, dx,
    (const double*) gkyl_array_cfetch(mom_fpo_vlasov->auxfields.D, linc), 
    f, out);
}

//
// Serendipity basis kernels
//

// FPO boundary corrections
GKYL_CU_D
static const gkyl_mom_bcorr_fpo_vlasov_kern_list ser_mom_bcorr_fpo_vlasov_kernels[] = {
  // 1x kernels
  { NULL, kernel_mom_bcorr_fpo_vlasov_1x3v_ser_p1, kernel_mom_bcorr_fpo_vlasov_1x3v_ser_p2 }, // 0
  // 2x kernels
  { NULL, kernel_mom_bcorr_fpo_vlasov_2x3v_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_mom_bcorr_fpo_vlasov_free(const struct gkyl_ref_count *ref);
