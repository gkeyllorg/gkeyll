#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_ref_count.h>
#include <gkyl_mom_bcorr_lbo_pkpm_kernels.h>

struct mom_type_bcorr_lbo_pkpm {
  struct gkyl_mom_type momt;
  double vBoundary[2];
  double mass;
};

// for use in kernel tables
typedef struct {
  momf_t kernels[3];
} gkyl_mom_bcorr_lbo_pkpm_kern_list;

GKYL_CU_DH
static void
kernel_mom_bcorr_lbo_pkpm_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_bcorr_lbo_pkpm *mom_pkpm = container_of(momt, struct mom_type_bcorr_lbo_pkpm, momt);
  enum gkyl_vel_edge edge = *(enum gkyl_vel_edge *)param;

  return mom_bcorr_lbo_pkpm_1x1v_ser_p1(idx, edge, mom_pkpm->vBoundary, dx, mom_pkpm->mass, f, out);  
}

GKYL_CU_DH
static void
kernel_mom_bcorr_lbo_pkpm_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_bcorr_lbo_pkpm *mom_pkpm = container_of(momt, struct mom_type_bcorr_lbo_pkpm, momt);
  enum gkyl_vel_edge edge = *(enum gkyl_vel_edge *)param;

  return mom_bcorr_lbo_pkpm_1x1v_ser_p2(idx, edge, mom_pkpm->vBoundary, dx, mom_pkpm->mass, f, out);  
}

GKYL_CU_DH
static void
kernel_mom_bcorr_lbo_pkpm_1x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_bcorr_lbo_pkpm *mom_pkpm = container_of(momt, struct mom_type_bcorr_lbo_pkpm, momt);
  enum gkyl_vel_edge edge = *(enum gkyl_vel_edge *)param;

  return mom_bcorr_lbo_pkpm_1x1v_tensor_p2(idx, edge, mom_pkpm->vBoundary, dx, mom_pkpm->mass, f, out);  
}

GKYL_CU_DH
static void
kernel_mom_bcorr_lbo_pkpm_2x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_bcorr_lbo_pkpm *mom_pkpm = container_of(momt, struct mom_type_bcorr_lbo_pkpm, momt);
  enum gkyl_vel_edge edge = *(enum gkyl_vel_edge *)param;

  return mom_bcorr_lbo_pkpm_2x1v_ser_p1(idx, edge, mom_pkpm->vBoundary, dx, mom_pkpm->mass, f, out);  
}

GKYL_CU_DH
static void
kernel_mom_bcorr_lbo_pkpm_2x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_bcorr_lbo_pkpm *mom_pkpm = container_of(momt, struct mom_type_bcorr_lbo_pkpm, momt);
  enum gkyl_vel_edge edge = *(enum gkyl_vel_edge *)param;

  return mom_bcorr_lbo_pkpm_2x1v_tensor_p2(idx, edge, mom_pkpm->vBoundary, dx, mom_pkpm->mass, f, out);  
}

GKYL_CU_DH
static void
kernel_mom_bcorr_lbo_pkpm_3x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_bcorr_lbo_pkpm *mom_pkpm = container_of(momt, struct mom_type_bcorr_lbo_pkpm, momt);
  enum gkyl_vel_edge edge = *(enum gkyl_vel_edge *)param;

  return mom_bcorr_lbo_pkpm_3x1v_ser_p1(idx, edge, mom_pkpm->vBoundary, dx, mom_pkpm->mass, f, out);  
}

// Moment boundary correction kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_bcorr_lbo_pkpm_kern_list ser_mom_bcorr_lbo_pkpm_kernels[] = {
  // 1x kernels
  { NULL, kernel_mom_bcorr_lbo_pkpm_1x1v_ser_p1, kernel_mom_bcorr_lbo_pkpm_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, kernel_mom_bcorr_lbo_pkpm_2x1v_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, kernel_mom_bcorr_lbo_pkpm_3x1v_ser_p1, NULL }, // 2
};

// Moment boundary correction kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_bcorr_lbo_pkpm_kern_list ten_mom_bcorr_lbo_pkpm_kernels[] = {
  // 1x kernels
  { NULL, kernel_mom_bcorr_lbo_pkpm_1x1v_ser_p1, kernel_mom_bcorr_lbo_pkpm_1x1v_tensor_p2 }, // 0
  // 2x kernels
  { NULL, kernel_mom_bcorr_lbo_pkpm_2x1v_ser_p1, kernel_mom_bcorr_lbo_pkpm_2x1v_tensor_p2 }, // 1
  // 3x kernels
  { NULL, kernel_mom_bcorr_lbo_pkpm_3x1v_ser_p1, NULL }, // 2
};

/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_mom_bcorr_lbo_pkpm_free(const struct gkyl_ref_count *ref);
