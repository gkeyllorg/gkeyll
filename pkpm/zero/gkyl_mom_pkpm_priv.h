#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_mom_type.h>
#include <gkyl_ref_count.h>
#include <gkyl_mom_pkpm_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct mom_type_pkpm {
  struct gkyl_mom_type momt;
  double mass; // mass of species
};

// for use in kernel tables
typedef struct {
  momf_t kernels[3];
} gkyl_mom_pkpm_kern_list;

GKYL_CU_DH
static void
kernel_mom_pkpm_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_pkpm *mom_pkpm = container_of(momt, struct mom_type_pkpm, momt);

  return mom_pkpm_1x1v_ser_p1(xc, dx, idx, mom_pkpm->mass, f, out);  
}

GKYL_CU_DH
static void
kernel_mom_pkpm_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_pkpm *mom_pkpm = container_of(momt, struct mom_type_pkpm, momt);

  return mom_pkpm_1x1v_ser_p2(xc, dx, idx, mom_pkpm->mass, f, out);  
}

GKYL_CU_DH
static void
kernel_mom_pkpm_1x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_pkpm *mom_pkpm = container_of(momt, struct mom_type_pkpm, momt);

  return mom_pkpm_1x1v_tensor_p2(xc, dx, idx, mom_pkpm->mass, f, out);  
}

GKYL_CU_DH
static void
kernel_mom_pkpm_2x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_pkpm *mom_pkpm = container_of(momt, struct mom_type_pkpm, momt);

  return mom_pkpm_2x1v_ser_p1(xc, dx, idx, mom_pkpm->mass, f, out);  
}

GKYL_CU_DH
static void
kernel_mom_pkpm_2x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_pkpm *mom_pkpm = container_of(momt, struct mom_type_pkpm, momt);

  return mom_pkpm_2x1v_tensor_p2(xc, dx, idx, mom_pkpm->mass, f, out);  
}

GKYL_CU_DH
static void
kernel_mom_pkpm_3x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_pkpm *mom_pkpm = container_of(momt, struct mom_type_pkpm, momt);

  return mom_pkpm_3x1v_ser_p1(xc, dx, idx, mom_pkpm->mass, f, out);  
}

// PKPM coupling moment kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_pkpm_kern_list ser_mom_pkpm_kernels[] = {
  // 1x kernels
  { NULL, kernel_mom_pkpm_1x1v_ser_p1, kernel_mom_pkpm_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, kernel_mom_pkpm_2x1v_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, kernel_mom_pkpm_3x1v_ser_p1, NULL }, // 2
};

// PKPM coupling moment kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_pkpm_kern_list ten_mom_pkpm_kernels[] = {
  // 1x kernels
  { NULL, kernel_mom_pkpm_1x1v_ser_p1, kernel_mom_pkpm_1x1v_tensor_p2 }, // 0
  // 2x kernels
  { NULL, kernel_mom_pkpm_2x1v_ser_p1, kernel_mom_pkpm_2x1v_tensor_p2 }, // 1
  // 3x kernels
  { NULL, kernel_mom_pkpm_3x1v_ser_p1, NULL }, // 2
};

GKYL_CU_DH
static void
kernel_mom_pkpm_diag_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_pkpm *mom_pkpm = container_of(momt, struct mom_type_pkpm, momt);

  return mom_pkpm_diag_1x1v_ser_p1(xc, dx, idx, mom_pkpm->mass, f, out);  
}

GKYL_CU_DH
static void
kernel_mom_pkpm_diag_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_pkpm *mom_pkpm = container_of(momt, struct mom_type_pkpm, momt);

  return mom_pkpm_diag_1x1v_ser_p2(xc, dx, idx, mom_pkpm->mass, f, out);  
}

GKYL_CU_DH
static void
kernel_mom_pkpm_diag_1x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_pkpm *mom_pkpm = container_of(momt, struct mom_type_pkpm, momt);

  return mom_pkpm_diag_1x1v_tensor_p2(xc, dx, idx, mom_pkpm->mass, f, out);  
}

GKYL_CU_DH
static void
kernel_mom_pkpm_diag_2x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_pkpm *mom_pkpm = container_of(momt, struct mom_type_pkpm, momt);

  return mom_pkpm_diag_2x1v_ser_p1(xc, dx, idx, mom_pkpm->mass, f, out);  
}

GKYL_CU_DH
static void
kernel_mom_pkpm_diag_2x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_pkpm *mom_pkpm = container_of(momt, struct mom_type_pkpm, momt);

  return mom_pkpm_diag_2x1v_tensor_p2(xc, dx, idx, mom_pkpm->mass, f, out);  
}

GKYL_CU_DH
static void
kernel_mom_pkpm_diag_3x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_pkpm *mom_pkpm = container_of(momt, struct mom_type_pkpm, momt);

  return mom_pkpm_diag_3x1v_ser_p1(xc, dx, idx, mom_pkpm->mass, f, out);  
}

// PKPM diagnostic kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_pkpm_kern_list ser_mom_pkpm_diag_kernels[] = {
  // 1x kernels
  { NULL, kernel_mom_pkpm_diag_1x1v_ser_p1, kernel_mom_pkpm_diag_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, kernel_mom_pkpm_diag_2x1v_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, kernel_mom_pkpm_diag_3x1v_ser_p1, NULL }, // 2
};

// PKPM diagnostic kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_pkpm_kern_list ten_mom_pkpm_diag_kernels[] = {
  // 1x kernels
  { NULL, kernel_mom_pkpm_diag_1x1v_ser_p1, kernel_mom_pkpm_diag_1x1v_tensor_p2 }, // 0
  // 2x kernels
  { NULL, kernel_mom_pkpm_diag_2x1v_ser_p1, kernel_mom_pkpm_diag_2x1v_tensor_p2 }, // 1
  // 3x kernels
  { NULL, kernel_mom_pkpm_diag_3x1v_ser_p1, NULL }, // 2
};

/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_mom_pkpm_free(const struct gkyl_ref_count *ref);
