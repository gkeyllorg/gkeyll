#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_mom_type.h>
#include <gkyl_ref_count.h>
#include <gkyl_mom_fpo_vlasov_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct mom_type_fpo_vlasov {
  struct gkyl_mom_type momt;
  struct gkyl_range phase_range; // velocity space range
  struct gkyl_mom_fpo_vlasov_auxfields auxfields; // Auxiliary fields.
};

// for use in kernel tables
typedef struct {
  momf_t kernels[3];
} gkyl_mom_fpo_vlasov_kern_list;

GKYL_CU_DH
static void
kernel_mom_fpo_vlasov_1x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_fpo_vlasov *mom_fpo_vlasov = container_of(momt, struct mom_type_fpo_vlasov, momt);

  long linp = gkyl_range_idx(&mom_fpo_vlasov->phase_range, idx);

  return mom_fpo_vlasov_1x3v_ser_p1(xc, dx, idx,
    (const double*) gkyl_array_cfetch(mom_fpo_vlasov->auxfields.a, linp),
    (const double*) gkyl_array_cfetch(mom_fpo_vlasov->auxfields.D, linp),
    f, out);
}

GKYL_CU_DH
static void
kernel_mom_fpo_vlasov_1x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_fpo_vlasov *mom_fpo_vlasov = container_of(momt, struct mom_type_fpo_vlasov, momt);

  long linp = gkyl_range_idx(&mom_fpo_vlasov->phase_range, idx);

  return mom_fpo_vlasov_1x3v_ser_p2(xc, dx, idx,
    (const double*) gkyl_array_cfetch(mom_fpo_vlasov->auxfields.a, linp),
    (const double*) gkyl_array_cfetch(mom_fpo_vlasov->auxfields.D, linp),
    f, out);
}

GKYL_CU_DH
static void
kernel_mom_fpo_vlasov_2x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_fpo_vlasov *mom_fpo_vlasov = container_of(momt, struct mom_type_fpo_vlasov, momt);

  long linp = gkyl_range_idx(&mom_fpo_vlasov->phase_range, idx);

  return mom_fpo_vlasov_2x3v_ser_p1(xc, dx, idx,
    (const double*) gkyl_array_cfetch(mom_fpo_vlasov->auxfields.a, linp),
    (const double*) gkyl_array_cfetch(mom_fpo_vlasov->auxfields.D, linp),
    f, out);
}

//
// Serendipity basis kernels
//

// FPO moments (a + div(D), a . v + div(D . v)) 4 components
GKYL_CU_D
static const gkyl_mom_fpo_vlasov_kern_list ser_mom_fpo_vlasov_kernels[] = {
  // 1x kernels
  { NULL, kernel_mom_fpo_vlasov_1x3v_ser_p1, kernel_mom_fpo_vlasov_1x3v_ser_p2 }, // 0
  // 2x kernels
  { NULL, kernel_mom_fpo_vlasov_2x3v_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_mom_fpo_vlasov_free(const struct gkyl_ref_count *ref);

