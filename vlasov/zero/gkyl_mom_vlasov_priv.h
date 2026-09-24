#pragma once

// Private header, not for direct use in user code

#include <gkyl_mom_type.h>
#include <gkyl_ref_count.h>
#include <gkyl_mom_vlasov_kernels.h>

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

typedef void (*vlasov_momf_t)(const double *xc, const double *dx,
  const int *idx, const double *fIn, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vlasov_momf_t kernels[3]; } gkyl_mom_kern_list;

// M0 kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_kern_list ser_m0_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M0_1x1v_ser_p1, vlasov_M0_1x1v_ser_p2 }, // 0
  { NULL, vlasov_M0_1x2v_ser_p1, vlasov_M0_1x2v_ser_p2 }, // 1
  { NULL, vlasov_M0_1x3v_ser_p1, vlasov_M0_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M0_2x2v_ser_p1, vlasov_M0_2x2v_ser_p2 }, // 3
  { NULL, vlasov_M0_2x3v_ser_p1, vlasov_M0_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_M0_3x3v_ser_p1, NULL                  }, // 5
};

// M0 kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_kern_list tensor_m0_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M0_1x1v_tensor_p1, vlasov_M0_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_M0_1x2v_tensor_p1, vlasov_M0_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_M0_1x3v_tensor_p1, vlasov_M0_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M0_2x2v_tensor_p1, vlasov_M0_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_M0_2x3v_tensor_p1, vlasov_M0_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_M0_3x3v_tensor_p1, NULL                  }, // 5
};

// M1i kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_kern_list ser_m1i_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M1i_1x1v_ser_p1, vlasov_M1i_1x1v_ser_p2 }, // 0
  { NULL, vlasov_M1i_1x2v_ser_p1, vlasov_M1i_1x2v_ser_p2 }, // 1
  { NULL, vlasov_M1i_1x3v_ser_p1, vlasov_M1i_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M1i_2x2v_ser_p1, vlasov_M1i_2x2v_ser_p2 }, // 3
  { NULL, vlasov_M1i_2x3v_ser_p1, vlasov_M1i_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_M1i_3x3v_ser_p1, NULL                   }, // 5
};

// M1i kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_kern_list tensor_m1i_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M1i_1x1v_tensor_p1, vlasov_M1i_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_M1i_1x2v_tensor_p1, vlasov_M1i_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_M1i_1x3v_tensor_p1, vlasov_M1i_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M1i_2x2v_tensor_p1, vlasov_M1i_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_M1i_2x3v_tensor_p1, vlasov_M1i_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_M1i_3x3v_tensor_p1, NULL                   }, // 5
};


// M2 kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_kern_list ser_m2_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M2_1x1v_ser_p1, vlasov_M2_1x1v_ser_p2 }, // 0
  { NULL, vlasov_M2_1x2v_ser_p1, vlasov_M2_1x2v_ser_p2 }, // 1
  { NULL, vlasov_M2_1x3v_ser_p1, vlasov_M2_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M2_2x2v_ser_p1, vlasov_M2_2x2v_ser_p2 }, // 3
  { NULL, vlasov_M2_2x3v_ser_p1, vlasov_M2_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_M2_3x3v_ser_p1, NULL                  }, // 5
};

// M2 kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_kern_list tensor_m2_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M2_1x1v_tensor_p1, vlasov_M2_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_M2_1x2v_tensor_p1, vlasov_M2_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_M2_1x3v_tensor_p1, vlasov_M2_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M2_2x2v_tensor_p1, vlasov_M2_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_M2_2x3v_tensor_p1, vlasov_M2_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_M2_3x3v_tensor_p1, NULL                  }, // 5
};

// M2ij kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_kern_list ser_m2ij_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M2ij_1x1v_ser_p1, vlasov_M2ij_1x1v_ser_p2 }, // 0
  { NULL, vlasov_M2ij_1x2v_ser_p1, vlasov_M2ij_1x2v_ser_p2 }, // 1
  { NULL, vlasov_M2ij_1x3v_ser_p1, vlasov_M2ij_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M2ij_2x2v_ser_p1, vlasov_M2ij_2x2v_ser_p2 }, // 3
  { NULL, vlasov_M2ij_2x3v_ser_p1, vlasov_M2ij_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_M2ij_3x3v_ser_p1, NULL                    }, // 5
};

// M2ij kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_kern_list tensor_m2ij_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M2ij_1x1v_tensor_p1, vlasov_M2ij_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_M2ij_1x2v_tensor_p1, vlasov_M2ij_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_M2ij_1x3v_tensor_p1, vlasov_M2ij_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M2ij_2x2v_tensor_p1, vlasov_M2ij_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_M2ij_2x3v_tensor_p1, vlasov_M2ij_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_M2ij_3x3v_tensor_p1, NULL                    }, // 5
};

// M3i kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_kern_list ser_m3i_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M3i_1x1v_ser_p1, vlasov_M3i_1x1v_ser_p2 }, // 0
  { NULL, vlasov_M3i_1x2v_ser_p1, vlasov_M3i_1x2v_ser_p2 }, // 1
  { NULL, vlasov_M3i_1x3v_ser_p1, vlasov_M3i_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M3i_2x2v_ser_p1, vlasov_M3i_2x2v_ser_p2 }, // 3
  { NULL, vlasov_M3i_2x3v_ser_p1, vlasov_M3i_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_M3i_3x3v_ser_p1, NULL                   }, // 5
};

// M3i kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_kern_list tensor_m3i_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M3i_1x1v_tensor_p1, vlasov_M3i_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_M3i_1x2v_tensor_p1, vlasov_M3i_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_M3i_1x3v_tensor_p1, vlasov_M3i_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M3i_2x2v_tensor_p1, vlasov_M3i_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_M3i_2x3v_tensor_p1, vlasov_M3i_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_M3i_3x3v_tensor_p1, NULL                   }, // 5
};

// M3ijk kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_kern_list ser_m3ijk_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M3ijk_1x1v_ser_p1, vlasov_M3ijk_1x1v_ser_p2 }, // 0
  { NULL, vlasov_M3ijk_1x2v_ser_p1, vlasov_M3ijk_1x2v_ser_p2 }, // 1
  { NULL, vlasov_M3ijk_1x3v_ser_p1, vlasov_M3ijk_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M3ijk_2x2v_ser_p1, vlasov_M3ijk_2x2v_ser_p2 }, // 3
  { NULL, vlasov_M3ijk_2x3v_ser_p1, vlasov_M3ijk_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_M3ijk_3x3v_ser_p1, NULL                     }, // 5
};

// M3ijk kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_kern_list tensor_m3ijk_kernels[] = {
  // 1x kernels
  { NULL, vlasov_M3ijk_1x1v_tensor_p1, vlasov_M3ijk_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_M3ijk_1x2v_tensor_p1, vlasov_M3ijk_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_M3ijk_1x3v_tensor_p1, vlasov_M3ijk_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_M3ijk_2x2v_tensor_p1, vlasov_M3ijk_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_M3ijk_2x3v_tensor_p1, vlasov_M3ijk_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_M3ijk_3x3v_tensor_p1, NULL                     }, // 5
};

// Five moments (Zeroth, First, and Second moment together) kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_kern_list ser_five_moments_kernels[] = {
  // 1x kernels
  { NULL, vlasov_five_moments_1x1v_ser_p1, vlasov_five_moments_1x1v_ser_p2 }, // 0
  { NULL, vlasov_five_moments_1x2v_ser_p1, vlasov_five_moments_1x2v_ser_p2 }, // 1
  { NULL, vlasov_five_moments_1x3v_ser_p1, vlasov_five_moments_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_five_moments_2x2v_ser_p1, vlasov_five_moments_2x2v_ser_p2 }, // 3
  { NULL, vlasov_five_moments_2x3v_ser_p1, vlasov_five_moments_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_five_moments_3x3v_ser_p1, NULL                     }, // 5
};

// Five moments (Zeroth, First, and Second moment together) kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_kern_list tensor_five_moments_kernels[] = {
  // 1x kernels
  { NULL, vlasov_five_moments_1x1v_tensor_p1, vlasov_five_moments_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_five_moments_1x2v_tensor_p1, vlasov_five_moments_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_five_moments_1x3v_tensor_p1, vlasov_five_moments_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_five_moments_2x2v_tensor_p1, vlasov_five_moments_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_five_moments_2x3v_tensor_p1, vlasov_five_moments_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_five_moments_3x3v_tensor_p1, NULL                     }, // 5
};

// Integrated moments kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_kern_list ser_int_five_moments_kernels[] = {
  // 1x kernels
  { NULL, vlasov_int_five_moments_1x1v_ser_p1, vlasov_int_five_moments_1x1v_ser_p2 }, // 0
  { NULL, vlasov_int_five_moments_1x2v_ser_p1, vlasov_int_five_moments_1x2v_ser_p2 }, // 1
  { NULL, vlasov_int_five_moments_1x3v_ser_p1, vlasov_int_five_moments_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_int_five_moments_2x2v_ser_p1, vlasov_int_five_moments_2x2v_ser_p2 }, // 3
  { NULL, vlasov_int_five_moments_2x3v_ser_p1, vlasov_int_five_moments_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_int_five_moments_3x3v_ser_p1, NULL                     }, // 5
};

// Integrated moments kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_kern_list tensor_int_five_moments_kernels[] = {
  // 1x kernels
  { NULL, vlasov_int_five_moments_1x1v_tensor_p1, vlasov_int_five_moments_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_int_five_moments_1x2v_tensor_p1, vlasov_int_five_moments_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_int_five_moments_1x3v_tensor_p1, vlasov_int_five_moments_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_int_five_moments_2x2v_tensor_p1, vlasov_int_five_moments_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_int_five_moments_2x3v_tensor_p1, vlasov_int_five_moments_2x3v_tensor_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_int_five_moments_3x3v_tensor_p1, NULL                     }, // 5
};

struct mom_type_vlasov {
  struct gkyl_mom_type momt;
  vlasov_momf_t kernel; // moment calculation kernel
};

/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_mom_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static void
kernel(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  return mom_vlasov->kernel(xc, dx, idx, f, out);
}

#ifdef GKYL_HAVE_CUDA
/**
 * Create new Vlasov moment type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_mom_type* 
gkyl_mom_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, enum gkyl_distribution_moments mom_type);

/**
 * Create new integrated Vlasov moment type object on NV-GPU:
 * see new() method above for documentation.
 */
struct gkyl_mom_type* 
gkyl_int_mom_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, enum gkyl_distribution_moments mom_type);
#endif
