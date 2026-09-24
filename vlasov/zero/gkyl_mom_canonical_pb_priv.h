#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_mom_type.h>
#include <gkyl_ref_count.h>
#include <gkyl_mom_canonical_pb_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct mom_type_canonical_pb {
  struct gkyl_mom_type momt;
  double mass; // mass of species
  struct gkyl_range phase_range; // phase space range
  struct gkyl_mom_canonical_pb_auxfields auxfields; // Auxiliary fields.
};

// The cv_index[cd].vdim[cd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct {
  momf_t kernels[3];
} gkyl_canonical_pb_mom_kern_list;

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_1x1v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_1x1v_ser_p2(dx, 
  (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_1x2v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_1x2v_ser_p2(dx, 
  (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_1x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_1x3v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_1x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_1x3v_ser_p2(dx, 
  (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_2x2v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_2x2v_ser_p2(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_2x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_2x3v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_2x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_2x3v_ser_p2(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}


GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_1x1v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_1x1v_ser_p2(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out); 
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_1x2v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_1x2v_ser_p2(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out); 
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_1x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_1x3v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_1x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_1x3v_ser_p2(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out); 
}


GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_2x2v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_2x2v_ser_p2(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_2x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_2x3v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_2x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_2x3v_ser_p2(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);
}


//
// Serendipity basis kernels
//

// MEnergy kernel list
GKYL_CU_D
static const gkyl_canonical_pb_mom_kern_list ser_menergy_kernels[] = {
  // 1x kernels
  { NULL, kernel_canonical_pb_MEnergy_1x1v_ser_p1, kernel_canonical_pb_MEnergy_1x1v_ser_p2 }, // 0
  { NULL, kernel_canonical_pb_MEnergy_1x2v_ser_p1, kernel_canonical_pb_MEnergy_1x2v_ser_p2 }, // 1
  { NULL, kernel_canonical_pb_MEnergy_1x3v_ser_p1, kernel_canonical_pb_MEnergy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_canonical_pb_MEnergy_2x2v_ser_p1, kernel_canonical_pb_MEnergy_2x2v_ser_p2 }, // 3
  { NULL, kernel_canonical_pb_MEnergy_2x3v_ser_p1, kernel_canonical_pb_MEnergy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                            }, // 5
};

// Integrated moment kernel list
GKYL_CU_D
static const gkyl_canonical_pb_mom_kern_list ser_int_five_moments_kernels[] = {
  // 1x kernels
  { NULL, kernel_canonical_pb_int_five_moments_1x1v_ser_p1, kernel_canonical_pb_int_five_moments_1x1v_ser_p2 }, // 0
  { NULL, kernel_canonical_pb_int_five_moments_1x2v_ser_p1, kernel_canonical_pb_int_five_moments_1x2v_ser_p2 }, // 1
  { NULL, kernel_canonical_pb_int_five_moments_1x3v_ser_p1, kernel_canonical_pb_int_five_moments_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_canonical_pb_int_five_moments_2x2v_ser_p1, kernel_canonical_pb_int_five_moments_2x2v_ser_p2 }, // 3
  { NULL, kernel_canonical_pb_int_five_moments_2x3v_ser_p1, kernel_canonical_pb_int_five_moments_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                                 }, // 5
};


GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_1x1v_tensor_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_1x1v_tensor_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_1x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_1x1v_tensor_p2(dx, 
  (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_1x2v_tensor_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_1x2v_tensor_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_1x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_1x2v_tensor_p2(dx, 
  (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_1x3v_tensor_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_1x3v_tensor_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_1x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_1x3v_tensor_p2(dx, 
  (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_2x2v_tensor_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_2x2v_tensor_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_2x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_2x2v_tensor_p2(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_2x3v_tensor_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_2x3v_tensor_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_MEnergy_3x3v_tensor_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_MEnergy_3x3v_tensor_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_1x1v_tensor_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_1x1v_tensor_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_1x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_1x1v_tensor_p2(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out); 
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_1x2v_tensor_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_1x2v_tensor_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_1x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_1x2v_tensor_p2(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out); 
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_1x3v_tensor_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_1x3v_tensor_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_1x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_1x3v_tensor_p2(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out); 
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_2x2v_tensor_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_2x2v_tensor_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_2x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_2x2v_tensor_p2(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_2x3v_tensor_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_2x3v_tensor_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);
}

GKYL_CU_DH
static void
kernel_canonical_pb_int_five_moments_3x3v_tensor_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_int_five_moments_3x3v_tensor_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);
}

//
// Tensor basis kernels
//

// MEnergy kernel list
GKYL_CU_D
static const gkyl_canonical_pb_mom_kern_list tensor_menergy_kernels[] = {
  // 1x kernels
  { NULL, kernel_canonical_pb_MEnergy_1x1v_tensor_p1, kernel_canonical_pb_MEnergy_1x1v_tensor_p2 }, // 0
  { NULL, kernel_canonical_pb_MEnergy_1x2v_tensor_p1, kernel_canonical_pb_MEnergy_1x2v_tensor_p2 }, // 1
  { NULL, kernel_canonical_pb_MEnergy_1x3v_tensor_p1, kernel_canonical_pb_MEnergy_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, kernel_canonical_pb_MEnergy_2x2v_tensor_p1, kernel_canonical_pb_MEnergy_2x2v_tensor_p2 }, // 3
  { NULL, kernel_canonical_pb_MEnergy_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, kernel_canonical_pb_MEnergy_3x3v_tensor_p1, NULL                            }, // 5
};

// Integrated moment kernel list
GKYL_CU_D
static const gkyl_canonical_pb_mom_kern_list tensor_int_five_moments_kernels[] = {
  // 1x kernels
  { NULL, kernel_canonical_pb_int_five_moments_1x1v_tensor_p1, kernel_canonical_pb_int_five_moments_1x1v_tensor_p2 }, // 0
  { NULL, kernel_canonical_pb_int_five_moments_1x2v_tensor_p1, kernel_canonical_pb_int_five_moments_1x2v_tensor_p2 }, // 1
  { NULL, kernel_canonical_pb_int_five_moments_1x3v_tensor_p1, kernel_canonical_pb_int_five_moments_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, kernel_canonical_pb_int_five_moments_2x2v_tensor_p1, kernel_canonical_pb_int_five_moments_2x2v_tensor_p2 }, // 3
  { NULL, kernel_canonical_pb_int_five_moments_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, kernel_canonical_pb_int_five_moments_3x3v_tensor_p1, NULL                                 }, // 5
};


GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_1x1v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_1x1v_ser_p2(dx, 
  (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_1x2v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_1x2v_ser_p2(dx, 
  (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_1x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_1x3v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_1x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_1x3v_ser_p2(dx, 
  (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_2x2v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_2x2v_ser_p2(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_2x3v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_2x3v_ser_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_2x3v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_2x3v_ser_p2(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}


//
// Serendipity basis kernels
//

// M1i_from_H kernel list
GKYL_CU_D
static const gkyl_canonical_pb_mom_kern_list ser_m1i_from_h_kernels[] = {
  // 1x kernels
  { NULL, kernel_canonical_pb_M1i_from_H_1x1v_ser_p1, kernel_canonical_pb_M1i_from_H_1x1v_ser_p2 }, // 0
  { NULL, kernel_canonical_pb_M1i_from_H_1x2v_ser_p1, kernel_canonical_pb_M1i_from_H_1x2v_ser_p2 }, // 1
  { NULL, kernel_canonical_pb_M1i_from_H_1x3v_ser_p1, kernel_canonical_pb_M1i_from_H_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, kernel_canonical_pb_M1i_from_H_2x2v_ser_p1, kernel_canonical_pb_M1i_from_H_2x2v_ser_p2 }, // 3
  { NULL, kernel_canonical_pb_M1i_from_H_2x3v_ser_p1, kernel_canonical_pb_M1i_from_H_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, NULL, NULL                            }, // 5
};

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_1x1v_tensor_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_1x1v_tensor_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_1x1v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_1x1v_tensor_p2(dx, 
  (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_1x2v_tensor_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_1x2v_tensor_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_1x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_1x2v_tensor_p2(dx, 
  (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_1x3v_tensor_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_1x3v_tensor_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_1x3v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_1x3v_tensor_p2(dx, 
  (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_2x2v_tensor_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_2x2v_tensor_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_2x2v_tensor_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_2x2v_tensor_p2(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_2x3v_tensor_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_2x3v_tensor_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

GKYL_CU_DH
static void
kernel_canonical_pb_M1i_from_H_3x3v_tensor_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  long pidx = gkyl_range_idx(&mom_can_pb->phase_range, idx);

  return canonical_pb_M1i_from_H_3x3v_tensor_p1(dx, 
    (const double*) gkyl_array_cfetch(mom_can_pb->auxfields.hamil, pidx), f, out);  
}

//
// Tensor basis kernels
//

// M1i_from_H kernel list
GKYL_CU_D
static const gkyl_canonical_pb_mom_kern_list tensor_m1i_from_h_kernels[] = {
  // 1x kernels
  { NULL, kernel_canonical_pb_M1i_from_H_1x1v_tensor_p1, kernel_canonical_pb_M1i_from_H_1x1v_tensor_p2 }, // 0
  { NULL, kernel_canonical_pb_M1i_from_H_1x2v_tensor_p1, kernel_canonical_pb_M1i_from_H_1x2v_tensor_p2 }, // 1
  { NULL, kernel_canonical_pb_M1i_from_H_1x3v_tensor_p1, kernel_canonical_pb_M1i_from_H_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, kernel_canonical_pb_M1i_from_H_2x2v_tensor_p1, kernel_canonical_pb_M1i_from_H_2x2v_tensor_p2 }, // 3
  { NULL, kernel_canonical_pb_M1i_from_H_2x3v_tensor_p1, NULL }, // 4
  // 3x kernels
  { NULL, kernel_canonical_pb_M1i_from_H_3x3v_tensor_p1, NULL                            }, // 5
};


/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_mom_can_pb_free(const struct gkyl_ref_count *ref);

#ifdef GKYL_HAVE_CUDA
/**
 * Create new canonical-pb moment type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_mom_type* gkyl_mom_canonical_pb_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range, enum gkyl_distribution_moments mom_type);

/**
 * Create new canonical-pb integrated moment type
 * object on NV-GPU: see new() method above for documentation.
 */
struct gkyl_mom_type* gkyl_int_mom_canonical_pb_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range, enum gkyl_distribution_moments mom_type);

/**
 * CUDA device function to set auxiliary fields needed in computing moments.
 * 
 * @param momt moment type.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_mom_canonical_pb_set_auxfields_cu(const struct gkyl_mom_type *momt,
  struct gkyl_mom_canonical_pb_auxfields auxin);

#endif
