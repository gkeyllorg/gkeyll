#pragma once

// private header for array_average updater, not for direct use by user.

#include <gkyl_util.h>
#include <math.h>
#include <gkyl_array_average.h>
#include <gkyl_array_average_kernels.h>
#include <gkyl_dg_bin_ops.h>
#include <assert.h>

// function pointer type for array_average kernels.
typedef void (*array_average_t)( const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);

// for use in kernel tables.
typedef struct { array_average_t kernels[2]; } array_average_kern_list;
typedef struct { array_average_kern_list list[7]; } dim_array_average_kern_list;

GKYL_CU_D
static const dim_array_average_kern_list gkyl_array_average_ker_list[] = {
  { // kernel list for 1x integration
    .list = 
      {
        {gkyl_array_average_1x_ser_p1_avgx, gkyl_array_average_1x_ser_p2_avgx},
        {NULL, NULL},
        {NULL, NULL},
        {NULL, NULL},
        {NULL, NULL},
        {NULL, NULL},
        {NULL, NULL},
      }
  },
  { // kernel list for 2x integration
    .list = 
      {
        {gkyl_array_average_2x_ser_p1_avgx, gkyl_array_average_2x_ser_p2_avgy},
        {gkyl_array_average_2x_ser_p1_avgy, gkyl_array_average_2x_ser_p2_avgy},
        {gkyl_array_average_2x_ser_p1_avgxy, gkyl_array_average_2x_ser_p2_avgxy},
        {NULL, NULL},
        {NULL, NULL},
        {NULL, NULL},
        {NULL, NULL},
      }
  },
  { // kernel list for 3x integration
    .list = 
      {
        {gkyl_array_average_3x_ser_p1_avgx, gkyl_array_average_3x_ser_p2_avgx},
        {gkyl_array_average_3x_ser_p1_avgy, gkyl_array_average_3x_ser_p2_avgy},
        {gkyl_array_average_3x_ser_p1_avgxy, gkyl_array_average_3x_ser_p2_avgxy},
        {gkyl_array_average_3x_ser_p1_avgz, gkyl_array_average_3x_ser_p2_avgz},
        {gkyl_array_average_3x_ser_p1_avgxz, gkyl_array_average_3x_ser_p2_avgxz},
        {gkyl_array_average_3x_ser_p1_avgyz, gkyl_array_average_3x_ser_p2_avgyz},
        {gkyl_array_average_3x_ser_p1_avgxyz, gkyl_array_average_3x_ser_p2_avgxyz},
      }
  }
};

// primary struct in this updater.
struct gkyl_array_average {
  // dimensionality of the full array
  int ndim;
  struct gkyl_basis basis;
  struct gkyl_basis basis_avg;
  // the updater stores the ranges of the input and output
  struct gkyl_range local;
  struct gkyl_range local_avg;

  // if we use gpu or not
  bool use_gpu;

  // array that indicates if the dimension is also a reduced dim
  // (i.e. if the dimension remains)
  int dim_remains[GKYL_MAX_CDIM];
  // array that indicates if the dimension is averaged
  int avg_dim[GKYL_MAX_CDIM];
  // number of averaged dimensions
  int num_avg_dim;
  // number of remaining dimensions
  int num_dim_remain;

  // array that maps the sub dimensions to the full one
  // examples:
  // - for 3x op_yz, then sub_dir = {1,2},
  // - for 2x op_x, then sub_dir = {0}
  int sub_dir[GKYL_MAX_CDIM]; 

  // Single-cell sub-volume element (length for 1D avg, area for 2D)
  double subvol;

  // Single cell average kernel.
  array_average_t kernel;  

  // Pointer to itself on device.
  struct gkyl_array_average *on_dev;
  uint32_t flags;

  // weighted average
  bool isweighted;
  struct gkyl_array *weight;
  struct gkyl_array *weight_avg;
  // inverse volume of the average domain
  double vol_avg_inv; 

  // memory for the weak division at the end of averaging
  gkyl_dg_bin_op_mem *div_mem;

};

GKYL_CU_D static
void gkyl_array_average_choose_kernel(struct gkyl_array_average *up)
{
  int ndim =  up->basis.ndim, poly_order = up->basis.poly_order;

  // We encode the average operations as a binary number 
  // (e.g. 011 = 3 = avgxy, 101 = 5 = avgxz, 111 = 7 = avgxyz)
  int op = -1; // -1 shifted to start with 0
  for (int d = 0; d < ndim; d++)
    op += pow(2,d) * up->avg_dim[d];

  up->kernel = gkyl_array_average_ker_list[ndim-1].list[op].kernels[poly_order-1];

}

#ifdef GKYL_HAVE_CUDA
// Device new functions
struct gkyl_array_average*
gkyl_array_average_cu_dev_new(struct gkyl_array_average *up);
// Device advance functions
void gkyl_array_average_advance_cu(const struct gkyl_array_average *up, 
const struct gkyl_array *fin, struct gkyl_array *avgout);
#endif
