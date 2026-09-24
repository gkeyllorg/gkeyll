#pragma once

// Private header for translate_dim updater, not for direct use in user code.

#include <gkyl_translate_dim.h>
#include <gkyl_translate_dim_kernels.h>
#include <gkyl_util.h>
#include <assert.h>

static void
translate_dim_range_check_conf_deflate(int dir, int cdim_do, int cdim_tar, int vdim,
  const struct gkyl_range *rng_do, const struct gkyl_range *rng_tar)
{
  int c = 0;
  for (int d=0; d<cdim_do; d++) {
    if (d != dir) {
      assert(rng_tar->lower[c] == rng_do->lower[d]);
      assert(rng_tar->upper[c] == rng_do->upper[d]);
      c++;
    }
  }
}

static void
translate_dim_range_check_conf_inflate(int dir, int cdim_do, int cdim_tar, int vdim,
  const struct gkyl_range *rng_do, const struct gkyl_range *rng_tar)
{
  for (int d=0; d<cdim_do-1; d++) {
    assert(rng_do->lower[d] == rng_tar->lower[d]);
    assert(rng_do->upper[d] == rng_tar->upper[d]);
  }
  assert(rng_do->lower[cdim_do-1] == rng_tar->lower[cdim_tar-1]);
  assert(rng_do->upper[cdim_do-1] == rng_tar->upper[cdim_tar-1]);
}

static void
translate_dim_range_check_phase_deflate(int dir, int cdim_do, int cdim_tar, int vdim,
  const struct gkyl_range *rng_do, const struct gkyl_range *rng_tar)
{
  translate_dim_range_check_conf_deflate(dir, cdim_do, cdim_tar, vdim,
    rng_do, rng_tar);
  for (int d=0; d<vdim; d++) {
    assert(rng_do->lower[cdim_do+d] == rng_tar->lower[cdim_tar+d]);
    assert(rng_do->upper[cdim_do+d] == rng_tar->upper[cdim_tar+d]);
  };
}

static void
translate_dim_range_check_phase_inflate(int dir, int cdim_do, int cdim_tar, int vdim,
  const struct gkyl_range *rng_do, const struct gkyl_range *rng_tar)
{
  translate_dim_range_check_conf_inflate(dir, cdim_do, cdim_tar, vdim,
    rng_do, rng_tar);
  for (int d=0; d<vdim; d++) {
    assert(rng_do->lower[cdim_do+d] == rng_tar->lower[cdim_tar+d]);
    assert(rng_do->upper[cdim_do+d] == rng_tar->upper[cdim_tar+d]);
  };
}

// Function pointer type for sheath reflection kernels.
typedef void (*translate_dim_t)(const double *fdo, double *ftar);

typedef struct {translate_dim_t kernels[3];} trans_dim_kern_list;  // For use in kernel tables.
typedef struct {trans_dim_kern_list list[9];} trans_dim_kern_list_updown;  // For use in kernel tables.

// Serendipity  kernels.
GKYL_CU_D
static const trans_dim_kern_list_updown trans_dim_kern_list_ser[] = {
  // 1x
  { .list = {
      {NULL, NULL, NULL },
      {NULL, NULL, NULL },
      {NULL, NULL, NULL },
      {translate_dim_1x_ser_p1_to_2x_p1, NULL, NULL },
      {NULL, NULL, NULL },
      {NULL, NULL, NULL },
      {NULL, NULL, NULL },
      {NULL, NULL, NULL },
      {NULL, NULL, NULL },
    },
  },
  // 2x
  { .list = {
      { translate_dim_2x_ser_p1_to_1x_p1_dirx_lo , NULL, NULL },
      { translate_dim_2x_ser_p1_to_1x_p1_dirx_mid, NULL, NULL },
      { translate_dim_2x_ser_p1_to_1x_p1_dirx_up , NULL, NULL },
      { translate_dim_2x_ser_p1_to_1x_p1_diry_lo , NULL, NULL },
      { translate_dim_2x_ser_p1_to_1x_p1_diry_mid, NULL, NULL },
      { translate_dim_2x_ser_p1_to_1x_p1_diry_up , NULL, NULL },
      { translate_dim_2x_ser_p1_to_3x_p1, NULL, NULL },
      {NULL, NULL, NULL },
      {NULL, NULL, NULL },
    },
  },
  // 3x
  { .list = {
      { translate_dim_3x_ser_p1_to_2x_p1_dirx_lo , NULL, NULL },
      { translate_dim_3x_ser_p1_to_2x_p1_dirx_mid, NULL, NULL },
      { translate_dim_3x_ser_p1_to_2x_p1_dirx_up , NULL, NULL },
      { translate_dim_3x_ser_p1_to_2x_p1_diry_lo , NULL, NULL },
      { translate_dim_3x_ser_p1_to_2x_p1_diry_mid, NULL, NULL },
      { translate_dim_3x_ser_p1_to_2x_p1_diry_up , NULL, NULL },
      { translate_dim_3x_ser_p1_to_2x_p1_dirz_lo , NULL, NULL },
      { translate_dim_3x_ser_p1_to_2x_p1_dirz_mid, NULL, NULL },
      { translate_dim_3x_ser_p1_to_2x_p1_dirz_up , NULL, NULL },
    },
  },
};

// GkHybrid kernels.
GKYL_CU_D
static const trans_dim_kern_list trans_dim_kern_list_gkhyb[] = {
  { translate_dim_gyrokinetic_2x2v_ser_p1_from_1x2v_p1, NULL, NULL },
  { translate_dim_gyrokinetic_3x2v_ser_p1_from_1x2v_p1, NULL, NULL },
  { translate_dim_gyrokinetic_3x2v_ser_p1_from_2x2v_p1, NULL, NULL },
};

struct gkyl_translate_dim_kernels {
  translate_dim_t translate;  // Kernel that translate the DG coefficients.
  void (*get_idx_do)(int cdim_tar, int vdim, const int *idx_tar,
    const struct gkyl_range *rng_do, int cdim_do, int *idx_do, int dir);
};

GKYL_CU_DH
static void
translate_dim_get_idx_do_gk(int cdim_tar, int vdim, const int *idx_tar,
  const struct gkyl_range *rng_do, int cdim_do, int *idx_do, int dir)
{
  for (int d=0; d<cdim_do-1; d++) idx_do[d] = idx_tar[d]; 
  idx_do[cdim_do-1] = idx_tar[cdim_tar-1]; 
  for (int d=0; d<vdim; d++) idx_do[cdim_do+d] = idx_tar[cdim_tar+d]; 
}

GKYL_CU_DH
static void
translate_dim_get_idx_do_conf_down(int cdim_tar, int vdim, const int *idx_tar,
  const struct gkyl_range *rng_do, int cdim_do, int *idx_do, int dir)
{
  int c = 0;
  for (int d=0; d<cdim_do; d++) {
    if (d != dir) {
      idx_do[d] = idx_tar[c]; 
      c++;
    }
  }
  idx_do[dir] = rng_do->lower[dir];
}

GKYL_CU_DH
static void
translate_dim_get_idx_do_conf_up(int cdim_tar, int vdim, const int *idx_tar,
  const struct gkyl_range *rng_do, int cdim_do, int *idx_do, int dir)
{
  for (int d=0; d<cdim_do; d++) idx_do[d] = idx_tar[d]; 
}

// Primary struct in this updater.
struct gkyl_translate_dim {
  int cdim_do; // Configuration space dimension of donor field.
  int vdim_do; // Velocity space dimension of donor field.
  int cdim_tar; // Configuration space dimension of target field.
  int vdim_tar; // Velocity space dimension of target field.
  int num_basis_do; // Number of polynomials in donor basis.
  int num_basis_tar; // Number of polynomials in target basis.
  int dir; // Direction to remove when going to lower dimensionality.
  bool use_gpu;
  struct gkyl_translate_dim_kernels *kernels;
  void (*range_check_func)(int dir, int cdim_do, int cdim_tar, int vdim,
    const struct gkyl_range *rng_do, const struct gkyl_range *rng_tar);
};

#ifdef GKYL_HAVE_CUDA
// Declaration of cuda device functions.
void trans_dim_choose_kernel_cu(struct gkyl_translate_dim_kernels *kernels,
  int cdim_do, struct gkyl_basis basis_do, int cdim_tar, struct gkyl_basis basis_tar,
  int dir, enum gkyl_edge_loc edge);

void gkyl_translate_dim_advance_cu(gkyl_translate_dim* up,
  const struct gkyl_range *rng_do, const struct gkyl_range *rng_tar,
  const struct gkyl_array *GKYL_RESTRICT fdo, int ncomp,
  struct gkyl_array *GKYL_RESTRICT ftar);
#endif

GKYL_CU_D
static void trans_dim_choose_kernel(struct gkyl_translate_dim_kernels *kernels, int cdim_do,
  struct gkyl_basis basis_do, int cdim_tar, struct gkyl_basis basis_tar, int dir, enum gkyl_edge_loc edge, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    trans_dim_choose_kernel_cu(kernels, cdim_do, basis_do, cdim_tar, basis_tar, dir, edge);
    return;
  }
#endif

  enum gkyl_basis_type basis_type = basis_tar.b_type;
  int poly_order = basis_tar.poly_order;
  int dir_idx = cdim_tar-1;
  int edge_idx = 0;
  if (cdim_tar < cdim_do) {
    dir_idx = dir;
    switch (edge) {
      case GKYL_LOWER_EDGE:
        edge_idx = 0;
        break;
      case GKYL_NO_EDGE:
        edge_idx = 1;
        break;
      case GKYL_UPPER_EDGE:
        edge_idx = 2;
        break;
    }
  }

  // Choose kernel that translates DG coefficients.
  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
      kernels->translate = trans_dim_kern_list_gkhyb[cdim_tar+cdim_do-3].kernels[poly_order-1];
      break;
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->translate = trans_dim_kern_list_ser[cdim_do-1].list[dir_idx*3+edge_idx].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break;
  }

  // Choose the function that populates the donor index.
  int vdim = basis_do.ndim - cdim_do;
  if (vdim > 0) {
    kernels->get_idx_do = translate_dim_get_idx_do_gk;
  }
  else {
    if (cdim_tar < cdim_do)
      kernels->get_idx_do = translate_dim_get_idx_do_conf_down;
    else
      kernels->get_idx_do = translate_dim_get_idx_do_conf_up;
  }
}
