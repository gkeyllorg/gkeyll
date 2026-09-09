#pragma once

// Private header for twistshift_dg updater, not for direct use in user code.

#include <gkyl_twistshift_dg.h>
#include <gkyl_twistshift_dg_gyrokinetic_kernels.h>
#include <assert.h>
#include <gkyl_mat.h>
#include <gkyl_math.h>
#include <gkyl_eval_on_nodes.h>
#include <string.h> // memcpy

// Function pointer type for twistshift kernels.
typedef void (*twistshift_xlimdg_t)(double sFac, const double *xLimLo,
  const double *xLimUp, double yLimLo, double yLimUp,
  double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);

typedef void (*twistshift_ylimdg_t)(double sFac, double xLimLo,
  double xLimUp, const double *yLimLo, const double *yLimUp,
  double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);

typedef void (*twistshift_fullcell_t)(double dyDo, double yOff,
  const double *ySh, struct gkyl_mat *tsmat);

typedef struct { twistshift_xlimdg_t kernels[3]; }   twistshift_xlimdg_kern_list;  // For use in kernel tables.
typedef struct { twistshift_ylimdg_t kernels[3]; }   twistshift_ylimdg_kern_list;  // For use in kernel tables.
typedef struct { twistshift_fullcell_t kernels[3]; } twistshift_fullcell_kern_list;  // For use in kernel tables.

// Serendipity  kernels.
// p=1 representation of the shift:
static const twistshift_xlimdg_kern_list ser_twistshift_xlimdg_list_0v_yShp1[] = {
  {NULL, twistshift_xlimdg_2x_ser_p1_yshift_p1, NULL,},
  {NULL, twistshift_xlimdg_3x_ser_p1_yshift_p1, NULL,},
};
static const twistshift_ylimdg_kern_list ser_twistshift_ylimdg_list_0v_yShp1[] = {
  {NULL, twistshift_ylimdg_2x_ser_p1_yshift_p1, NULL,},
  {NULL, twistshift_ylimdg_3x_ser_p1_yshift_p1, NULL,},
};
static const twistshift_fullcell_kern_list ser_twistshift_fullcell_list_0v_yShp1[] = {
  {NULL, twistshift_fullcell_2x_ser_p1_yshift_p1, NULL,},
  {NULL, twistshift_fullcell_3x_ser_p1_yshift_p1, NULL,},
};

static const twistshift_xlimdg_kern_list ser_twistshift_xlimdg_list_2v_yShp1[] = {
  {NULL, NULL, NULL,},
  {NULL, twistshift_xlimdg_3x2v_ser_p1_yshift_p1, NULL,},
};
static const twistshift_ylimdg_kern_list ser_twistshift_ylimdg_list_2v_yShp1[] = {
  {NULL, NULL, NULL,},
  {NULL, twistshift_ylimdg_3x2v_ser_p1_yshift_p1, NULL,},
};
static const twistshift_fullcell_kern_list ser_twistshift_fullcell_list_2v_yShp1[] = {
  {NULL, NULL, NULL,},
  {NULL, twistshift_fullcell_3x2v_ser_p1_yshift_p1, NULL,},
};
// p=2 representation of the shift:
static const twistshift_xlimdg_kern_list ser_twistshift_xlimdg_list_0v_yShp2[] = {
  {NULL, twistshift_xlimdg_2x_ser_p1_yshift_p2, NULL,},
  {NULL, twistshift_xlimdg_3x_ser_p1_yshift_p2, NULL,},
};
static const twistshift_ylimdg_kern_list ser_twistshift_ylimdg_list_0v_yShp2[] = {
  {NULL, twistshift_ylimdg_2x_ser_p1_yshift_p2, NULL,},
  {NULL, twistshift_ylimdg_3x_ser_p1_yshift_p2, NULL,},
};
static const twistshift_fullcell_kern_list ser_twistshift_fullcell_list_0v_yShp2[] = {
  {NULL, twistshift_fullcell_2x_ser_p1_yshift_p2, NULL,},
  {NULL, twistshift_fullcell_3x_ser_p1_yshift_p2, NULL,},
};

static const twistshift_xlimdg_kern_list ser_twistshift_xlimdg_list_2v_yShp2[] = {
  {NULL, NULL, NULL,},
  {NULL, NULL, NULL,},
};
static const twistshift_ylimdg_kern_list ser_twistshift_ylimdg_list_2v_yShp2[] = {
  {NULL, NULL, NULL,},
  {NULL, NULL, NULL,},
};
static const twistshift_fullcell_kern_list ser_twistshift_fullcell_list_2v_yShp2[] = {
  {NULL, NULL, NULL,},
  {NULL, NULL, NULL,},
};


struct gkyl_twistshift_dg_kernels {
  twistshift_xlimdg_t xlimdg;
  twistshift_ylimdg_t ylimdg;
  twistshift_fullcell_t fullcell;
};

struct ts_shift_dg_eval_ctx {
  struct gkyl_array *shift_dg; // DG representation of the shift function.
  struct gkyl_basis *shift_b; // Basis for the shift.
  struct gkyl_rect_grid *shear_grid; // shear grid along x.
  struct gkyl_range *shear_r; // Shear grid range.
};

// Primary struct in this updater.
struct gkyl_twistshift_dg {
  int bc_dir; // Direction of the BC is applied in.
  int shift_dir; // Direction of the shift.
  int shear_dir; // Direction the shift varies in (shear).
  enum gkyl_edge_loc edge; // Indicates if BC is for lowe/upper edge.
  struct gkyl_basis basis; // Basis the shifted field is defined with.
  struct gkyl_range local_bcdir_ext_r; // Local range.
  struct gkyl_rect_grid grid; // Grid the shifted field is defined in.
  evalf_t shift_func; // Function defining the shift.
  void *shift_func_ctx; // Context for shift_func.
  struct ts_shift_dg_eval_ctx shift_dg_eval_ctx; // Context for DG shift_func.
  bool use_gpu; // Whether to apply the BC on the GPU.

  struct gkyl_rect_grid shift_grid; // 1D grid in the direction of the shift.
  struct gkyl_range shift_r; // 1D range in the direction of the shift.

  struct gkyl_rect_grid shear_grid; // 1D grid in the direction of the shear.
  struct gkyl_range shear_r; // 1D range in the direction of the shear.

  struct gkyl_rect_grid ts_grid; // Grid the shift twistshift takes place in.
  struct gkyl_range ts_r; // Range the twistshift takes place in.
  int shift_dir_in_ts_grid; // Dimension the shift is in, in the TS grid.
  int shear_dir_in_ts_grid; // Dimension the shear is in, in the TS grid.

  int shift_poly_order; // Poly order of the DG representation of the shift.
  struct gkyl_basis shift_b; // 1D Basis for the DG shift.
  struct gkyl_array *shift_dg; // DG shift.

  int *num_do; // Number of donors at each cell in shear_dir;
  int *shift_dir_idx_do; // Indices of donor cells, in the direction of the
                         // shift, for each cell in the TS grid.

  struct gkyl_twistshift_dg_kernels *kernels;  // kernels for sub-cell integrals.

  // Projection object used in constructing the matrices.
  struct gkyl_eval_on_nodes *ev_on_nod1d;
  // Evaluations of a function at 1D nodes.
  struct gkyl_array *func_nod1d;

  struct gkyl_nmat *scimat; // Subcell integral matrices.
  struct gkyl_nmat *fmat; // Distribution function matrices.
  struct gkyl_nmat *mm_contr; // Contribution resulting from a mat-mat mult.

  long *num_numcol_fidx_do; // 1D indexer, from a index identitying the num-numcol
                            // plane (in the num-numcol-num_basis space), to a
                            // linear index into the donor distribution function f.

  long *num_numcol_fidx_tar; // 1D indexer, from a index identitying the num-numcol
                             // plane (in the num-numcol-num_basis space), to a
                             // linear index into the target distribution function f.

  int *num_do_cum; // Cumulative number of donors up to a give cell in shear_dir;
  struct gkyl_range permutted_ghost_r; // Ghost range to populate in the target
                                       // field, with some dimensions permutted.
  struct gkyl_range ghost_r; // Ghost range this BC fills.
};

#ifdef GKYL_HAVE_CUDA
/**
 * Apply the twist-shift on the NVIDIA GPU. It assumes that periodicity along bc_dir has been
 * applied to the donor field. Can be used in-place.
 *
 * @param up Twist-shift BC updater object.
 * @param fdo Donor field.
 * @param ftar Target field.
 */
void gkyl_twistshift_dg_advance_cu(struct gkyl_twistshift_dg *up, struct gkyl_array *fdo, struct gkyl_array *ftar);
#endif
