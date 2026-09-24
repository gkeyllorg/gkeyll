#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_dg_basis_ops_priv.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_ref_count.h>
#include <gkyl_util.h>
#include <float.h>

void
gkyl_dg_basis_ops_eval_array_at_coord_comp(const struct gkyl_array *arr, const double *coord,
  const struct gkyl_basis *basis, const struct gkyl_rect_grid *grid, const struct gkyl_range *rng,
  double *out)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) {
    gkyl_dg_basis_ops_eval_array_at_coord_comp_cu(arr, coord, basis, grid, rng, out);
    return;
  }
#endif

  // Get the index of cell that contains coord.
  int coord_idx[grid->ndim];
  gkyl_rect_grid_coord_idx(grid, coord, coord_idx);

  if (gkyl_range_contains_idx(rng, coord_idx)) {
    // Logical coordinate.
    double xc[grid->ndim], coord_log[grid->ndim];
    gkyl_rect_grid_cell_center(grid, coord_idx, xc);
    for (int d=0; d<grid->ndim; d++)
      coord_log[d] = (2.0/grid->dx[d])*(coord[d]-xc[d]);

    long linidx = gkyl_range_idx(rng, coord_idx);
    const double *arr_c = gkyl_array_cfetch(arr, linidx);
    out[0] = basis->eval_expand(coord_log, arr_c);
  }
  else
    out[0] = -DBL_MAX;
}

enum dg_basis_op_code {
  GKYL_DG_BASIS_OP_CUBIC_1D,
  GKYL_DG_BASIS_OP_CUBIC_2D
};

struct gkyl_dg_basis_op_mem {
  enum dg_basis_op_code opcode;
  union {
    // data for GKYL_DG_BASIS_OP_CUBIC_1D
    struct {
      struct gkyl_array *grad1dx;
    };
    // data for GKYL_DG_BASIS_OP_CUBIC_2D
    struct {
      struct gkyl_array *grad2dx, *grad2dy, *grad2dxy;
    };
  };
};

struct dg_basis_ops_evalf_ctx {
  int ndim; // number of dimensions
  int cells[2]; // cells in each direction
  double dx[2]; // cell-spacing in each direction
  struct gkyl_rect_grid grid; // grid on which cubic is define
  struct gkyl_range local, local_ext; // ranges for cubic
  struct gkyl_basis basis; // p=3 basis functions
  
  struct gkyl_array *cubic; // cubic DG representation
};

void
gkyl_dg_calc_cubic_1d(const double val[2], const double grad[2], double *coeff)
{
  coeff[0] = 0.7071067811865475*val[1]-0.2357022603955158*grad[1]+0.7071067811865475*val[0]+0.2357022603955158*grad[0]; 
  coeff[1] = 0.4898979485566357*val[1]-0.08164965809277261*grad[1]-0.4898979485566357*val[0]-0.08164965809277261*grad[0]; 
  coeff[2] = 0.105409255338946*grad[1]-0.105409255338946*grad[0]; 
  coeff[3] = (-0.05345224838248487*val[1])+0.05345224838248487*grad[1]+0.05345224838248487*val[0]+0.05345224838248487*grad[0];
}

void
gkyl_dg_calc_cubic_2d(const double f[4],
  const double fx[4], const double fy[4], const double fxy[4],
  double *coeff)
{
  coeff[0] = (-0.1666666666666667*fy[3])+0.05555555555555555*fxy[3]-0.1666666666666667*fx[3]+0.5*f[3]+0.1666666666666667*fy[2]-0.05555555555555555*fxy[2]-0.1666666666666667*fx[2]+0.5*f[2]-0.1666666666666667*fy[1]-0.05555555555555555*fxy[1]+0.1666666666666667*fx[1]+0.5*f[1]+0.1666666666666667*fy[0]+0.05555555555555555*fxy[0]+0.1666666666666667*fx[0]+0.5*f[0]; 
  coeff[1] = (-0.1154700538379252*fy[3])+0.01924500897298753*fxy[3]-0.05773502691896259*fx[3]+0.3464101615137755*f[3]+0.1154700538379252*fy[2]-0.01924500897298753*fxy[2]-0.05773502691896259*fx[2]+0.3464101615137755*f[2]+0.1154700538379252*fy[1]+0.01924500897298753*fxy[1]-0.05773502691896259*fx[1]-0.3464101615137755*f[1]-0.1154700538379252*fy[0]-0.01924500897298753*fxy[0]-0.05773502691896259*fx[0]-0.3464101615137755*f[0]; 
  coeff[2] = (-0.05773502691896259*fy[3])+0.01924500897298753*fxy[3]-0.1154700538379252*fx[3]+0.3464101615137755*f[3]-0.05773502691896259*fy[2]+0.01924500897298753*fxy[2]+0.1154700538379252*fx[2]-0.3464101615137755*f[2]-0.05773502691896259*fy[1]-0.01924500897298753*fxy[1]+0.1154700538379252*fx[1]+0.3464101615137755*f[1]-0.05773502691896259*fy[0]-0.01924500897298753*fxy[0]-0.1154700538379252*fx[0]-0.3464101615137755*f[0]; 
  coeff[3] = (-0.04*fy[3])+0.006666666666666667*fxy[3]-0.04*fx[3]+0.24*f[3]-0.04*fy[2]+0.006666666666666667*fxy[2]+0.04*fx[2]-0.24*f[2]+0.04*fy[1]+0.006666666666666667*fxy[1]-0.04*fx[1]-0.24*f[1]+0.04*fy[0]+0.006666666666666667*fxy[0]+0.04*fx[0]+0.24*f[0]; 
  coeff[4] = (-0.02484519974999766*fxy[3])+0.07453559924999298*fx[3]+0.02484519974999766*fxy[2]+0.07453559924999298*fx[2]+0.02484519974999766*fxy[1]-0.07453559924999298*fx[1]-0.02484519974999766*fxy[0]-0.07453559924999298*fx[0]; 
  coeff[5] = 0.07453559924999298*fy[3]-0.02484519974999766*fxy[3]-0.07453559924999298*fy[2]+0.02484519974999766*fxy[2]+0.07453559924999298*fy[1]+0.02484519974999766*fxy[1]-0.07453559924999298*fy[0]-0.02484519974999766*fxy[0]; 
  coeff[6] = (-0.008606629658238702*fxy[3])+0.05163977794943223*fx[3]-0.008606629658238702*fxy[2]-0.05163977794943223*fx[2]+0.008606629658238702*fxy[1]-0.05163977794943223*fx[1]+0.008606629658238702*fxy[0]+0.05163977794943223*fx[0]; 
  coeff[7] = 0.05163977794943223*fy[3]-0.008606629658238702*fxy[3]-0.05163977794943223*fy[2]+0.008606629658238702*fxy[2]-0.05163977794943223*fy[1]-0.008606629658238702*fxy[1]+0.05163977794943223*fy[0]+0.008606629658238702*fxy[0]; 
  coeff[8] = 0.01259881576697424*fy[3]-0.01259881576697424*fxy[3]+0.03779644730092272*fx[3]-0.03779644730092272*f[3]-0.01259881576697424*fy[2]+0.01259881576697424*fxy[2]+0.03779644730092272*fx[2]-0.03779644730092272*f[2]-0.01259881576697424*fy[1]-0.01259881576697424*fxy[1]+0.03779644730092272*fx[1]+0.03779644730092272*f[1]+0.01259881576697424*fy[0]+0.01259881576697424*fxy[0]+0.03779644730092272*fx[0]+0.03779644730092272*f[0]; 
  coeff[9] = 0.03779644730092272*fy[3]-0.01259881576697424*fxy[3]+0.01259881576697424*fx[3]-0.03779644730092272*f[3]+0.03779644730092272*fy[2]-0.01259881576697424*fxy[2]-0.01259881576697424*fx[2]+0.03779644730092272*f[2]+0.03779644730092272*fy[1]+0.01259881576697424*fxy[1]-0.01259881576697424*fx[1]-0.03779644730092272*f[1]+0.03779644730092272*fy[0]+0.01259881576697424*fxy[0]+0.01259881576697424*fx[0]+0.03779644730092272*f[0]; 
  coeff[10] = 0.01111111111111111*fxy[3]-0.01111111111111111*fxy[2]-0.01111111111111111*fxy[1]+0.01111111111111111*fxy[0]; 
  coeff[11] = 0.004364357804719848*fy[3]-0.004364357804719848*fxy[3]+0.02618614682831908*fx[3]-0.02618614682831908*f[3]+0.004364357804719848*fy[2]-0.004364357804719848*fxy[2]-0.02618614682831908*fx[2]+0.02618614682831908*f[2]-0.004364357804719848*fy[1]-0.004364357804719848*fxy[1]+0.02618614682831908*fx[1]+0.02618614682831908*f[1]-0.004364357804719848*fy[0]-0.004364357804719848*fxy[0]-0.02618614682831908*fx[0]-0.02618614682831908*f[0]; 
  coeff[12] = 0.02618614682831908*fy[3]-0.004364357804719848*fxy[3]+0.004364357804719848*fx[3]-0.02618614682831908*f[3]+0.02618614682831908*fy[2]-0.004364357804719848*fxy[2]-0.004364357804719848*fx[2]+0.02618614682831908*f[2]-0.02618614682831908*fy[1]-0.004364357804719848*fxy[1]+0.004364357804719848*fx[1]+0.02618614682831908*f[1]-0.02618614682831908*fy[0]-0.004364357804719848*fxy[0]-0.004364357804719848*fx[0]-0.02618614682831908*f[0]; 
  coeff[13] = (-0.00563436169819011*fy[3])+0.00563436169819011*fxy[3]+0.00563436169819011*fy[2]-0.00563436169819011*fxy[2]+0.00563436169819011*fy[1]+0.00563436169819011*fxy[1]-0.00563436169819011*fy[0]-0.00563436169819011*fxy[0]; 
  coeff[14] = 0.00563436169819011*fxy[3]-0.00563436169819011*fx[3]+0.00563436169819011*fxy[2]+0.00563436169819011*fx[2]-0.00563436169819011*fxy[1]+0.00563436169819011*fx[1]-0.00563436169819011*fxy[0]-0.00563436169819011*fx[0]; 
  coeff[15] = (-0.002857142857142857*fy[3])+0.002857142857142857*fxy[3]-0.002857142857142857*fx[3]+0.002857142857142857*f[3]-0.002857142857142857*fy[2]+0.002857142857142857*fxy[2]+0.002857142857142857*fx[2]-0.002857142857142857*f[2]+0.002857142857142857*fy[1]+0.002857142857142857*fxy[1]-0.002857142857142857*fx[1]-0.002857142857142857*f[1]+0.002857142857142857*fy[0]+0.002857142857142857*fxy[0]+0.002857142857142857*fx[0]+0.002857142857142857*f[0];
}

static inline double
calc_bilinear_grad_xy(double val[4], double dx[2])
{
  return ( (val[3]-val[2])/dx[1] - (val[1]-val[0])/dx[1])/dx[0];
}

gkyl_dg_basis_op_mem *
gkyl_dg_alloc_cubic_1d(int cells)
{
  struct gkyl_dg_basis_op_mem *mem = gkyl_malloc(sizeof(*mem));
  mem->opcode = GKYL_DG_BASIS_OP_CUBIC_1D;
  mem->grad1dx = gkyl_array_new(GKYL_DOUBLE, 1, cells+1);
  return mem;
}

gkyl_dg_basis_op_mem *
gkyl_dg_alloc_cubic_2d(int cells[2])
{
  struct gkyl_dg_basis_op_mem *mem = gkyl_malloc(sizeof(*mem));
  mem->opcode = GKYL_DG_BASIS_OP_CUBIC_2D;

  size_t ncells = (cells[0]+1)*(cells[1]+1);
  mem->grad2dx = gkyl_array_new(GKYL_DOUBLE, 1, ncells);
  mem->grad2dy = gkyl_array_new(GKYL_DOUBLE, 1, ncells);
  mem->grad2dxy = gkyl_array_new(GKYL_DOUBLE, 1, ncells);
  
  return mem;
}

void
gkyl_dg_basis_op_mem_release(gkyl_dg_basis_op_mem *mem)
{
  if (mem->opcode == GKYL_DG_BASIS_OP_CUBIC_1D) {
    gkyl_array_release(mem->grad1dx);
  }
  if (mem->opcode == GKYL_DG_BASIS_OP_CUBIC_2D) {
    gkyl_array_release(mem->grad2dx);
    gkyl_array_release(mem->grad2dy);
    gkyl_array_release(mem->grad2dxy);
  }  
  gkyl_free(mem);
}

void
gkyl_dg_calc_cubic_1d_from_nodal_vals(gkyl_dg_basis_op_mem *mem, int cells, double dx,
  const struct gkyl_array *nodal_vals, struct gkyl_array *cubic)
{
  enum { I, LL, L, R, RR, XE }; // i, i-2, i-1, i+1, i+2 nodes
  
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 1, (int[]) { cells });

  struct gkyl_range nc_range;
  gkyl_range_init_from_shape(&nc_range, 1, (int[]) { cells+1 });

  long offset[XE];
  offset[I] = 0; // i
  offset[LL] = gkyl_range_offset(&nc_range, (int[]) { -2 } ); // i-1  
  offset[L] = gkyl_range_offset(&nc_range, (int[]) { -1 } ); // i-1
  offset[R] = gkyl_range_offset(&nc_range, (int[]) { 1 } ); // i+1
  offset[RR] = gkyl_range_offset(&nc_range, (int[]) { 2 } ); // i+2

  struct gkyl_array *gradx = mem->grad1dx;

  int ilo = nc_range.lower[0], iup = nc_range.upper[0];

  // Step 1: compute gradients at nodes using differencing
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &nc_range); // loop over node range
  while (gkyl_range_iter_next(&iter)) {
    long nidx = gkyl_range_idx(&nc_range, iter.idx);

    if ((iter.idx[0] != ilo) && (iter.idx[0] != iup)) {
      // interior nodes
      const double *val_L = gkyl_array_cfetch(nodal_vals, nidx+offset[L]);
      const double *val_R = gkyl_array_cfetch(nodal_vals, nidx+offset[R]);
      
      double *grad_I = gkyl_array_fetch(gradx, nidx+offset[I]);
      grad_I[0] = (val_R[0]-val_L[0])/(2*dx);
    }

    if (iter.idx[0] == ilo) {
      // left boundary: use second-order one sided differencing
      const double *val_I = gkyl_array_cfetch(nodal_vals, nidx+offset[I]);
      const double *val_R = gkyl_array_cfetch(nodal_vals, nidx+offset[R]);
      const double *val_RR = gkyl_array_cfetch(nodal_vals, nidx+offset[RR]);
      
      double *grad_I = gkyl_array_fetch(gradx, nidx+offset[I]);
      grad_I[0] = -(val_RR[0]-4*val_R[0]+3*val_I[0])/(2*dx);
    }

    if (iter.idx[0] == iup) {
      // right boundary: use second-order one sided differencing
      const double *val_I = gkyl_array_cfetch(nodal_vals, nidx+offset[I]);
      const double *val_L = gkyl_array_cfetch(nodal_vals, nidx+offset[L]);
      const double *val_LL = gkyl_array_cfetch(nodal_vals, nidx+offset[LL]);
      
      double *grad_I = gkyl_array_fetch(gradx, nidx+offset[I]);
      grad_I[0] = (3*val_I[0]-4*val_L[0]+val_LL[0])/(2*dx);
    }
  }

  // Step 2: compute cubic expansions in each cell
  gkyl_range_iter_init(&iter, &range); // loop is over cells
  while (gkyl_range_iter_next(&iter)) {
    
    long nidx = gkyl_range_idx(&nc_range, iter.idx);

    const double *val_I = gkyl_array_cfetch(nodal_vals, nidx+offset[I]);
    const double *val_R = gkyl_array_cfetch(nodal_vals, nidx+offset[R]);
    const double *grad_I = gkyl_array_cfetch(gradx, nidx+offset[I]);
    const double *grad_R = gkyl_array_cfetch(gradx, nidx+offset[R]);

    double val[2] = { val_I[0], val_R[0] };
    double grad[2] = { grad_I[0]*dx/2, grad_R[0]*dx/2 };

    long cidx = gkyl_range_idx(&range, iter.idx);    
    gkyl_dg_calc_cubic_1d(val, grad, gkyl_array_fetch(cubic, cidx));
  }
}

void
gkyl_dg_calc_cubic_2d_from_nodal_vals(gkyl_dg_basis_op_mem *mem, int cells[2], double dx[2],
  const struct gkyl_array *nodal_vals, struct gkyl_array *cubic)
{
  enum {
    I, // (i,j)
    LL, L, R, RR, // (i-2,j) (i-1,j) (i+1,j) (i+2,j)
    BB, B, T, TT, // (i,j-2) (i,j-1) (i,j+1) (i,j+2)
    LT, RT, // (i-1,j+1), (i+1,j+1)
    LB, RB, // (i-1,j-1), (i+1,j-1)
    RRT, RRB, // (i+2,j+1) (i+2,j-1)
    LLT, LLB, // (i-2,j+1) (i-2,j-1)
    LTT, RTT, // (i-1,j+2), (i+1,j+2)
    LBB, RBB, // (i-1,j-2), (i+1,j-2)
    XE
  }; 
  
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, cells);

  struct gkyl_range nc_range;
  gkyl_range_init_from_shape(&nc_range, 2, (int[]) { cells[0]+1, cells[1]+1 });

  long offset[XE];
  offset[I] = 0; // i,j
  offset[LL] = gkyl_range_offset(&nc_range, (int[]) { -2,0 } ); // i-2,j  
  offset[L] = gkyl_range_offset(&nc_range, (int[]) { -1,0 } ); // i-1,j
  offset[R] = gkyl_range_offset(&nc_range, (int[]) { 1,0 } ); // i+1,j
  offset[RR] = gkyl_range_offset(&nc_range, (int[]) { 2,0 } ); // i+2,j
  
  offset[BB] = gkyl_range_offset(&nc_range, (int[]) { 0,-2 } ); // i,j-2  
  offset[B] = gkyl_range_offset(&nc_range, (int[]) { 0,-1 } ); // i,j-1
  offset[T] = gkyl_range_offset(&nc_range, (int[]) { 0,1 } ); // i,j+1
  offset[TT] = gkyl_range_offset(&nc_range, (int[]) { 0,2 } ); // i,j+2

  offset[LT] = gkyl_range_offset(&nc_range, (int[]) { -1,1 } ); // i-1,j+1
  offset[RT] = gkyl_range_offset(&nc_range, (int[]) { 1,1 } ); // i+1,j+1
  offset[LB] = gkyl_range_offset(&nc_range, (int[]) { -1,-1 } ); // i-1,j-1
  offset[RB] = gkyl_range_offset(&nc_range, (int[]) { 1,-1 } ); // i+1,j-1

  offset[RRT] = gkyl_range_offset(&nc_range, (int[]) { 2,1 } ); // i+2,j+1
  offset[RRB] = gkyl_range_offset(&nc_range, (int[]) { 2,-1 } ); // i+2,j-1

  offset[LLT] = gkyl_range_offset(&nc_range, (int[]) { -2,1 } ); // i-2,j+1
  offset[LLB] = gkyl_range_offset(&nc_range, (int[]) { -2,-1 } ); // i-2,j-1

  offset[LTT] = gkyl_range_offset(&nc_range, (int[]) { -1,2 } ); // i-1,j+2
  offset[RTT] = gkyl_range_offset(&nc_range, (int[]) { 1,2 } ); // i+1,j+2

  offset[LBB] = gkyl_range_offset(&nc_range, (int[]) { -1,-2 } ); // i-1,j-2
  offset[RBB] = gkyl_range_offset(&nc_range, (int[]) { 1,-2 } ); // i+1,j-2

  struct gkyl_array *gradx = mem->grad2dx;
  struct gkyl_array *grady = mem->grad2dy;
  struct gkyl_array *gradxy = mem->grad2dxy;

  int ilo = nc_range.lower[0], iup = nc_range.upper[0];
  int jlo = nc_range.lower[1], jup = nc_range.upper[1];

  // Step 1.0: compute gradients at interior nodes using differencing
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &nc_range); // loop over node range
  
  while (gkyl_range_iter_next(&iter)) {
    long nidx = gkyl_range_idx(&nc_range, iter.idx);

    if ((iter.idx[0] != ilo) && (iter.idx[0] != iup) && (iter.idx[1] != jlo) && (iter.idx[1] != jup)) {
      // interior nodes
      const double *val_L = gkyl_array_cfetch(nodal_vals, nidx+offset[L]);
      const double *val_R = gkyl_array_cfetch(nodal_vals, nidx+offset[R]);
      const double *val_T = gkyl_array_cfetch(nodal_vals, nidx+offset[T]);
      const double *val_B = gkyl_array_cfetch(nodal_vals, nidx+offset[B]);

      double *gradx_I = gkyl_array_fetch(gradx, nidx+offset[I]);
      gradx_I[0] = (val_R[0]-val_L[0])/(2*dx[0]);

      double *grady_I = gkyl_array_fetch(grady, nidx+offset[I]);
      grady_I[0] = (val_T[0]-val_B[0])/(2*dx[1]);

      const double *val_LT = gkyl_array_cfetch(nodal_vals, nidx+offset[LT]);
      const double *val_RT = gkyl_array_cfetch(nodal_vals, nidx+offset[RT]);
      const double *val_LB = gkyl_array_cfetch(nodal_vals, nidx+offset[LB]);
      const double *val_RB = gkyl_array_cfetch(nodal_vals, nidx+offset[RB]);

      double *gradxy_I = gkyl_array_fetch(gradxy, nidx+offset[I]);
      gradxy_I[0] = (
        (val_RT[0]-val_RB[0])/(2*dx[1]) - (val_LT[0]-val_LB[0])/(2*dx[1])
      )/(2*dx[0]);
    }
 
    if ((iter.idx[0] == ilo) && (iter.idx[1] != jlo) && (iter.idx[1] != jup)) {
      // left boundary (excluding corners)
      const double *val_I = gkyl_array_cfetch(nodal_vals, nidx+offset[I]);
      const double *val_R = gkyl_array_cfetch(nodal_vals, nidx+offset[R]);
      const double *val_RR = gkyl_array_cfetch(nodal_vals, nidx+offset[RR]);
      const double *val_T = gkyl_array_cfetch(nodal_vals, nidx+offset[T]);
      const double *val_B = gkyl_array_cfetch(nodal_vals, nidx+offset[B]);

      double *gradx_I = gkyl_array_fetch(gradx, nidx+offset[I]);
      gradx_I[0] = -(val_RR[0]-4*val_R[0]+3*val_I[0])/(2*dx[0]);
      
      double *grady_I = gkyl_array_fetch(grady, nidx+offset[I]);
      grady_I[0] = (val_T[0]-val_B[0])/(2*dx[1]);

      const double *val_RT = gkyl_array_cfetch(nodal_vals, nidx+offset[RT]);
      const double *val_RB = gkyl_array_cfetch(nodal_vals, nidx+offset[RB]);
      const double *val_RRT = gkyl_array_cfetch(nodal_vals, nidx+offset[RRT]);
      const double *val_RRB = gkyl_array_cfetch(nodal_vals, nidx+offset[RRB]);

      double gy_RR = (val_RRT[0]-val_RRB[0])/(2*dx[1]);
      double gy_R = (val_RT[0]-val_RB[0])/(2*dx[1]);
      double gy_I = (val_T[0]-val_B[0])/(2*dx[1]);

      double *gradxy_I = gkyl_array_fetch(gradxy, nidx+offset[I]);
      gradxy_I[0] = -(gy_RR-4*gy_R+3*gy_I)/(2*dx[0]);
    }

    if ((iter.idx[0] == iup) && (iter.idx[1] != jlo) && (iter.idx[1] != jup)) {
      // right boundary (excluding corners)
      const double *val_I = gkyl_array_cfetch(nodal_vals, nidx+offset[I]);
      const double *val_L = gkyl_array_cfetch(nodal_vals, nidx+offset[L]);
      const double *val_LL = gkyl_array_cfetch(nodal_vals, nidx+offset[LL]);
      const double *val_T = gkyl_array_cfetch(nodal_vals, nidx+offset[T]);
      const double *val_B = gkyl_array_cfetch(nodal_vals, nidx+offset[B]);

      double *gradx_I = gkyl_array_fetch(gradx, nidx+offset[I]);
      gradx_I[0] = (3*val_I[0]-4*val_L[0]+val_LL[0])/(2*dx[0]);
      
      double *grady_I = gkyl_array_fetch(grady, nidx+offset[I]);
      grady_I[0] = (val_T[0]-val_B[0])/(2*dx[1]);

      const double *val_LT = gkyl_array_cfetch(nodal_vals, nidx+offset[LT]);
      const double *val_LB = gkyl_array_cfetch(nodal_vals, nidx+offset[LB]);
      const double *val_LLT = gkyl_array_cfetch(nodal_vals, nidx+offset[LLT]);
      const double *val_LLB = gkyl_array_cfetch(nodal_vals, nidx+offset[LLB]);

      double gy_LL = (val_LLT[0]-val_LLB[0])/(2*dx[1]);
      double gy_L = (val_LT[0]-val_LB[0])/(2*dx[1]);
      double gy_I = (val_T[0]-val_B[0])/(2*dx[1]);

      double *gradxy_I = gkyl_array_fetch(gradxy, nidx+offset[I]);
      gradxy_I[0] = (3*gy_I-4*gy_L+gy_LL)/(2*dx[0]);
    }

    if ((iter.idx[0] != ilo) && (iter.idx[0] != iup) && (iter.idx[1] == jlo)) {
      // bottom boundary (excluding corners)
      const double *val_I = gkyl_array_cfetch(nodal_vals, nidx+offset[I]);
      const double *val_T = gkyl_array_cfetch(nodal_vals, nidx+offset[T]);
      const double *val_TT = gkyl_array_cfetch(nodal_vals, nidx+offset[TT]);
      const double *val_L = gkyl_array_cfetch(nodal_vals, nidx+offset[L]);
      const double *val_R = gkyl_array_cfetch(nodal_vals, nidx+offset[R]);

      double *gradx_I = gkyl_array_fetch(gradx, nidx+offset[I]);
      gradx_I[0] = (val_R[0]-val_L[0])/(2*dx[0]);
      
      double *grady_I = gkyl_array_fetch(grady, nidx+offset[I]);
      grady_I[0] = -(val_TT[0]-4*val_T[0]+3*val_I[0])/(2*dx[1]);

      const double *val_LT = gkyl_array_cfetch(nodal_vals, nidx+offset[LT]);
      const double *val_RT = gkyl_array_cfetch(nodal_vals, nidx+offset[RT]);
      const double *val_LTT = gkyl_array_cfetch(nodal_vals, nidx+offset[LTT]);
      const double *val_RTT = gkyl_array_cfetch(nodal_vals, nidx+offset[RTT]);

      double gx_TT = (val_RTT[0]-val_LTT[0])/(2*dx[0]);
      double gx_T = (val_RT[0]-val_LT[0])/(2*dx[0]);
      double gx_I = (val_R[0]-val_L[0])/(2*dx[0]);

      double *gradxy_I = gkyl_array_fetch(gradxy, nidx+offset[I]);
      gradxy_I[0] = -(gx_TT-4*gx_T+3*gx_I)/(2*dx[1]);
    }

    if ((iter.idx[0] != ilo) && (iter.idx[0] != iup) && (iter.idx[1] == jup)) {
      // top boundary (excluding corners)
      const double *val_I = gkyl_array_cfetch(nodal_vals, nidx+offset[I]);
      const double *val_B = gkyl_array_cfetch(nodal_vals, nidx+offset[B]);
      const double *val_BB = gkyl_array_cfetch(nodal_vals, nidx+offset[BB]);
      const double *val_L = gkyl_array_cfetch(nodal_vals, nidx+offset[L]);
      const double *val_R = gkyl_array_cfetch(nodal_vals, nidx+offset[R]);

      double *gradx_I = gkyl_array_fetch(gradx, nidx+offset[I]);
      gradx_I[0] = (val_R[0]-val_L[0])/(2*dx[0]);
      
      double *grady_I = gkyl_array_fetch(grady, nidx+offset[I]);
      grady_I[0] = (3*val_I[0]-4*val_B[0]+val_BB[0])/(2*dx[1]);

      const double *val_LB = gkyl_array_cfetch(nodal_vals, nidx+offset[LB]);
      const double *val_RB = gkyl_array_cfetch(nodal_vals, nidx+offset[RB]);
      const double *val_LBB = gkyl_array_cfetch(nodal_vals, nidx+offset[LBB]);
      const double *val_RBB = gkyl_array_cfetch(nodal_vals, nidx+offset[RBB]);

      double gx_BB = (val_RBB[0]-val_LBB[0])/(2*dx[0]);
      double gx_B = (val_RB[0]-val_LB[0])/(2*dx[0]);
      double gx_I = (val_R[0]-val_L[0])/(2*dx[0]);

      double *gradxy_I = gkyl_array_fetch(gradxy, nidx+offset[I]);
      gradxy_I[0] = (3*gx_I-4*gx_B+gx_BB)/(2*dx[1]);
    }

    if ((iter.idx[0] == ilo) && (iter.idx[1] == jlo)) {
      // lower-left corner
      const double *val_I = gkyl_array_cfetch(nodal_vals, nidx+offset[I]);
      const double *val_R = gkyl_array_cfetch(nodal_vals, nidx+offset[R]);
      const double *val_RR = gkyl_array_cfetch(nodal_vals, nidx+offset[RR]);
      const double *val_T = gkyl_array_cfetch(nodal_vals, nidx+offset[T]);
      const double *val_TT = gkyl_array_cfetch(nodal_vals, nidx+offset[TT]);
      const double *val_RT = gkyl_array_cfetch(nodal_vals, nidx+offset[RT]);

      double *gradx_I = gkyl_array_fetch(gradx, nidx+offset[I]);
      gradx_I[0] = -(val_RR[0]-4*val_R[0]+3*val_I[0])/(2*dx[0]);
      
      double *grady_I = gkyl_array_fetch(grady, nidx+offset[I]);
      grady_I[0] = -(val_TT[0]-4*val_T[0]+3*val_I[0])/(2*dx[1]);

      double *gradxy_I = gkyl_array_fetch(gradxy, nidx+offset[I]);
      double vxy[4] = { val_I[0], val_T[0], val_R[0], val_RT[0] };
      gradxy_I[0] = calc_bilinear_grad_xy(vxy, dx);
    }

    if ((iter.idx[0] == ilo) && (iter.idx[1] == jup)) {
      // upper-left corner
      const double *val_I = gkyl_array_cfetch(nodal_vals, nidx+offset[I]);
      const double *val_R = gkyl_array_cfetch(nodal_vals, nidx+offset[R]);
      const double *val_RR = gkyl_array_cfetch(nodal_vals, nidx+offset[RR]);
      const double *val_B = gkyl_array_cfetch(nodal_vals, nidx+offset[B]);
      const double *val_BB = gkyl_array_cfetch(nodal_vals, nidx+offset[BB]);
      const double *val_RB = gkyl_array_cfetch(nodal_vals, nidx+offset[RB]);

      double *gradx_I = gkyl_array_fetch(gradx, nidx+offset[I]);
      gradx_I[0] = -(val_RR[0]-4*val_R[0]+3*val_I[0])/(2*dx[0]);
      
      double *grady_I = gkyl_array_fetch(grady, nidx+offset[I]);
      grady_I[0] = (3*val_I[0]-4*val_B[0]+val_BB[0])/(2*dx[1]);

      double *gradxy_I = gkyl_array_fetch(gradxy, nidx+offset[I]);
      double vxy[4] = { val_B[0], val_I[0], val_RB[0], val_R[0] };
      gradxy_I[0] = calc_bilinear_grad_xy(vxy, dx);
    }

    if ((iter.idx[0] == iup) && (iter.idx[1] == jlo)) {
      // lower-right corner
      const double *val_I = gkyl_array_cfetch(nodal_vals, nidx+offset[I]);
      const double *val_L = gkyl_array_cfetch(nodal_vals, nidx+offset[L]);
      const double *val_LL = gkyl_array_cfetch(nodal_vals, nidx+offset[LL]);
      const double *val_T = gkyl_array_cfetch(nodal_vals, nidx+offset[T]);
      const double *val_TT = gkyl_array_cfetch(nodal_vals, nidx+offset[TT]);
      const double *val_LT = gkyl_array_cfetch(nodal_vals, nidx+offset[LT]);

      double *gradx_I = gkyl_array_fetch(gradx, nidx+offset[I]);
      gradx_I[0] = (3*val_I[0]-4*val_L[0]+val_LL[0])/(2*dx[0]);
      
      double *grady_I = gkyl_array_fetch(grady, nidx+offset[I]);
      grady_I[0] = -(val_TT[0]-4*val_T[0]+3*val_I[0])/(2*dx[1]);

      double *gradxy_I = gkyl_array_fetch(gradxy, nidx+offset[I]);
      double vxy[4] = { val_L[0], val_LT[0], val_I[0], val_T[0] };
      gradxy_I[0] = calc_bilinear_grad_xy(vxy, dx);
    }


    if ((iter.idx[0] == iup) && (iter.idx[1] == jup)) {
      // upper-right corner
      const double *val_I = gkyl_array_cfetch(nodal_vals, nidx+offset[I]);
      const double *val_L = gkyl_array_cfetch(nodal_vals, nidx+offset[L]);
      const double *val_LL = gkyl_array_cfetch(nodal_vals, nidx+offset[LL]);
      const double *val_B = gkyl_array_cfetch(nodal_vals, nidx+offset[B]);
      const double *val_BB = gkyl_array_cfetch(nodal_vals, nidx+offset[BB]);
      const double *val_LB = gkyl_array_cfetch(nodal_vals, nidx+offset[LB]);

      double *gradx_I = gkyl_array_fetch(gradx, nidx+offset[I]);
      gradx_I[0] = (3*val_I[0]-4*val_L[0]+val_LL[0])/(2*dx[0]);
      
      double *grady_I = gkyl_array_fetch(grady, nidx+offset[I]);
      grady_I[0] = (3*val_I[0]-4*val_B[0]+val_BB[0])/(2*dx[1]);

      double *gradxy_I = gkyl_array_fetch(gradxy, nidx+offset[I]);
      double vxy[4] = { val_LB[0], val_L[0], val_B[0], val_I[0] };
      gradxy_I[0] = calc_bilinear_grad_xy(vxy, dx);
    }
  }

  // Step 2: compute cubic expansions in each cell
  gkyl_range_iter_init(&iter, &range); // loop is over cells
  while (gkyl_range_iter_next(&iter)) {
    
    long nidx = gkyl_range_idx(&nc_range, iter.idx);

    // four corners nodes of cell
    const double *val_I = gkyl_array_cfetch(nodal_vals, nidx+offset[I]);
    const double *val_R = gkyl_array_cfetch(nodal_vals, nidx+offset[R]);
    const double *val_T = gkyl_array_cfetch(nodal_vals, nidx+offset[T]);
    const double *val_RT = gkyl_array_cfetch(nodal_vals, nidx+offset[RT]);
    
    const double *gradx_I = gkyl_array_cfetch(gradx, nidx+offset[I]);
    const double *gradx_R = gkyl_array_cfetch(gradx, nidx+offset[R]);
    const double *gradx_T = gkyl_array_cfetch(gradx, nidx+offset[T]);
    const double *gradx_RT = gkyl_array_cfetch(gradx, nidx+offset[RT]);

    const double *grady_I = gkyl_array_cfetch(grady, nidx+offset[I]);
    const double *grady_R = gkyl_array_cfetch(grady, nidx+offset[R]);
    const double *grady_T = gkyl_array_cfetch(grady, nidx+offset[T]);
    const double *grady_RT = gkyl_array_cfetch(grady, nidx+offset[RT]);

    const double *gradxy_I = gkyl_array_cfetch(gradxy, nidx+offset[I]);
    const double *gradxy_R = gkyl_array_cfetch(gradxy, nidx+offset[R]);
    const double *gradxy_T = gkyl_array_cfetch(gradxy, nidx+offset[T]);
    const double *gradxy_RT = gkyl_array_cfetch(gradxy, nidx+offset[RT]);

    double val[4] = { val_I[0], val_T[0], val_R[0], val_RT[0] };
    double gradx[4] = { gradx_I[0]*dx[0]/2, gradx_T[0]*dx[0]/2, gradx_R[0]*dx[0]/2, gradx_RT[0]*dx[0]/2 };
    double grady[4] = { grady_I[0]*dx[1]/2, grady_T[0]*dx[1]/2, grady_R[0]*dx[1]/2, grady_RT[0]*dx[1]/2 };
    double gradxy[4] = { gradxy_I[0]*dx[0]/2*dx[1]/2, gradxy_T[0]*dx[0]/2*dx[1]/2, gradxy_R[0]*dx[0]/2*dx[1]/2, gradxy_RT[0]*dx[0]/2*dx[1]/2 };

    long cidx = gkyl_range_idx(&range, iter.idx);
    double *coeff = gkyl_array_fetch(cubic, cidx);
    gkyl_dg_calc_cubic_2d(val, gradx, grady, gradxy, coeff);
  }
}


static double
eval_laplacian_expand_2d_tensor_p3(int dir, const double *z, const double *f )
{
  const double z0 = z[0];
  const double z1 = z[1];
  if (dir == 0)
    return 131.25*f[15]*z0*z1*z1*z1+22.18529918662356*f[14]*z1*z1*z1+66.55589755987069*f[13]*z0*z1*z1+11.25*f[10]*z1*z1-78.75*f[15]*z0*z1+ 34.3693177121688*f[11]*z0*z1-13.31117951197414*f[14]*z1+5.809475019311125*f[6]*z1-22.18529918662356*f[13]*z0+ 19.84313483298443*f[8]*z0-3.75*f[10]+3.354101966249685*f[4];

  if (dir == 1)
    return 131.25*f[15]*z0*z0*z0*z1+66.55589755987069*f[14]*z0*z0*z1-78.75*f[15]*z0*z1+34.3693177121688*f[12]*z0*z1-22.18529918662356*f[14]*z1+19.84313483298443*f[9]*z1+22.18529918662356*f[13]*z0*z0*z0+11.25*f[10]*z0*z0-13.31117951197414*f[13]*z0+5.809475019311125*f[7]*z0-3.75*f[10]+3.354101966249685*f[5];

  return 0.0; // can't happen, suppresses warning

}

static double
eval_mixedpartial_expand_2d_tensor_p3(const double *z, const double *f )
{
  const double z0 = z[0];
  const double z1 = z[1];
  return 196.875*f[15]*z0*z0*z1*z1+66.55589755987069*f[14]*z0*z1*z1-39.375*f[15]*z1*z1+17.1846588560844*f[12]*z1*z1+66.55589755987069*f[13]*z0*z0*z1+22.5*f[10]*z0*z1-13.31117951197414*f[13]*z1+5.809475019311125*f[7]*z1-39.375*f[15]*z0*z0+17.1846588560844*f[11]*z0*z0-13.31117951197414*f[14]*z0+5.809475019311125*f[6]*z0+7.875*f[15]-3.43693177121688*f[12]-3.43693177121688*f[11]+1.5*f[3];

}

static void
evalf_free(const struct gkyl_ref_count* rc)
{
  struct gkyl_basis_ops_evalf *evf = container_of(rc, struct gkyl_basis_ops_evalf, ref_count);

  struct dg_basis_ops_evalf_ctx *ctx = evf->ctx;
  gkyl_array_release(ctx->cubic);
  gkyl_free(ctx);
  gkyl_free(evf);
}

// function for computing cubic at a specified coordinate
static void
eval_cubic(double t, const double *xn, double *fout, void *ctx)
{
  struct dg_basis_ops_evalf_ctx *ectx = ctx;
  
  int idx[GKYL_MAX_DIM];  
  gkyl_rect_grid_coord_idx(&ectx->grid, xn, idx);
  for (int d=0; d<ectx->ndim; ++d) {
    idx[d] = GKYL_MIN2(ectx->local.upper[d], idx[d]);
    idx[d] = GKYL_MAX2(ectx->local.lower[d], idx[d]);
  }

  double xc[GKYL_MAX_DIM];
  gkyl_rect_grid_cell_center(&ectx->grid, idx, xc);

  double eta[GKYL_MAX_DIM];
  for (int d=0; d<ectx->ndim; ++d)
    eta[d] = 2.0*(xn[d]-xc[d])/ectx->grid.dx[d];
  
  long lidx = gkyl_range_idx(&ectx->local, idx);
  const double *fdg = gkyl_array_cfetch(ectx->cubic, lidx);
  
  fout[0] = ectx->basis.eval_expand(eta, fdg);
}

// function for computing cubic at a specified coordinate
static void
eval_cubic_wgrad(double t, const double *xn, double *fout, void *ctx)
{
  struct dg_basis_ops_evalf_ctx *ectx = ctx;
  
  int idx[GKYL_MAX_DIM];  
  gkyl_rect_grid_coord_idx(&ectx->grid, xn, idx);
  for (int d=0; d<ectx->ndim; ++d) {
    idx[d] = GKYL_MIN2(ectx->local.upper[d], idx[d]);
    idx[d] = GKYL_MAX2(ectx->local.lower[d], idx[d]);
  }

  double xc[GKYL_MAX_DIM];
  gkyl_rect_grid_cell_center(&ectx->grid, idx, xc);

  double eta[GKYL_MAX_DIM];
  for (int d=0; d<ectx->ndim; ++d)
    eta[d] = 2.0*(xn[d]-xc[d])/ectx->grid.dx[d];
  
  long lidx = gkyl_range_idx(&ectx->local, idx);
  const double *fdg = gkyl_array_cfetch(ectx->cubic, lidx);
  
  fout[0] = ectx->basis.eval_expand(eta, fdg);
  fout[1] = ectx->basis.eval_grad_expand(0, eta, fdg)*2/ectx->grid.dx[0];
  if (ectx->ndim > 1)
    fout[2] = ectx->basis.eval_grad_expand(1, eta, fdg)*2/ectx->grid.dx[1];
}

// function for computing cubic at a specified coordinate
static void
eval_cubic_wgrad2(double t, const double *xn, double *fout, void *ctx)
{
  struct dg_basis_ops_evalf_ctx *ectx = ctx;

  int idx[GKYL_MAX_DIM];
  gkyl_rect_grid_coord_idx(&ectx->grid, xn, idx);
  for (int d=0; d<ectx->ndim; ++d) {
    idx[d] = GKYL_MIN2(ectx->local.upper[d], idx[d]);
    idx[d] = GKYL_MAX2(ectx->local.lower[d], idx[d]);
  }

  double xc[GKYL_MAX_DIM];
  gkyl_rect_grid_cell_center(&ectx->grid, idx, xc);

  double eta[GKYL_MAX_DIM];
  for (int d=0; d<ectx->ndim; ++d)
    eta[d] = 2.0*(xn[d]-xc[d])/ectx->grid.dx[d];

  long lidx = gkyl_range_idx(&ectx->local, idx);
  const double *fdg = gkyl_array_cfetch(ectx->cubic, lidx);

  fout[0] = ectx->basis.eval_expand(eta, fdg);
  if (ectx->ndim > 1) {
    fout[1] = eval_laplacian_expand_2d_tensor_p3(0, eta, fdg)*2/ectx->grid.dx[0]*2/ectx->grid.dx[0];
    fout[2] = eval_laplacian_expand_2d_tensor_p3(1, eta, fdg)*2/ectx->grid.dx[1]*2/ectx->grid.dx[1];
    fout[3] = eval_mixedpartial_expand_2d_tensor_p3(eta, fdg)*2/ectx->grid.dx[0]*2/ectx->grid.dx[1];
  }
}


struct gkyl_basis_ops_evalf*
gkyl_dg_basis_ops_evalf_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_array *nodal_vals)
{
  if (grid->ndim > 2) return 0;

  struct dg_basis_ops_evalf_ctx *ctx = gkyl_malloc(sizeof(*ctx));
  int ndim = ctx->ndim = grid->ndim;
  
  int cells[2];
  double dx[2];
  size_t vol = 1;  
  for (int d=0; d<ndim; ++d) {
    vol *= grid->cells[d];
    dx[d] = ctx->dx[d] = grid->dx[d];    
    cells[d] = ctx->cells[d] = grid->cells[d];
  }

  ctx->grid = *grid;
  int nghost[GKYL_MAX_CDIM] = { 0 };
  gkyl_create_grid_ranges(grid, nghost, &ctx->local_ext, &ctx->local);

  gkyl_cart_modal_tensor(&ctx->basis, ndim, 3);
  ctx->cubic = gkyl_array_new(GKYL_DOUBLE, ctx->basis.num_basis, vol);

  gkyl_dg_basis_op_mem *mem = 0;
  if (ndim == 1) {
    mem = gkyl_dg_alloc_cubic_1d(cells[0]);
    gkyl_dg_calc_cubic_1d_from_nodal_vals(mem, grid->cells[0], dx[0], nodal_vals, ctx->cubic);
  }
  if (ndim == 2) {
    mem = gkyl_dg_alloc_cubic_2d(cells);
    gkyl_dg_calc_cubic_2d_from_nodal_vals(mem, cells, dx, nodal_vals, ctx->cubic);
  }
  gkyl_dg_basis_op_mem_release(mem);

  struct gkyl_basis_ops_evalf *evf = gkyl_malloc(sizeof(*evf));
  evf->ctx = ctx;
  evf->eval_cubic = eval_cubic;
  evf->eval_cubic_wgrad = eval_cubic_wgrad;
  evf->eval_cubic_wgrad2 = eval_cubic_wgrad2;
  evf->eval_cubic_laplacian = eval_laplacian_expand_2d_tensor_p3;
  evf->eval_cubic_mixedpartial = eval_mixedpartial_expand_2d_tensor_p3;
  evf->ref_count = (struct gkyl_ref_count) { evalf_free, 1 };
  
  return evf;
}

bool
gkyl_dg_basis_ops_evalf_write_cubic(const struct gkyl_basis_ops_evalf *evf, const char *fname)
{
  struct dg_basis_ops_evalf_ctx *ectx = evf->ctx;

  struct gkyl_msgpack_data *mdata = gkyl_msgpack_create(2,
    (struct gkyl_msgpack_map_elem []) {
      { .key = "polyOrder", .elem_type = GKYL_MP_INT, .ival = 3 },
      { .key = "basisType", .elem_type = GKYL_MP_STRING, .cval = ectx->basis.id }
    }
  );

  enum gkyl_array_rio_status status =
    gkyl_grid_sub_array_write(&ectx->grid, &ectx->local, mdata, ectx->cubic, fname);
  gkyl_msgpack_data_release(mdata);  
  
  return status == GKYL_ARRAY_RIO_SUCCESS;
}

struct gkyl_basis_ops_evalf *
gkyl_dg_basis_ops_evalf_acquire(const struct gkyl_basis_ops_evalf *evf)
{
  gkyl_ref_count_inc(&evf->ref_count);
  return (struct gkyl_basis_ops_evalf*) evf;
}

void
gkyl_dg_basis_ops_evalf_release(struct gkyl_basis_ops_evalf *evf)
{
  gkyl_ref_count_dec(&evf->ref_count);
}
