#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_mat.h>
#include <gkyl_twistshift_dg.h>
#include <gkyl_twistshift_dg_priv.h>
#include <gkyl_util.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>

// Notes:
//   a) Hard-coded parameters:
//     - wrap_to_range: eps.
//     - find_donors: delta_frac, num_test_pt.
//     - find_intersect: tol, max_iter, num_steps.
//     - calc_mats: shift_dir_idx_tar.
//     - tol_xi: Minimum allowed spacing between the lower and
//               upper xi (logical x) limits of subcell integral.
//   b) Unlike the procedures described in M. Francisquez, et al. CPC 298
//   (2024) 109109, all subcell integrals are now done with variable y limits.
//   This is possible once we realize that figure 4 is not drawn accurately;
//   the blue lines should be separated by Delta y at all points.
//   c) This updater only works on 5D distributions. Likely only minor changes
//   are needed to make it work in other dimensions.
//   d) 99% of the code is written to support a BC, a shift and shear in any
//   direction. Maybe the only thing that needs to change is the permutted
//   range and its use.
//
// List of functions used in computing sub-cell integrals (scimat).
//   - ts_grid_cell_boundary_in_dir: cell boundary coordinate in given dir.
//   - ts_grid_cell_boundaries: get all cell boundary coords.
//   - ts_p2l: physical to logical transform.
//   - ts_interval_dx_and_xc: compute length and center of an interval.
//   - ts_grid_length_in_dir: length of the grid in given dir.
//   - ts_wrap_to_range: wrap a number to a range assuming periodicity.
//   - ts_shift_dir_idx_do_linidx: linear index to first donor of a given target
//                                cell in shift_dir_idx_do.
//   - ts_check_shifted_test_point: evaluate a shifted point's cell as a
//                                  potential donor cell.
//   - ts_find_donors: find and record the donor cells for each target.
//   - ts_root_find: Finds the root of a given function.
//   - ts_shifted_coord_loss_func: Loss function used to find where yTar-S
//                                intersects yDo.
//   - ts_sign: return the sign of a double.
//   - ts_ts_donor_target_offset: offset between donor and target cells.
//   - ts_find_intersect: finds the intersection of yTar-S and yDo.
//   - ts_comp_to_phys: transform a computational to a physical coord.
//   - ts_nod2mod_proj_1d: evaluate a 1D function at nodes and do a n2m transform
//                         to get the coefficients of the DG representation.
//   - ts_integral_xlimdg: subcell integral with variable x limits.
//   - ts_integral_ylimdg: subcell integral with variable y limits.
//   - ts_integral_fullcelllimdg: integral over the whole cell.
//   - ts_one: return 1 (for projections).
//   - ts_minus_one: return -1 (for projections).
//   - ts_shift_coord_shifted_log: coordinate in shift_dir shifted and transformed
//                                 to logical space.
//   - ts_subcellint_sNi_sNii: subcell integral sNi or sNii.
//   - ts_subcellint_si_sii: subcell integral si or sii.
//   - ts_subcellint_siii_siv: subcell integral siii or siv.
//   - ts_subcellint_sv_svi: subcell integral sv or svi.
//   - ts_subcellint_svii_sviii: subcell integral svii or sviii.
//   - ts_subcellint_six_sx: subcell integral six or sx.
//   - ts_subcellint_sxi_sxii: subcell integral sxi or sxii.
//   - ts_subcellint_sxiii_sxiv: subcell integral sxiii or sxiv.
//   - ts_subcellint_sxv_sxvi: subcell integral sxv or sxvi.
//   - ts_calc_mats: create scimat with the result of the subcell integrals.
//
//  Two additional helper functions:
//   - ts_calc_num_numcol_fidx_do: index map to populate fmat with donors.
//   - ts_calc_num_numcol_fidx_tar: index map to get mat-mat mult results.

// Option to use the user-provided function describing the shift
// or a DG representation of it:
//   = 0 DG representation (default and preferred).
//   = 1 user-provided shift function.
// Note: the kernels that ultimately perform the integrals
// always use the DG representation.
#define shift_func_op 0

// Minimum allowed spacing between the lower and
// upper xi (logical x) limits of subcell integral.
#define tol_xi 1.0e-15

// Indices in 4-element cell boundary array.
#define cellb_lo(dir) (2*dir)
#define cellb_up(dir) (2*dir+1)

double
ts_grid_cell_boundary_in_dir(struct gkyl_rect_grid *grid, const int *idx, enum gkyl_edge_loc edge, int dir)
{
  // Get the coordinate of the cell boundary in specified direction.
  double xc[grid->ndim];
  gkyl_rect_grid_cell_center(grid, idx, xc);
  return edge == GKYL_LOWER_EDGE? xc[dir]-0.5*grid->dx[dir] : xc[dir]+0.5*grid->dx[dir];
}

void
ts_grid_cell_boundaries(struct gkyl_rect_grid *grid, const int *idx, double *cell_bounds)
{
  // Get the cell boundaries in every dimension. The array cell_bounds
  // must be a 2*grid->ndim array.
  for (int d=0; d<grid->ndim; d++) {
    cell_bounds[d*2] = ts_grid_cell_boundary_in_dir(grid, idx, GKYL_LOWER_EDGE, d);
    cell_bounds[d*2+1] = ts_grid_cell_boundary_in_dir(grid, idx, GKYL_UPPER_EDGE, d);
  }
}

static inline double
ts_p2l(double coord, double cell_center, double dx)
{
  // Transform a physical coordinate (coord) to the [-1,1] logical
  // space in a cell centered at cell_center and with length dx. 
  return 2.0*(coord - cell_center)/dx;
}

// Evaluation of the shift through the DG representation.
static inline void
ts_shift_dg_eval(double t, const double *coord, double *fout, void *ctx)
{
  struct ts_shift_dg_eval_ctx *tsectx = ctx;

  int cell_idx[GKYL_MAX_DIM];
  gkyl_rect_grid_coord_idx(tsectx->shear_grid, coord, cell_idx);
  // Ensure that we do not go outside of the range 
  // (it does sometimes if x=x_max,x_min).
  cell_idx[0] = fmin(cell_idx[0], tsectx->shear_r->upper[0]);
  cell_idx[0] = fmax(cell_idx[0], tsectx->shear_r->lower[0]);

  double xc[GKYL_MAX_DIM];
  gkyl_rect_grid_cell_center(tsectx->shear_grid, cell_idx, xc);

  long shift_loc = gkyl_range_idx(tsectx->shear_r, cell_idx);
  double *shift_c = (double *) gkyl_array_fetch(tsectx->shift_dg, shift_loc);
  double xp = ts_p2l(coord[0], xc[0], tsectx->shear_grid->dx[0]);

  fout[0] = tsectx->shift_b->eval_expand(&(double) {xp}, shift_c);
}

void
ts_interval_dx_and_xc(const double *interval, double *dx, double *xc)
{
  // Compute the lenth (dx) and center (xc) of [interval[0], interval[1]].
  double lo = interval[0], up = interval[1];
  dx[0] = up - lo;
  xc[0] = 0.5*(up + lo);
}

static inline double
ts_grid_length_in_dir(struct gkyl_rect_grid *grid, int dir)
{
  return grid->upper[dir] - grid->lower[dir];
}

double
ts_wrap_to_range(double val, double lower, double upper, bool pick_upper)
{
  // Wrap a number to range [lower,upper]. If pickUpper=true, output upper when
  // val is a multiple of upper. Otherwise multiples of upper wrap to lower.
  double L        = upper - lower;
  double disp     = fmod(val - lower, L);
  double vwrapped = lower + fmod(L + disp, L);
  double eps      = 1.e-12;
  if ( (lower-eps < vwrapped && vwrapped < lower + eps) ||
       (upper-eps < vwrapped && vwrapped < upper + eps) ) {
    if (pick_upper) 
      return upper;
    else
      return lower;
  }
  else
    return vwrapped;
}

long
ts_shift_dir_idx_do_linidx(const int *num_do, int shear_dir_idx, int shift_dir_idx,
  int shift_dir_num_cells, int shear_r_lower)
{
  // Return the linear index to the first donor for the idx=(i,j) target cell,
  // in the shift_dir_idx_do array. We assume shift_dir_idx_do (whose dimensions
  // are Nx,Ny,num_do(i)) is in row-major order, and that it has num_do donors
  // at each cell in the shear_dir_in_ts_grid direction.
  long linc = 0;
  // Count the number of donors in cells with an idx in the shear dir lower
  // than this one. NOTE: the -1 here is because the idx is often 1-index
  // (since ghost cells are the 0th index) but num_do is only defined on the
  // local range. 
  for (int i=0; i<shear_dir_idx-shear_r_lower; i++)
    linc += num_do[i] * shift_dir_num_cells;

  // Count the number of donors in cells with the same idx in the shear dir.
  // NOTE: the -1 here is because the idx is often 1-index
  // (since ghost cells are the 0th index) but num_do is only defined on the
  // local range.
  linc += num_do[shear_dir_idx-shear_r_lower] * (shift_dir_idx-1);

  return linc;
}

void
ts_check_shifted_test_point(struct gkyl_twistshift_dg *up, const double *test_pt, const double *xc,
  const double *dx, double *shift_c, const int *idx, bool pick_lower, int *num_do_curr,
  gkyl_mem_buff shift_dir_idx_do_buff)
{
  // Shift the test point, find the cell that contains it, and if we
  // haven't included it yet, add it to our list of donors.

  int shear_idx[] = {idx[up->shear_dir_in_ts_grid]};
  int shift_idx[] = {idx[up->shift_dir_in_ts_grid]};

  int *shift_dir_idx_do_buff_ptr = (int *) gkyl_mem_buff_data(shift_dir_idx_do_buff);

  // Evaluate the shift at this test point.
  double test_pt_in_shear_dir = test_pt[up->shear_dir_in_ts_grid];
  double xc_in_shear_dir = xc[up->shear_dir_in_ts_grid];
  double dx_in_shear_dir = dx[up->shear_dir_in_ts_grid];
  double shift_at_pt = up->shift_b.eval_expand(
    &(double) {ts_p2l(test_pt_in_shear_dir, xc_in_shear_dir, dx_in_shear_dir)}, shift_c);

  // Find the index of the cell that owns the shifted point.
  double shifted_test_pt[] = { ts_wrap_to_range(test_pt[up->shift_dir_in_ts_grid] - shift_at_pt,
    up->ts_grid.lower[up->shift_dir_in_ts_grid], up->ts_grid.upper[up->shift_dir_in_ts_grid],
    false) }; // Shifted test point.
  int shift_dir_idx_test_pt[1];
  gkyl_rect_grid_find_cell(&up->shift_grid, shifted_test_pt, (bool[]) {pick_lower}, (int[]) {-1}, shift_dir_idx_test_pt);

  // Get the linear index to the list of donors for this target.
  long linidx = ts_shift_dir_idx_do_linidx(up->num_do,
    shear_idx[0], shift_idx[0], up->ts_grid.cells[up->shift_dir_in_ts_grid], up->shear_r.lower[0]);
  // If this donor is not in our list of donors, include it.
  bool donor_not_found = true;
  for (int k=0; k<num_do_curr[0]; k++) {
    if (shift_dir_idx_do_buff_ptr[linidx+k] == shift_dir_idx_test_pt[0]) {
      donor_not_found = false;
      break;
    }
  }
  if (donor_not_found) {
    if (num_do_curr[0] > 0) {
      // Insert one more int. Only if num_do_curr>0 because we already
      // allocated space for the first donor.
      size_t new_buff_sz = gkyl_mem_buff_size(shift_dir_idx_do_buff) + sizeof(int);
      shift_dir_idx_do_buff = gkyl_mem_buff_resize(shift_dir_idx_do_buff, new_buff_sz);
    }

    // Get the pointer again in case it changed.
    shift_dir_idx_do_buff_ptr = (int *) gkyl_mem_buff_data(shift_dir_idx_do_buff);
    shift_dir_idx_do_buff_ptr[linidx+num_do_curr[0]] = shift_dir_idx_test_pt[0];

    num_do_curr[0] += 1;
  }
}

void
ts_find_donors(struct gkyl_twistshift_dg *up)
{
  // Find the donor cells for each target cell in the TS grid.

  double delta_frac = 1.e-4; // Distance away from the boundary, as fraction of cell length.
  // Must be larger than the eps=1e-6 tolerance in is_in_cell (used by
  // gkyl_rect_grid_find_cell) so that shifted test points are never
  // ambiguously on a cell boundary.
  int num_test_pt[2] = {10, 10}; // Number of test points taken along each side of the cell.

  double step_sz[2] = {0.0}; // Size of the step between test points.
  double delta[2] = {0.0}; // Space between cell boundary and test points.
  for (int d=0; d<2; d++) {
    delta[d] = delta_frac*up->ts_grid.dx[d];
    step_sz[d] = (up->ts_grid.dx[d] - 2.0*delta[d])/(num_test_pt[d]-1);
  }

  // Number of donors at each cell of the shear direction.
  up->num_do = (int*) gkyl_malloc(up->shear_r.volume * sizeof(int));
  for (int i=0; i<up->shear_r.volume; i++)
    up->num_do[i] = -1;

  // Temporary buffer to store donors at (resized below).
  size_t curr_buff_sz = up->ts_r.volume * sizeof(int);
  gkyl_mem_buff shift_dir_idx_do_buff = gkyl_mem_buff_new(curr_buff_sz);

  int idx[] = {up->shear_r.lower[0]};
  long linidx = gkyl_range_idx(&up->shear_r, idx);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->ts_r);
  while (gkyl_range_iter_next(&iter)) {

    // Get the cell boundaries and cell center.
    double cell_b[4] = {0.0}; // Cell boundaries, x lo and up, y lo and up;
    double xc[2] = {0.0}; // Cell center.
    ts_grid_cell_boundaries(&up->ts_grid, iter.idx, cell_b);
    gkyl_rect_grid_cell_center(&up->ts_grid, iter.idx, xc);

    int shear_idx[] = {iter.idx[up->shear_dir_in_ts_grid]};
    int shift_idx[] = {iter.idx[up->shift_dir_in_ts_grid]};
    long shift_loc = gkyl_range_idx(&up->shear_r, shear_idx);
    double *shift_c = (double *) gkyl_array_fetch(up->shift_dg, shift_loc);

    int num_do_curr = 0;

    for (int dC=0; dC<2; dC++) { // dC=0: x=const, dC=1: y=const (boundaries).
      for (int xS=0; xS<2; xS++) { // xS=0: lower, xS=1 upper (boundary).

        double test_pt[2] = {0.0}; // Test point to shift.
        for (int d=0; d<2; d++)
          test_pt[d] = cell_b[2*d]+delta[d];

        // Search first shifted point. Use pick_lower=false in find_cell unless
        // searching for points along a x=const line near the upper y-boundary. 
        bool pick_lower = dC==0 && xS==1;
        test_pt[dC] += xS*(up->ts_grid.dx[dC]-2.0*delta[dC]);

        // Shift the test point, find the cell that contains it, and if we
        // haven't included it yet, add it to our list of donors.
        ts_check_shifted_test_point(up, test_pt, xc, up->ts_grid.dx,
          shift_c, iter.idx, pick_lower, &num_do_curr, shift_dir_idx_do_buff);

        // Search for other shifted points along this line. 
        int step_dim = (dC+1) % 2;
        for (int sI=1; sI<num_test_pt[step_dim]; sI++) {
          test_pt[step_dim] += step_sz[step_dim];

          // Shift the test point, find the cell that contains it, and if we
          // haven't included it yet, add it to our list of donors.
          ts_check_shifted_test_point(up, test_pt, xc, up->ts_grid.dx,
            shift_c, iter.idx, pick_lower, &num_do_curr, shift_dir_idx_do_buff);
        }
      }
    }

    up->num_do[shear_idx[0]-up->shear_r.lower[0]] = num_do_curr;
  }

  // Copy the donor list to the persistent object and release the buffer.
  size_t buff_sz = gkyl_mem_buff_size(shift_dir_idx_do_buff);
  up->shift_dir_idx_do = (int *) gkyl_malloc(buff_sz);
  int *shift_dir_idx_do_buff_ptr = (int *) gkyl_mem_buff_data(shift_dir_idx_do_buff);
  memcpy(up->shift_dir_idx_do, shift_dir_idx_do_buff_ptr, buff_sz);
  gkyl_mem_buff_release(shift_dir_idx_do_buff);
}

struct gkyl_qr_res
ts_root_find(double (*func)(double,void*), void *ctx, const double *lims, int max_iter, double tol)
{
  // Use a Ridder's root finder to find the root of func in the interval
  // [lims[0],lims[1]] down to a tolerance 'tol'. Return the interval limit
  // if the function is smaller than the tolerance there. Return nil if the
  // function does not change sign in the interval (interval doesn't contain the root).
  double funcLo = func(lims[0], ctx), funcUp = func(lims[1], ctx);
//  if (fabs(funcLo) < tol)
//    return (struct gkyl_qr_res) {.res=lims[0], .status=0, .nevals=2};
//  else if (fabs(funcUp) < tol)                                     
//    return (struct gkyl_qr_res) {.res=lims[1], .status=0, .nevals=2};
//  else {
//    if (funcLo*funcUp < 0)
//      return gkyl_ridders(func, ctx, lims[0], lims[1], funcLo, funcUp, max_iter, tol);
//    else
//      return (struct gkyl_qr_res) {.status=1, .nevals=2};
//  }
  if (fabs(funcLo) > tol && fabs(funcUp) > tol) {
    if (funcLo*funcUp < 0)
      return gkyl_ridders(func, ctx, lims[0], lims[1], funcLo, funcUp, max_iter, tol);
    else
      return (struct gkyl_qr_res) {.status=1, .nevals=2};
  }
  else if (fabs(funcLo) < tol && fabs(funcUp) < tol)
    return (struct gkyl_qr_res) {.status=1, .nevals=2};
  else if (fabs(funcLo) < tol)
    return (struct gkyl_qr_res) {.res=lims[0], .status=0, .nevals=2};
  else if (fabs(funcUp) < tol)                                     
    return (struct gkyl_qr_res) {.res=lims[1], .status=0, .nevals=2};
  return (struct gkyl_qr_res) {.status=1, .nevals=2};
}

struct ts_shifted_coord_loss_func_ctx {
  double shiftCoordTar; // Target coordinate in shift_dir.
  double shiftCoordDo; // Donor coordinate in shift_dir.
  double shiftDirL; // Length of the domain in shift_dir.
  int periodicCopyIdx; // Used to search a periodic copy of the domain (signed).
  evalf_t shift_func; // Function defining the shift.
  void *shift_func_ctx; // Context for shift_func.
};

double ts_shifted_coord_loss_func(double shearCoord, void *ctx)
{
  // Loss function used to find the shear coord.
  struct ts_shifted_coord_loss_func_ctx *tsctx = ctx;

  double shift;
  tsctx->shift_func(0.0, (double[]){shearCoord}, &shift, tsctx->shift_func_ctx);

  return tsctx->shiftCoordTar - shift
    - (tsctx->shiftCoordDo - tsctx->periodicCopyIdx * tsctx->shiftDirL);
}

int static inline
ts_sign(double a)
{
  if (a < 0.0)
    return -1;
  else if (a > 0.0)
    return 1;
  else
    return 0;
}

double
ts_donor_target_offset(struct gkyl_twistshift_dg *up, const double *xc_do, const double *xc_tar) {
  // y-offset between the donor and the target cell (yDo-yTar), in the direction of the shift.
  //   xc_do:  cell center coordinates of donor cell.
  //   xc_tar: cell center coordinates of target cell.
  int shift_dir = up->shift_dir_in_ts_grid;
  double x_eval = xc_do[up->shear_dir];
  double shift;
  up->shift_func(0.0, (double[]){x_eval}, &shift, up->shift_func_ctx);

  double shift_dir_L = up->ts_grid.upper[up->shift_dir] - up->ts_grid.lower[up->shift_dir];

  // The idea here is that we keep shifting the donor cell center until it is in a
  // periodic copy of our domain which overlaps with the shifted target cell center.
  double xs_shifted_do = xc_do[shift_dir];
  double xs_shifted_tar = xc_tar[shift_dir] - shift;
  // Step toward the target.
  int shift_sign = xs_shifted_tar < xs_shifted_do? 1 : -1;
  bool keep_shifting = true;
  while (keep_shifting) {
    double xs_shifted_dolo = xs_shifted_do - shift_dir_L/2.0;
    double xs_shifted_doup = xs_shifted_do + shift_dir_L/2.0;
    if (xs_shifted_dolo <= xs_shifted_tar && xs_shifted_tar <= xs_shifted_doup) {
      keep_shifting = false;
      break;
    }
    else
      xs_shifted_do = xs_shifted_do - shift_sign*shift_dir_L;
  }
  return xc_tar[shift_dir] - xs_shifted_do;
}

struct gkyl_qr_res
ts_find_intersect(struct gkyl_twistshift_dg *up, double shiftCoordTar, double shiftCoordDo,
  const double *shearDirBounds, const double *shiftDirLimits, int nP_primary)
{
  // Given a y-coordinate of the target cell (yTar), and a y-coordinate
  // of the donor cell (yDo), find the x-coordinate where yTar-yShift(x)=yDo
  // for the periodic copy identified by nP_primary, i.e. the root of
  // yTar - yShift(x) - (yDo - nP_primary*Ly) = 0 in [shearDirBounds[0], shearDirBounds[1]].
  // Returns status=1 (not found) if no root exists for that specific nP_primary.
  //
  // nP_primary must be the same for all 4 inter_pts of a donor-target pair: all four
  // intersection conditions describe corners of a single physical overlap diamond and
  // therefore belong to the same periodic copy of the domain.
  double tol = 1.e-13;
  int max_iter = 100;

  double shiftDirL = shiftDirLimits[1] - shiftDirLimits[0];

  struct ts_shifted_coord_loss_func_ctx func_ctx = {
    .shiftCoordTar = shiftCoordTar,
    .shiftCoordDo = shiftCoordDo,
    .shiftDirL = shiftDirL,
    .periodicCopyIdx = nP_primary,
    .shift_func = up->shift_func,
    .shift_func_ctx = up->shift_func_ctx,
  };
  return ts_root_find(ts_shifted_coord_loss_func, &func_ctx, shearDirBounds, max_iter, tol);
}

struct ts_val_found {
  bool status; // =0 not found, =1 found.
  double value; // value found.
};

static inline void
ts_comp_to_phys(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d]+xc[d];
}

void
ts_nod2mod_proj_1d(struct gkyl_twistshift_dg *up, evalf_t func, void *func_ctx, const double *interval, double *out)
{
  // Project 'func' onto 1D DG basis in [interval[0], interval[1]]
  // evaluating at nodes and transforming to modal representation.
  //   func: 1D scalar function to be projected.
  //   interval: limits of the interval in which to project the function.
  //   out: DG field output.
  double dx[] = {0.0}, xc[] = {0.0}, xmu[] = {0.0};
  ts_interval_dx_and_xc(interval, dx, xc);

  for (int i=0; i<up->shift_b.num_basis; ++i) {
    ts_comp_to_phys(1, gkyl_eval_on_nodes_fetch_node(up->ev_on_nod1d, i),
      dx, xc, xmu);
    func(0.0, xmu, (double *)gkyl_array_fetch(up->func_nod1d,i), func_ctx);
  }
  gkyl_eval_on_nodes_nod2mod(up->ev_on_nod1d, up->func_nod1d, out); 
}

void
ts_integral_xlimdg(struct gkyl_twistshift_dg *up, double sFac, const double *xLimLo,
  const double *xLimUp, double yLimLo, double yLimUp, double dyDo, double yOff,
  const double *ySh, struct gkyl_mat *mat_do) {
  // Populate a matrix (mat_do) with a sub-cell integral that has variably x limits
  // represented by a DG polynomial, and a y-integral that goes from yLimLo to yLimUp.
  //   up: BC updater.
  //   sFac: +/-1 factor to add or subtract this subcell integral.
  //   xLimLo: DG representation of the lower x-limit.
  //   xLimUp: DG representation of the upper x-limit.
  //   yLimLo: lower y-limit.
  //   yLimUp: upper y-limit.
  //   dyDo: Cell length along y.
  //   yOff: Offset along y.
  //   ySh: DG representation of the y shift.
  //   mat_do: donor matrix.
  up->kernels->xlimdg(sFac, xLimLo, xLimUp, yLimLo, yLimUp, dyDo, yOff, ySh, mat_do);
}

void
ts_integral_ylimdg(struct gkyl_twistshift_dg *up, double sFac, double xLimLo, double xLimUp,
  const double *yLimLo, const double *yLimUp, double dyDo, double yOff,
  const double *ySh, struct gkyl_mat *mat_do) {
  // Populate a matrix (mat_do) with a sub-cell integral that has variable y limits
  // represented by a DG polynomial, and a x-integral that goes from xLimLo to xLimUp.
  //   up: BC updater.
  //   sFac: +/-1 factor to add or subtract this subcell integral.
  //   xLimLo: lower x-limit.
  //   xLimUp: upper x-limit.
  //   yLimLo: DG representation of the lower y-limit.
  //   yLimUp: DG representation of the upper y-limit.
  //   dyDo: Cell length along y.
  //   yOff: Offset along y.
  //   ySh: DG representation of the y shift.
  //   mat_do: donor matrix.
  up->kernels->ylimdg(sFac, xLimLo, xLimUp, yLimLo, yLimUp, dyDo, yOff, ySh, mat_do);
}

void
ts_integral_fullcelllimdg(struct gkyl_twistshift_dg *up, double dyDo, double yOff,
  const double *ySh, struct gkyl_mat *mat_do) {
  // Populate a matrix (mat_do) with the full-cell integral.
  //   up: BC updater.
  //   sFac: +/-1 factor to add or subtract this subcell integral.
  //   dyDo: Cell length along y.
  //   yOff: Offset along y.
  //   ySh: DG representation of the y shift.
  //   mat_do: donor matrix.
  up->kernels->fullcell(dyDo, yOff, ySh, mat_do);
}

static inline void
ts_one(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = 1.0;
}

static inline void
ts_minus_one(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = -1.0;
}

struct ts_shift_coord_shifted_log_ctx {
  double shift_coord_tar; // Target coordinate in shift_dir.
  int shift_sign_fac;
  const double *xc_do, *xc_tar; // Cell centers (donor and target).
  double *dx; // Cell lengths.
  bool pick_upper;
  int shear_dir, shift_dir; // Shear and shift directions.
  double shift_dir_bounds[2]; // Domain boundaries in shift_dir.
  evalf_t shift_func; // Function defining the shift.
  void *shift_func_ctx; // Context for shift_func.
};

void
ts_shift_coord_shifted_log(double t, const double *xn, double *fout, void *ctx)
{
  // Given a logical space x coordinate (xi) and a (physical) y-coordinate in the target cell,
  // compute the shifted y-coordinate in the logical space of the donor cell (eta \in [-1,1]).
  //   xi:    logical space x coordinate.
  //   yTar:  physical y-coordinate in target cell.
  //   pmSh:  factor multiplying the y-shift (+/- 1).
  //   xcDo:  cell center coordinates of donor cell.
  //   xcTar: cell center coordinates of target cell.
  //   dx:    cell lengths.
  //   pickUpper: boolean indicating if wrapping function should return upper/lower boundary.

  double xi = xn[0];

  struct ts_shift_coord_shifted_log_ctx *tsctx = ctx;
  double shift_coord_tar = tsctx->shift_coord_tar;
  int shift_sign_fac = tsctx->shift_sign_fac;
  const double *xc_do = tsctx->xc_do, *xc_tar = tsctx->xc_tar;
  double *dx = tsctx->dx;
  bool pick_upper = tsctx->pick_upper;
  int shear_dir = tsctx->shear_dir, shift_dir = tsctx->shift_dir;
  double *shift_dir_bounds = tsctx->shift_dir_bounds;

  double shear_coord_phys = xc_tar[shear_dir] + 0.5*dx[shear_dir]*xi;
  double shift;
  tsctx->shift_func(0.0, (double[]){shear_coord_phys}, &shift, tsctx->shift_func_ctx);

  double shift_coord_shifted = shift_coord_tar - shift_sign_fac * shift;
  shift_coord_shifted = ts_wrap_to_range(shift_coord_shifted, shift_dir_bounds[0], shift_dir_bounds[1], pick_upper);

  fout[0] = ts_p2l(shift_coord_shifted, xc_do[shift_dir], dx[shift_dir]);
}

void
ts_subcellint_sNi_sNii(struct gkyl_twistshift_dg *up, struct ts_val_found *inter_pts,
  const double *xc_do, const double *xc_tar, const double *cellb_do, const double *cellb_tar,
  bool is_upper_shift_dir_cell, const double *shift_c, struct gkyl_mat *mat_do)
{
  // Perform sNi or sNii subcell integrals.
  //   inter_pts: intersections y_{j_tar-/+1/2}-yShift and y_{j_do-/+1/2} (lower/upper y-boundaries of donor cell).
  //   xc_do: cell center of donor cell.
  //   xc_tar: cell center of target cell.
  //   cellb_do: boundaries of target cell.
  //   cellb_tar: boundaries of target cell.
  //   is_upper_shift_dir_cell: is the donor the upper cell in the shift dir?
  //   shift_c: DG coefficients of the shift.
  //   mat_do: current donor matrix.

  bool is_sNi = inter_pts[2].value < inter_pts[0].value;

  struct ts_shift_coord_shifted_log_ctx eta_lims_ctx = {
    .shift_sign_fac = 1,
    .xc_do = xc_do,
    .xc_tar = xc_tar,
    .dx = up->ts_grid.dx,
    .pick_upper = is_upper_shift_dir_cell,
    .shear_dir = up->shear_dir_in_ts_grid,
    .shift_dir = up->shift_dir_in_ts_grid,
    .shift_dir_bounds = {up->ts_grid.lower[up->shift_dir_in_ts_grid], 
                         up->ts_grid.upper[up->shift_dir_in_ts_grid]},
    .shift_func = up->shift_func,
    .shift_func_ctx = up->shift_func_ctx,
  };

  double xi_b[2]; // Limits of xi integral.
  evalf_t eta_lims[2]; // Table of functions definting the limits of eta integral.
  double etalo_xi[up->shift_b.num_basis], etaup_xi[up->shift_b.num_basis];

  // Offset between cell centers in direction of the shift.
  double xs_off = ts_donor_target_offset(up, xc_do, xc_tar);

  if (is_sNi) {
    // sNi
    // 1) Add the contribution of the left portion.
    xi_b[0] = ts_p2l(inter_pts[2].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
    xi_b[1] = ts_p2l(inter_pts[3].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);

    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_up(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_minus_one;
    eta_lims[1] = ts_shift_coord_shifted_log;

    ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
    ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);

    if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
      up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
        up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);

    // 2) Add the contribution of the right portion.
    xi_b[0] = ts_p2l(inter_pts[0].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
    xi_b[1] = ts_p2l(inter_pts[1].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);

    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_lo(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_shift_coord_shifted_log;
    eta_lims[1] = ts_one;

    ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
    ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);

    if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
      up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
        up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);
  }
  else {
    // sNii
    // 1) Add the contribution of the left portion.
    xi_b[0] = ts_p2l(inter_pts[1].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
    xi_b[1] = ts_p2l(inter_pts[0].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);

    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_lo(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_shift_coord_shifted_log;
    eta_lims[1] = ts_one;

    ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
    ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);

    if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
      up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
        up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);

    // 2) Add the contribution of the right portion.
    xi_b[0] = ts_p2l(inter_pts[3].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
    xi_b[1] = ts_p2l(inter_pts[2].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);

    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_up(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_minus_one;
    eta_lims[1] = ts_shift_coord_shifted_log;

    ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
    ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);

    if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
      up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
        up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);
  }
}

void
ts_subcellint_si_sii(struct gkyl_twistshift_dg *up, struct ts_val_found *inter_pts,
  const double *xc_do, const double *xc_tar, const double *cellb_do, const double *cellb_tar,
  bool is_upper_shift_dir_cell, const double *shift_c, struct gkyl_mat *mat_do)
{
  // Perform subcell integral si or sii, using fixed x-limits and variable y limits.
  //   inter_pts: intersections y_{j_tar-/+1/2}-yShift and y_{j_do-/+1/2} (lower/upper y-boundaries of donor cell).
  //   xc_do: cell center of donor cell.
  //   xc_tar: cell center of target cell.
  //   cellb_do: boundaries of target cell.
  //   cellb_tar: boundaries of target cell.
  //   is_upper_shift_dir_cell: is the donor the upper cell in the shift dir?
  //   shift_c: DG coefficients of the shift.
  //   mat_do: current donor matrix.

  double x_up = cellb_tar[cellb_up(up->shear_dir_in_ts_grid)];
  double x_lo = cellb_tar[cellb_lo(up->shear_dir_in_ts_grid)];
  double shift_lo, shift_up;

  up->shift_func(0.0, (double[]){x_lo}, &shift_lo, up->shift_func_ctx);
  up->shift_func(0.0, (double[]){x_up}, &shift_up, up->shift_func_ctx);

  bool is_si = -shift_lo < -shift_up;

  double xi_b[2];  // Limits of xi integral.
  if (is_si) {
    // si integral.
    xi_b[0] = -1.0;
    xi_b[1] = ts_p2l(inter_pts[1].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);;
  }
  else {
    // sii integral.
    xi_b[0] = ts_p2l(inter_pts[1].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);;
    xi_b[1] = 1.0;
  }

  evalf_t eta_lims[2]; // Table of functions definting the limits of eta integral.
  eta_lims[0] = ts_shift_coord_shifted_log;
  eta_lims[1] = ts_one;

  struct ts_shift_coord_shifted_log_ctx eta_lims_ctx = {
    .shift_coord_tar = cellb_tar[cellb_lo(up->shift_dir_in_ts_grid)],
    .shift_sign_fac = 1,
    .xc_do = xc_do,
    .xc_tar = xc_tar,
    .dx = up->ts_grid.dx,
    .pick_upper = is_upper_shift_dir_cell,
    .shear_dir = up->shear_dir_in_ts_grid,
    .shift_dir = up->shift_dir_in_ts_grid,
    .shift_dir_bounds = {up->ts_grid.lower[up->shift_dir_in_ts_grid], 
                         up->ts_grid.upper[up->shift_dir_in_ts_grid]},
    .shift_func = up->shift_func,
    .shift_func_ctx = up->shift_func_ctx,
  };

  double etalo_xi[up->shift_b.num_basis], etaup_xi[up->shift_b.num_basis];
  ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
  ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);
  // Offset between cell centers in direction of the shift.
  double xs_off = ts_donor_target_offset(up, xc_do, xc_tar);

  if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
    up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
      up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);
}

void
ts_subcellint_siii_siv(struct gkyl_twistshift_dg *up, struct ts_val_found *inter_pts,
  const double *xc_do, const double *xc_tar, const double *cellb_do, const double *cellb_tar,
  bool is_upper_shift_dir_cell, const double *shift_c, struct gkyl_mat *mat_do)
{
  // Perform subcell integral siii or siv, using fixed x-limits and variable y limits.
  //   inter_pts: intersections y_{j_tar-/+1/2}-yShift and y_{j_do-/+1/2} (lower/upper y-boundaries of donor cell).
  //   xc_do: cell center of donor cell.
  //   xc_tar: cell center of target cell.
  //   cellb_do: boundaries of target cell.
  //   cellb_tar: boundaries of target cell.
  //   is_upper_shift_dir_cell: is the donor the upper cell in the shift dir?
  //   shift_c: DG coefficients of the shift.
  //   mat_do: current donor matrix.
  
  double x_lo = cellb_tar[cellb_lo(up->shear_dir_in_ts_grid)];
  double x_up = cellb_tar[cellb_up(up->shear_dir_in_ts_grid)];
  double shift_lo, shift_up;

  up->shift_func(0.0, (double[]){x_lo}, &shift_lo, up->shift_func_ctx);
  up->shift_func(0.0, (double[]){x_up}, &shift_up, up->shift_func_ctx);

  bool is_siii = -shift_lo > -shift_up;

  double xi_b[2];  // Limits of xi integral.
  if (is_siii) {
    // siii integral.
    xi_b[0] = -1.0;
    xi_b[1] = ts_p2l(inter_pts[2].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);;
  }
  else {
    // siv integral.
    xi_b[0] = ts_p2l(inter_pts[2].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);;
    xi_b[1] = 1.0;
  }

  evalf_t eta_lims[2]; // Table of functions definting the limits of eta integral.
  eta_lims[0] = ts_minus_one;
  eta_lims[1] = ts_shift_coord_shifted_log;

  struct ts_shift_coord_shifted_log_ctx eta_lims_ctx = {
    .shift_coord_tar = cellb_tar[cellb_up(up->shift_dir_in_ts_grid)],
    .shift_sign_fac = 1,
    .xc_do = xc_do,
    .xc_tar = xc_tar,
    .dx = up->ts_grid.dx,
    .pick_upper = is_upper_shift_dir_cell,
    .shear_dir = up->shear_dir_in_ts_grid,
    .shift_dir = up->shift_dir_in_ts_grid,
    .shift_dir_bounds = {up->ts_grid.lower[up->shift_dir_in_ts_grid], 
                         up->ts_grid.upper[up->shift_dir_in_ts_grid]},
    .shift_func = up->shift_func,
    .shift_func_ctx = up->shift_func_ctx,
  };

  double etalo_xi[up->shift_b.num_basis], etaup_xi[up->shift_b.num_basis];
  ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
  ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);
  // Offset between cell centers in direction of the shift.
  double xs_off = ts_donor_target_offset(up, xc_do, xc_tar);

  if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
    up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
      up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);
}

void
ts_subcellint_sv_svi(struct gkyl_twistshift_dg *up, struct ts_val_found *inter_pts,
  const double *xc_do, const double *xc_tar, const double *cellb_do, const double *cellb_tar,
  bool is_upper_shift_dir_cell, const double *shift_c, struct gkyl_mat *mat_do)
{
  // Perform sv or svi subcell integrals.
  //   inter_pts: intersections y_{j_tar-/+1/2}-yShift and y_{j_do-/+1/2} (lower/upper y-boundaries of donor cell).
  //   xc_do: cell center of donor cell.
  //   xc_tar: cell center of target cell.
  //   cellb_do: boundaries of target cell.
  //   cellb_tar: boundaries of target cell.
  //   is_upper_shift_dir_cell: is the donor the upper cell in the shift dir?
  //   shift_c: DG coefficients of the shift.
  //   mat_do: current donor matrix.

  bool is_sv = inter_pts[3].value < inter_pts[1].value;

  struct ts_shift_coord_shifted_log_ctx eta_lims_ctx = {
    .shift_sign_fac = 1,
    .xc_do = xc_do,
    .xc_tar = xc_tar,
    .dx = up->ts_grid.dx,
    .pick_upper = is_upper_shift_dir_cell,
    .shear_dir = up->shear_dir_in_ts_grid,
    .shift_dir = up->shift_dir_in_ts_grid,
    .shift_dir_bounds = {up->ts_grid.lower[up->shift_dir_in_ts_grid], 
                         up->ts_grid.upper[up->shift_dir_in_ts_grid]},
    .shift_func = up->shift_func,
    .shift_func_ctx = up->shift_func_ctx,
  };

  double xi_b[2]; // Limits of xi integral.
  evalf_t eta_lims[2]; // Table of functions definting the limits of eta integral.
  double etalo_xi[up->shift_b.num_basis], etaup_xi[up->shift_b.num_basis];

  // Offset between cell centers in direction of the shift.
  double xs_off = ts_donor_target_offset(up, xc_do, xc_tar);

  if (is_sv) {
    // sv
    // 1) Add the contribution of the left portion.
    xi_b[0] = -1.0;
    xi_b[1] = ts_p2l(inter_pts[3].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);

    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_up(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_minus_one;
    eta_lims[1] = ts_shift_coord_shifted_log;

    ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
    ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);

    if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
      up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
        up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);

    // 2) Add the contribution of the right portion.
    xi_b[0] = ts_p2l(inter_pts[0].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
    xi_b[1] = ts_p2l(inter_pts[1].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);

    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_lo(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_shift_coord_shifted_log;
    eta_lims[1] = ts_one;

    ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
    ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);

    if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
      up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
        up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);
  }
  else {
    // svi
    // 1) Add the contribution of the left portion.
    xi_b[0] = ts_p2l(inter_pts[1].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
    xi_b[1] = ts_p2l(inter_pts[0].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);

    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_lo(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_shift_coord_shifted_log;
    eta_lims[1] = ts_one;

    ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
    ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);

    if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
      up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
        up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);

    // 2) Add the contribution of the right portion.
    xi_b[0] = ts_p2l(inter_pts[3].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
    xi_b[1] = 1.0;

    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_up(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_minus_one;
    eta_lims[1] = ts_shift_coord_shifted_log;

    ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
    ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);

    if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
      up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
        up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);
  }
}

void
ts_subcellint_svii_sviii(struct gkyl_twistshift_dg *up, struct ts_val_found *inter_pts,
  const double *xc_do, const double *xc_tar, const double *cellb_do, const double *cellb_tar,
  bool is_upper_shift_dir_cell, const double *shift_c, struct gkyl_mat *mat_do)
{
  // Perform svii or sviii subcell integrals.
  //   inter_pts: intersections y_{j_tar-/+1/2}-yShift and y_{j_do-/+1/2} (lower/upper y-boundaries of donor cell).
  //   xc_do: cell center of donor cell.
  //   xc_tar: cell center of target cell.
  //   cellb_do: boundaries of target cell.
  //   cellb_tar: boundaries of target cell.
  //   is_upper_shift_dir_cell: is the donor the upper cell in the shift dir?
  //   shift_c: DG coefficients of the shift.
  //   mat_do: current donor matrix.

  bool is_svii = inter_pts[0].value < inter_pts[2].value;

  struct ts_shift_coord_shifted_log_ctx eta_lims_ctx = {
    .shift_sign_fac = 1,
    .xc_do = xc_do,
    .xc_tar = xc_tar,
    .dx = up->ts_grid.dx,
    .pick_upper = is_upper_shift_dir_cell,
    .shear_dir = up->shear_dir_in_ts_grid,
    .shift_dir = up->shift_dir_in_ts_grid,
    .shift_dir_bounds = {up->ts_grid.lower[up->shift_dir_in_ts_grid], 
                         up->ts_grid.upper[up->shift_dir_in_ts_grid]},
    .shift_func = up->shift_func,
    .shift_func_ctx = up->shift_func_ctx,
  };

  double xi_b[2]; // Limits of xi integral.
  evalf_t eta_lims[2]; // Table of functions definting the limits of eta integral.
  double etalo_xi[up->shift_b.num_basis], etaup_xi[up->shift_b.num_basis];

  // Offset between cell centers in direction of the shift.
  double xs_off = ts_donor_target_offset(up, xc_do, xc_tar);

  if (is_svii) {
    // svii
    // 1) Add the contribution of the left portion.
    xi_b[0] = -1.0;
    xi_b[1] = ts_p2l(inter_pts[0].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);

    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_lo(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_shift_coord_shifted_log;
    eta_lims[1] = ts_one;

    ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
    ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);

    if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
      up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
        up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);

    // 2) Add the contribution of the right portion.
    xi_b[0] = ts_p2l(inter_pts[3].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
    xi_b[1] = ts_p2l(inter_pts[2].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);

    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_up(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_minus_one;
    eta_lims[1] = ts_shift_coord_shifted_log;

    ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
    ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);

    if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
      up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
        up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);
  }
  else {
    // sviii
    // 1) Add the contribution of the left portion.
    xi_b[0] = ts_p2l(inter_pts[2].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
    xi_b[1] = ts_p2l(inter_pts[3].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);

    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_up(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_minus_one;
    eta_lims[1] = ts_shift_coord_shifted_log;

    ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
    ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);

    if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
      up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
        up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);

    // 2) Add the contribution of the right portion.
    xi_b[0] = ts_p2l(inter_pts[0].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
    xi_b[1] = 1.0;

    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_lo(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_shift_coord_shifted_log;
    eta_lims[1] = ts_one;

    ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
    ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);

    if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
      up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
        up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);
  }
}

void
ts_subcellint_six_sx(struct gkyl_twistshift_dg *up, struct ts_val_found *inter_pts,
  const double *xc_do, const double *xc_tar, const double *cellb_do, const double *cellb_tar,
  bool is_upper_shift_dir_cell, const double *shift_c, struct gkyl_mat *mat_do)
{
  // Perform six or sx subcell integrals.
  //   inter_pts: intersections y_{j_tar-/+1/2}-yShift and y_{j_do-/+1/2} (lower/upper y-boundaries of donor cell).
  //   xc_do: cell center of donor cell.
  //   xc_tar: cell center of target cell.
  //   cellb_do: boundaries of target cell.
  //   cellb_tar: boundaries of target cell.
  //   is_upper_shift_dir_cell: is the donor the upper cell in the shift dir?
  //   shift_c: DG coefficients of the shift.
  //   mat_do: current donor matrix.

  bool is_six = inter_pts[0].value < inter_pts[1].value;

  // Limits of xi integral.
  double xi_b[2];
  if (is_six) {
    // six
    xi_b[0] = ts_p2l(inter_pts[0].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
    xi_b[1] = ts_p2l(inter_pts[1].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
  }
  else {
    // sx
    xi_b[0] = ts_p2l(inter_pts[1].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
    xi_b[1] = ts_p2l(inter_pts[0].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
  }

  struct ts_shift_coord_shifted_log_ctx eta_lims_ctx = {
    .shift_coord_tar = cellb_tar[cellb_lo(up->shift_dir_in_ts_grid)],
    .shift_sign_fac = 1,
    .xc_do = xc_do,
    .xc_tar = xc_tar,
    .dx = up->ts_grid.dx,
    .pick_upper = is_upper_shift_dir_cell,
    .shear_dir = up->shear_dir_in_ts_grid,
    .shift_dir = up->shift_dir_in_ts_grid,
    .shift_dir_bounds = {up->ts_grid.lower[up->shift_dir_in_ts_grid], 
                         up->ts_grid.upper[up->shift_dir_in_ts_grid]},
    .shift_func = up->shift_func,
    .shift_func_ctx = up->shift_func_ctx,
  };

  evalf_t eta_lims[2]; // Table of functions definting the limits of eta integral.
  eta_lims[0] = ts_shift_coord_shifted_log;
  eta_lims[1] = ts_one;

  double etalo_xi[up->shift_b.num_basis], etaup_xi[up->shift_b.num_basis];
  ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
  ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);
  // Offset between cell centers in direction of the shift.
  double xs_off = ts_donor_target_offset(up, xc_do, xc_tar);

  if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
    up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
      up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);
}

void
ts_subcellint_sxi_sxii(struct gkyl_twistshift_dg *up, struct ts_val_found *inter_pts,
  const double *xc_do, const double *xc_tar, const double *cellb_do, const double *cellb_tar,
  bool is_upper_shift_dir_cell, const double *shift_c, struct gkyl_mat *mat_do)
{
  // Perform sxi or sxii subcell integrals.
  //   inter_pts: intersections y_{j_tar-/+1/2}-yShift and y_{j_do-/+1/2} (lower/upper y-boundaries of donor cell).
  //   xc_do: cell center of donor cell.
  //   xc_tar: cell center of target cell.
  //   cellb_do: boundaries of target cell.
  //   cellb_tar: boundaries of target cell.
  //   is_upper_shift_dir_cell: is the donor the upper cell in the shift dir?
  //   shift_c: DG coefficients of the shift.
  //   mat_do: current donor matrix.

  bool is_sxi = inter_pts[3].value < inter_pts[2].value;

  // Limits of xi integral.
  double xi_b[2];
  if (is_sxi) {
    // six
    xi_b[0] = ts_p2l(inter_pts[3].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
    xi_b[1] = ts_p2l(inter_pts[2].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
  }
  else {
    // sx
    xi_b[0] = ts_p2l(inter_pts[2].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
    xi_b[1] = ts_p2l(inter_pts[3].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
  }

  struct ts_shift_coord_shifted_log_ctx eta_lims_ctx = {
    .shift_coord_tar = cellb_tar[cellb_up(up->shift_dir_in_ts_grid)],
    .shift_sign_fac = 1,
    .xc_do = xc_do,
    .xc_tar = xc_tar,
    .dx = up->ts_grid.dx,
    .pick_upper = is_upper_shift_dir_cell,
    .shear_dir = up->shear_dir_in_ts_grid,
    .shift_dir = up->shift_dir_in_ts_grid,
    .shift_dir_bounds = {up->ts_grid.lower[up->shift_dir_in_ts_grid], 
                         up->ts_grid.upper[up->shift_dir_in_ts_grid]},
    .shift_func = up->shift_func,
    .shift_func_ctx = up->shift_func_ctx,
  };

  evalf_t eta_lims[2]; // Table of functions definting the limits of eta integral.
  eta_lims[0] = ts_minus_one;
  eta_lims[1] = ts_shift_coord_shifted_log;

  double etalo_xi[up->shift_b.num_basis], etaup_xi[up->shift_b.num_basis];
  ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
  ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);
  // Offset between cell centers in direction of the shift.
  double xs_off = ts_donor_target_offset(up, xc_do, xc_tar);

  if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
    up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
      up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);
}

void
ts_subcellint_sxiii_sxiv(struct gkyl_twistshift_dg *up, struct ts_val_found *inter_pts,
  const double *xc_do, const double *xc_tar, const double *cellb_do, const double *cellb_tar,
  bool is_upper_shift_dir_cell, const double *shift_c, struct gkyl_mat *mat_do)
{
  // Perform sxiii or sxiv subcell integrals.
  //   inter_pts: intersections y_{j_tar-/+1/2}-yShift and y_{j_do-/+1/2} (lower/upper y-boundaries of donor cell).
  //   xc_do: cell center of donor cell.
  //   xc_tar: cell center of target cell.
  //   cellb_do: boundaries of target cell.
  //   cellb_tar: boundaries of target cell.
  //   is_upper_shift_dir_cell: is the donor the upper cell in the shift dir?
  //   shift_c: DG coefficients of the shift.
  //   mat_do: current donor matrix.

  double x_lo = cellb_tar[cellb_lo(up->shear_dir_in_ts_grid)];
  double x_up = cellb_tar[cellb_up(up->shear_dir_in_ts_grid)];
  double shift_lo, shift_up;

  up->shift_func(0.0, (double[]){x_lo}, &shift_lo, up->shift_func_ctx);
  up->shift_func(0.0, (double[]){x_up}, &shift_up, up->shift_func_ctx);

  bool is_sxiii = -shift_lo < -shift_up;

  struct ts_shift_coord_shifted_log_ctx eta_lims_ctx = {
    .shift_sign_fac = 1,
    .xc_do = xc_do,
    .xc_tar = xc_tar,
    .dx = up->ts_grid.dx,
    .pick_upper = is_upper_shift_dir_cell,
    .shear_dir = up->shear_dir_in_ts_grid,
    .shift_dir = up->shift_dir_in_ts_grid,
    .shift_dir_bounds = {up->ts_grid.lower[up->shift_dir_in_ts_grid], 
                         up->ts_grid.upper[up->shift_dir_in_ts_grid]},
    .shift_func = up->shift_func,
    .shift_func_ctx = up->shift_func_ctx,
  };

  double xi_b[2]; // Limits of xi integral.
  evalf_t eta_lims[2]; // Table of functions definting the limits of eta integral.
  double etalo_xi[up->shift_b.num_basis], etaup_xi[up->shift_b.num_basis];

  // Offset between cell centers in direction of the shift.
  double xs_off = ts_donor_target_offset(up, xc_do, xc_tar);

  // 1) Add the contribution of the left portion.
  xi_b[0] = -1.0;
  xi_b[1] = ts_p2l(inter_pts[3].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);

  if (is_sxiii) {
    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_up(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_minus_one;
    eta_lims[1] = ts_shift_coord_shifted_log;
  }
  else {
    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_lo(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_shift_coord_shifted_log;
    eta_lims[1] = ts_one;
  }

  ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
  ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);

  if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
    up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
      up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);

  // 2) Add the contribution of the right portion.
  xi_b[0] = ts_p2l(inter_pts[0].value, xc_do[up->shear_dir_in_ts_grid], up->ts_grid.dx[up->shear_dir_in_ts_grid]);
  xi_b[1] = 1.0;

  if (is_sxiii) {
    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_lo(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_shift_coord_shifted_log;
    eta_lims[1] = ts_one;
  }
  else {
    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_up(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_minus_one;
    eta_lims[1] = ts_shift_coord_shifted_log;
  }

  ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
  ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);

  if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
    up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
      up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);
}

void
ts_subcellint_sxv_sxvi(struct gkyl_twistshift_dg *up, struct ts_val_found *inter_pts,
  const double *xc_do, const double *xc_tar, const double *cellb_do, const double *cellb_tar,
  bool is_upper_shift_dir_cell, const double *shift_c, struct gkyl_mat *mat_do)
{
  // Perform subcell integral sxv or sxvi, using fixed x-limits and variable y limits.
  //   inter_pts: intersections y_{j_tar-/+1/2}-yShift and y_{j_do-/+1/2} (lower/upper y-boundaries of donor cell).
  //   xc_do: cell center of donor cell.
  //   xc_tar: cell center of target cell.
  //   cellb_do: boundaries of target cell.
  //   cellb_tar: boundaries of target cell.
  //   is_upper_shift_dir_cell: is the donor the upper cell in the shift dir?
  //   shift_c: DG coefficients of the shift.
  //   mat_do: current donor matrix.

  double xi_b[] = {-1.0, 1.0};  // Limits of xi integral.

  double shift_dir_bounds[] = {up->ts_grid.lower[up->shift_dir_in_ts_grid],
                               up->ts_grid.upper[up->shift_dir_in_ts_grid]};

  double x_eval = xc_do[up->shear_dir_in_ts_grid];
  double shift;

  up->shift_func(0.0, (double[]){x_eval}, &shift, up->shift_func_ctx);

  double shifted_coord = cellb_tar[cellb_lo(up->shift_dir_in_ts_grid)] - shift;
  shifted_coord = ts_wrap_to_range(shifted_coord, shift_dir_bounds[0], shift_dir_bounds[1],
    is_upper_shift_dir_cell);

  struct ts_shift_coord_shifted_log_ctx eta_lims_ctx = {
    .shift_sign_fac = 1,
    .xc_do = xc_do,
    .xc_tar = xc_tar,
    .dx = up->ts_grid.dx,
    .pick_upper = is_upper_shift_dir_cell,
    .shear_dir = up->shear_dir_in_ts_grid,
    .shift_dir = up->shift_dir_in_ts_grid,
    .shift_dir_bounds = {up->ts_grid.lower[up->shift_dir_in_ts_grid], 
                         up->ts_grid.upper[up->shift_dir_in_ts_grid]},
    .shift_func = up->shift_func,
    .shift_func_ctx = up->shift_func_ctx,
  };

  evalf_t eta_lims[2]; // Table of functions definting the limits of eta integral.
  if ( cellb_do[cellb_lo(up->shift_dir_in_ts_grid)] <= shifted_coord &&
       shifted_coord <= cellb_do[cellb_up(up->shift_dir_in_ts_grid)] ) {
    // sxv integral.
    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_lo(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_shift_coord_shifted_log;
    eta_lims[1] = ts_one;
  }
  else {
    // sxvi integral.
    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_up(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_minus_one;
    eta_lims[1] = ts_shift_coord_shifted_log;
  }

  double etalo_xi[up->shift_b.num_basis], etaup_xi[up->shift_b.num_basis];
  ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
  ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);
  // Offset between cell centers in direction of the shift.
  double xs_off = ts_donor_target_offset(up, xc_do, xc_tar);

  if (fabs(xi_b[1] - xi_b[0]) > tol_xi)
    up->kernels->ylimdg(1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
      up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);
}

struct gkyl_nmat *
ts_calc_mats(struct gkyl_twistshift_dg *up)
{

  // Allocate matrices containing the discrete subcell integrals.
  int num_do_tot = 0;
  for (int i=0; i<up->shear_r.volume; i++)
    num_do_tot += up->num_do[i];

  struct gkyl_nmat *matsdo = gkyl_nmat_new(num_do_tot, up->basis.num_basis, up->basis.num_basis);
  for (int n=0; n<matsdo->num; ++n) {
    struct gkyl_mat mat = gkyl_nmat_get(matsdo, n);
    for (int j=0; j<matsdo->nc; ++j)
      for (int i=0; i<matsdo->nr; ++i)
        gkyl_mat_set(&mat, i, j, 0.0);
  }

  // y-index of the reference target used to precalc matrices. For positive(negative)
  // yShift idx=1(last) might be better, but ideally it shouldn't matter.
  int shift_dir_idx_tar = 1;

  double shift_dir_lims[] = {up->ts_grid.lower[up->shift_dir_in_ts_grid],
                             up->ts_grid.upper[up->shift_dir_in_ts_grid]};

  // Create an eval_on_nodes updater to use its nodes and functions (but not
  // the whole advance method).
  up->ev_on_nod1d = gkyl_eval_on_nodes_new(&up->shear_grid, &up->shift_b, 1, ts_one, NULL);
  // Create an array to store evaluations of a function at 1D nodes.
  up->func_nod1d = gkyl_array_new(GKYL_DOUBLE, 1, up->shift_b.num_basis);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->shear_r);
  while (gkyl_range_iter_next(&iter)) {

    // Get the cell boundaries and cell center.
    int idx_tar[2]; // Target index.
    idx_tar[up->shift_dir_in_ts_grid] = shift_dir_idx_tar;
    idx_tar[up->shear_dir_in_ts_grid] = iter.idx[0];
    double cellb_tar[4] = {0.0}; // Cell boundaries, x lo and up, y lo and up;
    double xc_tar[2] = {0.0}; // Cell center.
    ts_grid_cell_boundaries(&up->ts_grid, idx_tar, cellb_tar);
    gkyl_rect_grid_cell_center(&up->ts_grid, idx_tar, xc_tar);

    long shift_loc = gkyl_range_idx(&up->shear_r, iter.idx);
    double *shift_c = (double *) gkyl_array_fetch(up->shift_dg, shift_loc);

    long linidx_do = ts_shift_dir_idx_do_linidx(up->num_do, iter.idx[0], shift_dir_idx_tar,
      up->ts_grid.cells[up->shift_dir_in_ts_grid], up->shear_r.lower[0]);
    int *shift_dir_idx_do_ptr = &up->shift_dir_idx_do[linidx_do];

    long linidx_mats_do = 0;
    for (int i=0; i<iter.idx[0]-up->shear_r.lower[0]; i++)
      linidx_mats_do += up->num_do[i];

    // Check that the shift variation within this x-cell < Ly.
    // The algorithm assumes at most one x-intersection per (y-boundary, y-boundary) pair,
    // which breaks when |S(x_up) - S(x_lo)| >= Ly.
    double x_lo = cellb_tar[cellb_lo(up->shear_dir_in_ts_grid)];
    double x_up = cellb_tar[cellb_up(up->shear_dir_in_ts_grid)];
    double Ly = shift_dir_lims[1] - shift_dir_lims[0];
    double S_lo, S_up, S_c;
    up->shift_func(0.0, (double[]){x_lo}, &S_lo, up->shift_func_ctx);
    up->shift_func(0.0, (double[]){x_up}, &S_up, up->shift_func_ctx);
    up->shift_func(0.0, (double[]){xc_tar[up->shear_dir_in_ts_grid]}, &S_c, up->shift_func_ctx);

    if (fabs(S_up - S_lo) >= Ly) {
      fprintf(stderr, "bc_twistshift: shift variation |S(x_up)-S(x_lo)| = %g across a single x-cell"
        " exceeds Ly = %g (cell ix=%d). Increase Nx, reduce the shear, or increase Ly.\n",
        fabs(S_up - S_lo), Ly, iter.idx[0]);
      assert(false);
    }

    for (int iC=0; iC<up->num_do[iter.idx[0]-up->shear_r.lower[0]]; iC++){
      int idx_do[2]; // Target index.
      idx_do[up->shift_dir_in_ts_grid] = shift_dir_idx_do_ptr[iC];
      idx_do[up->shear_dir_in_ts_grid] = iter.idx[0];

      double cellb_do[4] = {0.0}; // Cell boundaries, x lo and up, y lo and up;
      double xc_do[2] = {0.0}; // Cell center.
      ts_grid_cell_boundaries(&up->ts_grid, idx_do, cellb_do);
      gkyl_rect_grid_cell_center(&up->ts_grid, idx_do, xc_do);

      // Get the matrix we are presently assigning.
      struct gkyl_mat mat_do = gkyl_nmat_get(matsdo, linidx_mats_do+iC);

      // Periodic copy in which to find the target for this donor-target pair: the integer nP such that
      // S_c \approx y_tar_c - y_do_c + nP*Ly. All 4 inter_pts must use this same nP so that only
      // roots from the physical intersection are accepted (not roots from other periodic copies).
      int nP_primary = (int)round((S_c - (xc_tar[up->shift_dir_in_ts_grid] - xc_do[up->shift_dir_in_ts_grid])) / Ly);

      // Find the points where y_{j_tar-/+1/2}-yShift intersect the y=y_{j_do-/+1/2} lines.
      // Also record the number and indices of points found/not found.
      struct ts_val_found inter_pts[4] = {};
      int num_inter_pts_found = 0, num_inter_pts_not_found = 4;
      int inter_pts_found_idxs[4], inter_pts_not_found_idxs[4];
      for (int i=0; i<2; i++) { // Loop over j_tar-/+1/2
        for (int j=0; j<2; j++) { // Loop over j_do-/+1/2
          double shift_dir_coord_tar = cellb_tar[2*up->shift_dir_in_ts_grid+i];
          double shift_dir_coord_do = cellb_do[2*up->shift_dir_in_ts_grid+j];
          struct gkyl_qr_res inter_res = ts_find_intersect(up, shift_dir_coord_tar, shift_dir_coord_do,
            (double[]) {cellb_tar[cellb_lo(up->shear_dir_in_ts_grid)],cellb_tar[cellb_up(up->shear_dir_in_ts_grid)]},
            shift_dir_lims, nP_primary);
          int ip_linc = i*2+j;
          inter_pts[ip_linc].status = inter_res.status == 0;
          if (inter_res.status == 0) {
            inter_pts[ip_linc].value = inter_res.res;
            inter_pts_found_idxs[num_inter_pts_found] = ip_linc;
            num_inter_pts_found++;
          }
          else {
            inter_pts_not_found_idxs[num_inter_pts_not_found] = ip_linc;
            num_inter_pts_not_found--;
          }
        }
      }

      bool is_upper_shift_dir_cell = idx_do[up->shift_dir_in_ts_grid] == up->ts_grid.cells[up->shift_dir_in_ts_grid];

      if (num_inter_pts_found == 4) {
        // sN: all intersections are found at this cell.
        ts_subcellint_sNi_sNii(up, inter_pts, xc_do, xc_tar, cellb_do, cellb_tar, is_upper_shift_dir_cell, shift_c, &mat_do);
      }
      else if (num_inter_pts_found == 1) {
        if (inter_pts[1].status) {
          // si:   y_{j_tar-1/2}-yShift intersects x_{i-1/2}.
          // sii:  y_{j_tar-1/2}-yShift intersects x_{i+1/2}.
          ts_subcellint_si_sii(up, inter_pts, xc_do, xc_tar, cellb_do, cellb_tar, is_upper_shift_dir_cell, shift_c, &mat_do);
        }
        else {
          // siii: y_{j_tar+1/2}-yShift intersects x_{i-1/2}.
          // siv:  y_{j_tar+1/2}-yShift intersects x_{i+1/2}.
          ts_subcellint_siii_siv(up, inter_pts, xc_do, xc_tar, cellb_do, cellb_tar, is_upper_shift_dir_cell, shift_c, &mat_do);
        }
      }
      else if (num_inter_pts_found == 3) {
        if (!inter_pts[2].status) {
          // sv:    y_{j_tar+1/2}-yShift doesn't intersect y_{j_do-1/2} & intersects x_{i-1/2}.
          // svi:   y_{j_tar+1/2}-yShift doesn't intersect y_{j_do-1/2} & intersects x_{i+1/2}.
          ts_subcellint_sv_svi(up, inter_pts, xc_do, xc_tar, cellb_do, cellb_tar, is_upper_shift_dir_cell, shift_c, &mat_do);
        }
        else {
          // svii:  y_{j_tar-1/2}-yShift doesn't intersect y_{j_do+1/2} & intersects x_{i-1/2}.
          // sviii: y_{j_tar-1/2}-yShift doesn't intersect y_{j_do+1/2} & intersects x_{i+1/2}.
          ts_subcellint_svii_sviii(up, inter_pts, xc_do, xc_tar, cellb_do, cellb_tar, is_upper_shift_dir_cell, shift_c, &mat_do);
        }
      }
      else if (num_inter_pts_found == 2) {
        if (inter_pts[0].status && inter_pts[1].status) {
          // six:   y_{j_tar-1/2}-yShift crosses y_{j_do-/+1/2} (increasing yShift).
          // sx:    y_{j_tar-1/2}-yShift crosses y_{j_do-/+1/2} (decreasing yShift).
          ts_subcellint_six_sx(up, inter_pts, xc_do, xc_tar, cellb_do, cellb_tar, is_upper_shift_dir_cell, shift_c, &mat_do);
        }
        else if (inter_pts[2].status && inter_pts[3].status) {
          // sxi:   y_{j_tar+1/2}-yShift crosses y_{j_do-/+1/2} (decreasing yShift).
          // sxii:  y_{j_tar+1/2}-yShift crosses y_{j_do-/+1/2} (increasing yShift).
          ts_subcellint_sxi_sxii(up, inter_pts, xc_do, xc_tar, cellb_do, cellb_tar, is_upper_shift_dir_cell, shift_c, &mat_do);
        }
        else {
          // sxiii: y_{j_tar-1/2}-yShift crosses y_{j_do-1/2} & y_{j_tar+1/2}-yShift crosses y_{j_do+1/2} (increasing yShift).
          // sxiv:  y_{j_tar-1/2}-yShift crosses y_{j_do-1/2} & y_{j_tar+1/2}-yShift crosses y_{j_do+1/2} (decreasing yShift).
          ts_subcellint_sxiii_sxiv(up, inter_pts, xc_do, xc_tar, cellb_do, cellb_tar, is_upper_shift_dir_cell, shift_c, &mat_do);
        }
      }
      else if (num_inter_pts_found == 0) {
        // sxv:  y_{j_tar-1/2}-yShift crosses x_{i-/+1/2}.
        // sxvi: y_{j_tar+1/2}-yShift crosses x_{i-/+1/2}.
        ts_subcellint_sxv_sxvi(up, inter_pts, xc_do, xc_tar, cellb_do, cellb_tar, is_upper_shift_dir_cell, shift_c, &mat_do);
      }
      else {
        // An error occurred. This shouldn't happen.
        assert(false);
      }
    }

  }

  gkyl_array_release(up->func_nod1d);
  gkyl_eval_on_nodes_release(up->ev_on_nod1d);

  struct gkyl_nmat *matsdo_out = up->use_gpu? gkyl_nmat_cu_dev_new(matsdo->num, matsdo->nr, matsdo->nc)
                                            : gkyl_nmat_acquire(matsdo); 
  gkyl_nmat_copy(matsdo_out, matsdo);
  gkyl_nmat_release(matsdo);

  return matsdo_out;
}

long *
ts_calc_num_numcol_fidx_do(struct gkyl_twistshift_dg *up)
{
  // Calculate the linear indices into the donor distribution function gkyl_array
  // for each num-numcol plane (in the num-numcol-num_basis) space.

  long *num_numcol_fidx_do_ho = (long*) gkyl_malloc(up->fmat->num * up->fmat->nc * sizeof(long));

  // Location in the direction of the BC from which to take donor
  // distributions. We assume that the user filled the ghost cell with the skin
  // on the other side (i.e. applied periodicity first).
  int bc_dir_loc_do = up->edge == GKYL_LOWER_EDGE? up->local_bcdir_ext_r.lower[up->bc_dir]
                                                 : up->local_bcdir_ext_r.upper[up->bc_dir];

  // Range over directions other than shear and bc dirs.
  struct gkyl_range shearbc_perp_r;
  int remove[GKYL_MAX_DIM] = {0}, loc_in_dir[GKYL_MAX_DIM] = {0};;
  remove[up->shear_dir] = remove[up->bc_dir] = 1;
  loc_in_dir[up->shear_dir] = up->local_bcdir_ext_r.lower[up->shear_dir];
  loc_in_dir[up->bc_dir] = bc_dir_loc_do;
  gkyl_range_deflate(&shearbc_perp_r, &up->local_bcdir_ext_r, remove, loc_in_dir);

  int shift_dir_in_shearbc_perp_r;
  if (up->shift_dir < up->shear_dir && up->shift_dir < up->bc_dir)
    shift_dir_in_shearbc_perp_r = up->shift_dir;
  else if (up->shift_dir > up->shear_dir && up->shift_dir > up->bc_dir)
    shift_dir_in_shearbc_perp_r = up->shift_dir-2;
  else
    shift_dir_in_shearbc_perp_r = up->shift_dir-1;

  int prev_shift_dir_idx = 0;
  int donor_count = 0;
  int do_idx[up->local_bcdir_ext_r.ndim];

  // Loop over directions perpendicular to shear and BC dirs.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &shearbc_perp_r);
  while (gkyl_range_iter_next(&iter)) {
    int ic = 0;
    for (int d=0; d<up->local_bcdir_ext_r.ndim; d++) {
      if (d != up->bc_dir && d != up->shear_dir) {
        do_idx[d] = iter.idx[ic];
        ic++;
      }
    }
    do_idx[up->bc_dir] = bc_dir_loc_do;

    struct gkyl_range_iter shear_dir_iter;
    gkyl_range_iter_init(&shear_dir_iter, &up->shear_r);
    while (gkyl_range_iter_next(&shear_dir_iter)) {

      int shear_dir_idx = shear_dir_iter.idx[0];

      long linidx_do = ts_shift_dir_idx_do_linidx(up->num_do, shear_dir_idx,
        iter.idx[shift_dir_in_shearbc_perp_r], up->ts_grid.cells[up->shift_dir_in_ts_grid], up->shear_r.lower[0]);

      for (int i = 0; i < up->num_do[shear_dir_idx-up->shear_r.lower[0]]; i++) {
        do_idx[up->shear_dir] = shear_dir_idx;
        do_idx[up->shift_dir] = up->shift_dir_idx_do[linidx_do+i];

        long loc = gkyl_range_idx(&up->local_bcdir_ext_r, do_idx);
        num_numcol_fidx_do_ho[donor_count] = loc;

        donor_count += 1;
      }
    }
    prev_shift_dir_idx = iter.idx[shift_dir_in_shearbc_perp_r];
  }

  long *num_numcol_fidx_do;
  if (!up->use_gpu) {
    num_numcol_fidx_do = (long*) gkyl_malloc(up->fmat->num * up->fmat->nc * sizeof(long));
    memcpy(num_numcol_fidx_do, num_numcol_fidx_do_ho, up->fmat->num * up->fmat->nc * sizeof(long));
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    num_numcol_fidx_do = (long*) gkyl_cu_malloc(up->fmat->num * up->fmat->nc * sizeof(long));
    gkyl_cu_memcpy(num_numcol_fidx_do, num_numcol_fidx_do_ho, up->fmat->num * up->fmat->nc * sizeof(long), GKYL_CU_MEMCPY_H2D);
  }
#endif

  gkyl_free(num_numcol_fidx_do_ho);

  return num_numcol_fidx_do;
}

long *
ts_calc_num_numcol_fidx_tar(struct gkyl_twistshift_dg *up)
{

  long *num_numcol_fidx_tar_ho = (long*) gkyl_malloc(up->fmat->num * up->fmat->nc * sizeof(long));

  // Location in the direction of the BC in which to place the target
  // distributions. We assume that the local_bcdir_ext_r is a range extended in z (it
  // includes the z ghose cell).
  int bc_dir_loc_tar = up->edge == GKYL_LOWER_EDGE? up->local_bcdir_ext_r.lower[up->bc_dir]
                                                  : up->local_bcdir_ext_r.upper[up->bc_dir];

  // Range over directions other than shear and bc dirs.
  struct gkyl_range shearbc_perp_r;
  int remove[GKYL_MAX_DIM] = {0}, loc_in_dir[GKYL_MAX_DIM] = {0};;
  remove[up->shear_dir] = remove[up->bc_dir] = 1;
  loc_in_dir[up->shear_dir] = up->local_bcdir_ext_r.lower[up->shear_dir];
  loc_in_dir[up->bc_dir] = bc_dir_loc_tar;
  gkyl_range_deflate(&shearbc_perp_r, &up->local_bcdir_ext_r, remove, loc_in_dir);

  int tar_idx[up->local_bcdir_ext_r.ndim];
  int tar_count = 0;

  // Loop over directions perpendicular to shear and BC dirs.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &shearbc_perp_r);
  while (gkyl_range_iter_next(&iter)) {

    int ic = 0;
    for (int d=0; d<up->local_bcdir_ext_r.ndim; d++) {
      if (d != up->bc_dir && d != up->shear_dir) {
        tar_idx[d] = iter.idx[ic];
        ic++;
      }
    }
    tar_idx[up->bc_dir] = bc_dir_loc_tar;


    struct gkyl_range_iter shear_dir_iter;
    gkyl_range_iter_init(&shear_dir_iter, &up->shear_r);
    while (gkyl_range_iter_next(&shear_dir_iter)) {
      int shear_dir_idx = shear_dir_iter.idx[0];
      tar_idx[up->shear_dir] = shear_dir_idx;

      long loc = gkyl_range_idx(&up->local_bcdir_ext_r, tar_idx);
      num_numcol_fidx_tar_ho[tar_count] = loc;
      tar_count += 1;
    }
  }

  long *num_numcol_fidx_tar;
  if (!up->use_gpu) {
    num_numcol_fidx_tar = (long*) gkyl_malloc(up->fmat->num * up->fmat->nc * sizeof(long));
    memcpy(num_numcol_fidx_tar, num_numcol_fidx_tar_ho, up->fmat->num * up->fmat->nc * sizeof(long));
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    num_numcol_fidx_tar = (long*) gkyl_cu_malloc(up->fmat->num * up->fmat->nc * sizeof(long));
    gkyl_cu_memcpy(num_numcol_fidx_tar, num_numcol_fidx_tar_ho, up->fmat->num * up->fmat->nc * sizeof(long), GKYL_CU_MEMCPY_H2D);
  }
#endif

  gkyl_free(num_numcol_fidx_tar_ho);

  return num_numcol_fidx_tar;
}

void
gkyl_twistshift_dg_choose_kernels(struct gkyl_basis basis, int cdim, int shift_poly_order,
  struct gkyl_twistshift_dg_kernels *kers)
{
  int dim = basis.ndim;
  int vdim = dim - cdim;
  enum gkyl_basis_type basis_type = basis.b_type;
  int poly_order = basis.poly_order;
  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (shift_poly_order == 1) {
        kers->xlimdg   = vdim==0? ser_twistshift_xlimdg_list_0v_yShp1[cdim-2].kernels[poly_order]
                                : ser_twistshift_xlimdg_list_2v_yShp1[cdim-2].kernels[poly_order];
        kers->ylimdg   = vdim==0? ser_twistshift_ylimdg_list_0v_yShp1[cdim-2].kernels[poly_order]
                                : ser_twistshift_ylimdg_list_2v_yShp1[cdim-2].kernels[poly_order];
        kers->fullcell = vdim==0? ser_twistshift_fullcell_list_0v_yShp1[cdim-2].kernels[poly_order]
                                : ser_twistshift_fullcell_list_2v_yShp1[cdim-2].kernels[poly_order];
      }
      else if (shift_poly_order == 2) {
        assert(false); // MF 2025/09/20: removed 3x2v kernel because it's 8.5 MB.
        kers->xlimdg   = vdim==0? ser_twistshift_xlimdg_list_0v_yShp2[cdim-2].kernels[poly_order]
                                : ser_twistshift_xlimdg_list_2v_yShp2[cdim-2].kernels[poly_order];
        kers->ylimdg   = vdim==0? ser_twistshift_ylimdg_list_0v_yShp2[cdim-2].kernels[poly_order]
                                : ser_twistshift_ylimdg_list_2v_yShp2[cdim-2].kernels[poly_order];
        kers->fullcell = vdim==0? ser_twistshift_fullcell_list_0v_yShp2[cdim-2].kernels[poly_order]
                                : ser_twistshift_fullcell_list_2v_yShp2[cdim-2].kernels[poly_order];
      }
      return;
    default:
      assert(false);
      break;
  }
}

struct gkyl_twistshift_dg*
gkyl_twistshift_dg_inew(const struct gkyl_twistshift_dg_inp *inp)
{

  // Allocate space for new updater.
  struct gkyl_twistshift_dg *up = gkyl_malloc(sizeof(struct gkyl_twistshift_dg));

  up->bc_dir = inp->bc_dir;
  up->shift_dir = inp->shift_dir;
  up->shear_dir = inp->shear_dir;
  up->edge = inp->edge;
  up->basis = *inp->basis;
  up->grid = *inp->grid;
  up->use_gpu = inp->use_gpu;
  up->local_bcdir_ext_r = *inp->bcdir_ext_update_r;

  // Assume the poly order of the DG shift is the same as that of the field,
  // unless requested otherwise.
  up->shift_poly_order = inp->basis->poly_order;
  if (inp->shift_poly_order)
    up->shift_poly_order = inp->shift_poly_order;

  const int ndim = inp->bcdir_ext_update_r->ndim;
  // Check that it is being used for 3D or 5D. Likely only small changes are
  // needed to make it work in other dimensions.
  assert(ndim == 3 || ndim == 5); 

  double lo1d[1], up1d[1];  int cells1d[1];

  // Create 1D grid and range in the direction of the shear.
  gkyl_range_init(&up->shear_r, 1, (int[]) {up->local_bcdir_ext_r.lower[inp->shear_dir]},
                                   (int[]) {up->local_bcdir_ext_r.upper[inp->shear_dir]});
  lo1d[0] = inp->grid->lower[up->shear_dir];
  up1d[0] = inp->grid->upper[up->shear_dir];
  cells1d[0] = inp->grid->cells[up->shear_dir];
  gkyl_rect_grid_init(&up->shear_grid, 1, lo1d, up1d, cells1d);
  int idx[] = {up->shear_r.lower[0]};
  long linidx = gkyl_range_idx(&up->shear_r, idx);

  // Create 1D grid and range in the diretion of the shift.
  gkyl_range_init(&up->shift_r, 1, (int[]) {up->local_bcdir_ext_r.lower[inp->shift_dir]},
                                   (int[]) {up->local_bcdir_ext_r.upper[inp->shift_dir]});
  lo1d[0] = inp->grid->lower[up->shift_dir];
  up1d[0] = inp->grid->upper[up->shift_dir];
  cells1d[0] = inp->grid->cells[up->shift_dir];
  gkyl_rect_grid_init(&up->shift_grid, 1, lo1d, up1d, cells1d);

  // Create 2D grid (and range) the twist-shift takes place in.
  int dimlo, dimup;
  if (up->shift_dir < up->shear_dir) {
    dimlo = up->shift_dir;
    dimup = up->shear_dir;
    up->shift_dir_in_ts_grid = 0;
    up->shear_dir_in_ts_grid = 1;
  }
  else {
    dimlo = up->shear_dir;
    dimup = up->shift_dir;
    up->shift_dir_in_ts_grid = 1;
    up->shear_dir_in_ts_grid = 0;
  }
  gkyl_range_init(&up->ts_r, 2, (int[]) {up->local_bcdir_ext_r.lower[dimlo], up->local_bcdir_ext_r.lower[dimup]},
                                (int[]) {up->local_bcdir_ext_r.upper[dimlo], up->local_bcdir_ext_r.upper[dimup]});
  double lo2d[] = {inp->grid->lower[dimlo], inp->grid->lower[dimup]};
  double up2d[] = {inp->grid->upper[dimlo], inp->grid->upper[dimup]};
  int cells2d[] = {inp->grid->cells[dimlo], inp->grid->cells[dimup]};
  gkyl_rect_grid_init(&up->ts_grid, 2, lo2d, up2d, cells2d);

  // Project the shift onto the shift basis.
  gkyl_cart_modal_serendip(&up->shift_b, 1, up->shift_poly_order);
  if (inp->shift_func) {
    up->shift_dg = gkyl_array_new(GKYL_DOUBLE, up->shift_b.num_basis, up->shear_r.volume);
    gkyl_eval_on_nodes *evup = gkyl_eval_on_nodes_new(&up->shear_grid, &up->shift_b, 1,
      inp->shift_func, inp->shift_func_ctx);
    gkyl_eval_on_nodes_advance(evup, 0.0, &up->shear_r, up->shift_dg);
    gkyl_eval_on_nodes_release(evup);
  }
  else {
    // The shift must be discretized on the shear cells this updater indexes.
    assert(inp->shift_dg->size == up->shear_r.volume);
    assert(inp->shift_dg->ncomp == up->shift_b.num_basis);
    up->shift_dg = gkyl_array_acquire(inp->shift_dg);
  }

  // Function defining the shift (and its context).
  if (shift_func_op == 0) {
    up->shift_func = ts_shift_dg_eval;
    up->shift_dg_eval_ctx.shift_dg = up->shift_dg;
    up->shift_dg_eval_ctx.shift_b = &up->shift_b;
    up->shift_dg_eval_ctx.shear_grid = &up->shear_grid;
    up->shift_dg_eval_ctx.shear_r = &up->shear_r;
    up->shift_func_ctx = &up->shift_dg_eval_ctx;
  }
  else if (shift_func_op == 1) {
    up->shift_func     = inp->shift_func;
    up->shift_func_ctx = inp->shift_func_ctx;
  }
  else {
    fprintf(stderr, "Twist-shift function option not recognized. Exiting...\n");
    assert(false);
  }

  // Find the donor cells for each target. Store the number of donors for each
  // shear_dir idx (num_do) & the shift_dir idx of each donor (shift_dir_idx_do).
  // i.e. allocates and assigns num_do and up->shift_dir_idx_do.
  ts_find_donors(up);

  // Array of cummulative number of donors at given shear_dir cell.
  const int num_do_cum_sz = up->grid.cells[up->shear_dir]+1;
  int num_do_cum_ho[num_do_cum_sz];
  for (int i=0; i<num_do_cum_sz; i++)
    num_do_cum_ho[i] = 0;
  for (int i=up->shear_r.lower[0]; i<up->shear_r.upper[0]+1; i++)
    num_do_cum_ho[i] = num_do_cum_ho[i-1] + up->num_do[i-up->shear_r.lower[0]];

  if (!up->use_gpu) {
    up->num_do_cum = gkyl_malloc(num_do_cum_sz * sizeof(int));
    memcpy(up->num_do_cum, num_do_cum_ho, num_do_cum_sz * sizeof(int));
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    up->num_do_cum = gkyl_cu_malloc(num_do_cum_sz * sizeof(int));
    gkyl_cu_memcpy(up->num_do_cum, num_do_cum_ho, num_do_cum_sz * sizeof(int), GKYL_CU_MEMCPY_H2D);
  }
#endif

  // Choose the kernels that do the subcell and full cell integrals
  up->kernels = gkyl_malloc(sizeof(struct gkyl_twistshift_dg_kernels));
  gkyl_twistshift_dg_choose_kernels(*inp->basis, inp->cdim, up->shift_poly_order, up->kernels);

  // The BC is applied as a set of matrix-matrix multiplications
  //   f_i = sum_{q}^{N_do(i)} A_q,i B_q,i
  // where i indicates the shear_dir cell index, A_q is a
  // num_basis x num_basis matrix containing the discretization
  // of subcell integrals, B_q is a num_basis x (Ny * Nvpar * Nmu) matrix with
  // the DG coefficients of f common to a given A_q matrix, and thus where f_i
  // is a num_basis x (Ny*Nvpar*Nmu) matrix.
  //
  // Naming scheme:
  //   A_q: scimat (subscell integral matrices).
  //   B_q: fmat (distribution function matrices).
  //   A_q . B_q: mm_contr (contributions from mat-mat multiplication).

  // Calculate the entries in the matrices used to apply the BC.
  up->scimat = ts_calc_mats(up);

  // Number of colums in fmat.
  int fmat_num_col = 1;
  for (int d=0; d<ndim; d++) {
    if (d != up->bc_dir && d != up->shear_dir)
      fmat_num_col *= up->local_bcdir_ext_r.upper[d] - up->local_bcdir_ext_r.lower[d] + 1;
  }

  if (!up->use_gpu) {
    up->fmat = gkyl_nmat_new(up->scimat->num, up->scimat->nr, fmat_num_col);
    up->mm_contr = gkyl_nmat_new(up->scimat->num, up->scimat->nr, fmat_num_col);
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    up->fmat = gkyl_nmat_cu_dev_new(up->scimat->num, up->scimat->nr, fmat_num_col);
    up->mm_contr = gkyl_nmat_cu_dev_new(up->scimat->num, up->scimat->nr, fmat_num_col);
  }
#endif

  // Index translation from num-numcol plane index to linear index into the
  // donor distribution function gkyl_array.
  up->num_numcol_fidx_do = ts_calc_num_numcol_fidx_do(up);

  // Index translation from num-numcol plane index to linear index into the
  // tar distribution function gkyl_array.
  up->num_numcol_fidx_tar = ts_calc_num_numcol_fidx_tar(up);

  // Permutted ghost range, for indexing into the target field.
  // Order: Shift direction, redundant directions, shear direction.
  int lo4D[ndim-1], up4D[ndim-1];
  lo4D[0] = up->local_bcdir_ext_r.lower[up->shift_dir];
  up4D[0] = up->local_bcdir_ext_r.upper[up->shift_dir];
  int ic = 1;
  for (int d=0; d<ndim; d++) {
    if (d != up->bc_dir && d != up->shear_dir && d != up->shift_dir) {
      lo4D[ic] = up->local_bcdir_ext_r.lower[d];
      up4D[ic] = up->local_bcdir_ext_r.upper[d];
      ic++;
    }
  }
  lo4D[ndim-2] = up->local_bcdir_ext_r.lower[up->shear_dir];
  up4D[ndim-2] = up->local_bcdir_ext_r.upper[up->shear_dir];
  gkyl_range_init(&up->permutted_ghost_r, ndim-1, lo4D, up4D);

  // Create a ghost range, to clear it before adding contributions from TS BC.
  if (inp->edge == GKYL_LOWER_EDGE)
    gkyl_range_shorten_from_above(&up->ghost_r, &up->local_bcdir_ext_r, inp->bc_dir, inp->num_ghost[inp->bc_dir]);
  else
    gkyl_range_shorten_from_below(&up->ghost_r, &up->local_bcdir_ext_r, inp->bc_dir, inp->num_ghost[inp->bc_dir]);

  return up;
}

struct gkyl_twistshift_dg*
gkyl_twistshift_dg_new(int bc_dir, int shift_dir, int shear_dir,
  enum gkyl_edge_loc edge, int cdim, const struct gkyl_range *bcdir_ext_update_r, const int *num_ghost,
  const struct gkyl_basis *basis, const struct gkyl_rect_grid *grid, evalf_t shift_func, void *shift_func_ctx,
  struct gkyl_array *shift_dg, int shift_poly_order, bool use_gpu)
{
  struct gkyl_twistshift_dg_inp inp = {
    .bc_dir              = bc_dir            ,
    .shift_dir           = shift_dir         ,
    .shear_dir           = shear_dir         ,
    .edge                = edge              ,
    .cdim                = cdim              ,
    .bcdir_ext_update_r  = bcdir_ext_update_r,
    .num_ghost           = num_ghost         ,
    .basis               = basis             ,
    .grid                = grid              ,
    .shift_func          = shift_func        ,
    .shift_func_ctx      = shift_func_ctx    ,
    .shift_dg            = shift_dg          ,
    .use_gpu             = use_gpu           ,
    .shift_poly_order    = shift_poly_order  ,
  };
  return gkyl_twistshift_dg_inew(&inp);
}

void
gkyl_twistshift_dg_advance(struct gkyl_twistshift_dg *up, struct gkyl_array *fdo, struct gkyl_array *ftar)
{

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_twistshift_dg_advance_cu(up, fdo, ftar);
    return;
  }
#endif

  // Assign the distribution matrices.
  // This assumes that fdo->ncomp = fmat->nr.
  // Recall:
  //   fmat->num = sum_i^Nx num_do(i)
  //   fmat->nr = num_basis = ncomp
  //   fmat->nc = Ny*Nvpar*Nmu
  for (size_t i=0; i<up->fmat->num * up->fmat->nr * up->fmat->nc; i++) {

    long nc_idx = i / up->fmat->nr; // num-num_col index: current num-num_col plane.
    int row_idx = i % up->fmat->nr; // row index: current DG coeff.

    // This if-statement may only be needed in GPU kernel, not for CPUs.
    if ((nc_idx < up->fmat->num * up->fmat->nc) && (row_idx < fdo->ncomp)) {
      const double *fdo_c = (const double*) gkyl_array_cfetch(fdo, up->num_numcol_fidx_do[nc_idx]);
      struct gkyl_mat mcurr = gkyl_nmat_get(up->fmat, nc_idx % up->fmat->num);

      gkyl_mat_set(&mcurr, row_idx, nc_idx/up->fmat->num, fdo_c[row_idx]);
    }
  }

  // Perform the mat-mat multiplications.
  gkyl_nmat_mm(1.0, 0.0, GKYL_NO_TRANS, up->scimat, GKYL_NO_TRANS, up->fmat, up->mm_contr);

  // Clear the ghost range.
  gkyl_array_clear_range(ftar, 0.0, &up->ghost_r);

  // Perform reduction over num_do contributions from mat-mat mults (mm_contr).
  int num_cells_skin = (up->shear_r.upper[0]-up->shear_r.lower[0]+1) * up->fmat->nc;
  for (size_t i=0; i<ftar->ncomp * num_cells_skin; i++) {

    long linidx_tar = i / ftar->ncomp;
    int row_idx = i % ftar->ncomp;

    // This if-statement may only be needed in GPU kernel, not for CPUs.
    if ((linidx_tar < num_cells_skin) && (row_idx < ftar->ncomp)) {
      double *ftar_c = (double*) gkyl_array_fetch(ftar, up->num_numcol_fidx_tar[linidx_tar]);

      int idx[GKYL_MAX_DIM] = {1};
      gkyl_sub_range_inv_idx(&up->permutted_ghost_r, linidx_tar, idx);

      int ac[GKYL_MAX_DIM] = {1};
      for (int d=2; d<up->grid.ndim-1; d++)
        ac[d-2] = up->grid.cells[d+1];
      ac[up->permutted_ghost_r.ndim-2] = up->mm_contr->num;
  
      int start = 0;
      for (int d=0; d<up->permutted_ghost_r.ndim-1; d++)
        start = (start + (idx[d]-1)) * ac[d];

      int shear_idx = idx[up->permutted_ghost_r.ndim-1];

      int do_start = up->num_do_cum[shear_idx-1];
      int do_end   = up->num_do_cum[shear_idx-1+1];
      for (int j=do_start; j<do_end; j++) { // Only loop over num_do[i] elements.
        long linidx_mm_contr = start + j;
        struct gkyl_mat mat = gkyl_nmat_get(up->mm_contr, linidx_mm_contr % up->mm_contr->num);
        ftar_c[row_idx] += gkyl_mat_get(&mat, row_idx, linidx_mm_contr / up->mm_contr->num);
      }
    }
  }
}

struct gkyl_array*
gkyl_twistshift_dg_get_shift_objects(struct gkyl_twistshift_dg *up, struct gkyl_rect_grid *shear_grid,
  struct gkyl_range *shear_r, struct gkyl_basis *shift_b)
{
  *shear_grid = up->shear_grid;
  *shear_r    = up->shear_r   ;
  *shift_b    = up->shift_b   ;
  return gkyl_array_acquire(up->shift_dg);
};

void
gkyl_twistshift_dg_release(struct gkyl_twistshift_dg *up) {
  // Release memory associated with this updater.
  if (!up->use_gpu) {
    gkyl_free(up->num_do_cum);
    gkyl_free(up->num_numcol_fidx_do);
    gkyl_free(up->num_numcol_fidx_tar);
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->num_do_cum);
    gkyl_cu_free(up->num_numcol_fidx_do);
    gkyl_cu_free(up->num_numcol_fidx_tar);
  }
#endif

  gkyl_nmat_release(up->fmat);
  gkyl_nmat_release(up->mm_contr);

  gkyl_nmat_release(up->scimat);

  gkyl_free(up->kernels);

  gkyl_array_release(up->shift_dg);

  gkyl_free(up->num_do);
  gkyl_free(up->shift_dir_idx_do);

  gkyl_free(up);
}
