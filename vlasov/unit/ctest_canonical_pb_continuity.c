#include <acutest.h>

#include <gkyl_app.h>
#include <gkyl_app_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_array_rio.h>
#include <gkyl_canonical_pb_kernels.h> 
#include <gkyl_dg_calc_canonical_pb_vars.h>
#include <gkyl_dg_calc_canonical_pb_vars_priv.h>
#include <gkyl_eqn_type.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_lte_correct.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>
#include <gkyl_util.h>
#include <math.h>


void
info_h_ij_inv_1x1v(double t, const double* xn, double* fout, void* ctx)
{
  fout[0] = 1;
}

void
info_det_h_1x1v(double t, const double* xn, double* fout, void* ctx)
{
  fout[0] = 1;
}

void
info_hamil_1x1v(double t, const double* xn, double* fout, void* ctx)
{
  double x = xn[0], v = xn[1];
  fout[0] = 0.5*v*v;
}


void 
test_1x1v(int poly_order, enum gkyl_basis_type b_type)
{
  double pi = 3.14159265359;
  double lower[] = { 0.0, -5.0}, upper[] = {1.0, 5.0};
  int cells[] = {4, 4};
  int vdim = 1, cdim = 1;
  int pdim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1]}, velUpper[] = {upper[1]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    gkyl_cart_modal_hybrid(&basis, cdim, vdim);
  } 
  else if (b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  }
  else {
    gkyl_cart_modal_tensor(&basis, pdim, poly_order);;
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext; 
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0};
  struct gkyl_range local, local_ext; 
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Allocate arrays for specified metric inverse, hamiltonian and metric determinant
  struct gkyl_array *h_ij_inv, *det_h, *hamil;
  h_ij_inv = mkarr(false, confBasis.num_basis*vdim*(vdim+1)/2, local_ext.volume);
  det_h = mkarr(false, confBasis.num_basis, local_ext.volume);
  hamil = mkarr(false, basis.num_basis, local_ext.volume);

  // Evaluate specified inverse metric function and det. at nodes to insure continuity
  struct gkyl_eval_on_nodes* h_ij_inv_proj = gkyl_eval_on_nodes_new(&confGrid, &confBasis, vdim*(vdim+1)/2, info_h_ij_inv_1x1v, 0);
  struct gkyl_eval_on_nodes* det_h_proj = gkyl_eval_on_nodes_new(&confGrid, &confBasis, 1, info_det_h_1x1v, 0);
  struct gkyl_eval_on_nodes* hamil_proj = gkyl_eval_on_nodes_new(&grid, &basis, 1, info_hamil_1x1v, 0);
  gkyl_eval_on_nodes_advance(h_ij_inv_proj, 0.0, &confLocal, h_ij_inv);
  gkyl_eval_on_nodes_advance(det_h_proj, 0.0, &confLocal, det_h);
  gkyl_eval_on_nodes_advance(hamil_proj, 0.0, &local_ext, hamil);

  // Need to figure out size of alpha_surf and sgn_alpha_surf by finding size of surface basis set 
  struct gkyl_basis surf_basis, surf_quad_basis;
  if (basis.b_type == GKYL_BASIS_MODAL_HYBRID) {
    // NOTE: If we are hybrid, allocate more memory than we need to avoid surface basis are different
    // sizes in each direction.
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, 2);
    gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, 2);
  } 
  else {
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, poly_order);
    gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, poly_order);
  }

  // always vdim + cdim
  int alpha_surf_sz = (vdim + cdim)*surf_basis.num_basis; 
  int sgn_alpha_surf_sz = (vdim + cdim)*surf_quad_basis.num_basis; // sign(alpha) is store at quadrature points

  // allocate arrays to store fields: 
  // 1. alpha_surf (surface phase space velocity)
  // 2. sgn_alpha_surf (sign(alpha_surf) at quadrature points)
  // 3. const_sgn_alpha (boolean for if sign(alpha_surf) is a constant, either +1 or -1)
  struct gkyl_array *alpha_surf = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *sgn_alpha_surf = mkarr(false, sgn_alpha_surf_sz, local_ext.volume);
  struct gkyl_array *const_sgn_alpha = mk_int_arr(false, (vdim + cdim), local_ext.volume);

  // Pre-compute alpha_surf, sgn_alpha_surf, const_sgn_alpha, and cot_vec since they are time-independent
  struct gkyl_dg_calc_canonical_pb_vars *calc_vars = gkyl_dg_calc_canonical_pb_vars_new(&grid, 
    &confBasis, &basis, false);
  gkyl_dg_calc_canonical_pb_vars_alpha_surf(calc_vars, &confLocal, &local, &local_ext, hamil,
    alpha_surf, sgn_alpha_surf, const_sgn_alpha);

  // Check continuity (Directly via Kernels, Option 2)
  double w_edge_not_used[2] = { 0.0, 0.0 };
  double dxv[2] = {confGrid.dx[0], vel_grid.dx[0] };
  struct gkyl_array *alpha_surf_comp_L = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *alpha_surf_comp_R = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *sgn_alpha_surf_comp_not_used = mkarr(false, sgn_alpha_surf_sz, local_ext.volume);

  // Loop over cells, get expansions in cells
  struct gkyl_range_iter piter;
  gkyl_range_iter_init(&piter, &local_ext);
  while(gkyl_range_iter_next(&piter)){
    long pidx = gkyl_range_idx(&local_ext, piter.idx);
    const double *hamil_local = gkyl_array_cfetch(hamil,pidx);
    double *alpha_surf_comp_L_local = gkyl_array_fetch(alpha_surf_comp_L,pidx);
    double *alpha_surf_comp_R_local = gkyl_array_fetch(alpha_surf_comp_R,pidx);
    double *sgn_alpha_surf_comp_not_used_local = gkyl_array_fetch(sgn_alpha_surf_comp_not_used,pidx);

    if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_1x1v_ser_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_1x1v_ser_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfvx_1x1v_ser_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
    } 
    else if (b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_1x1v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_1x1v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfvx_1x1v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
    }
    else if ((b_type == GKYL_BASIS_MODAL_TENSOR) && (poly_order == 1)) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_1x1v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_1x1v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfvx_1x1v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
    } 
    else if ((b_type == GKYL_BASIS_MODAL_TENSOR) && (poly_order == 2)) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_1x1v_tensor_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_1x1v_tensor_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfvx_1x1v_tensor_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
    } 
    else {
      assert(false);
    }
  }

  // Compare suface wise element 
  // Loop over Number of surface coeff
  gkyl_range_iter_init(&piter, &local);
  int pidxl[GKYL_MAX_DIM];
  while(gkyl_range_iter_next(&piter)){
    long pidx = gkyl_range_idx(&local, piter.idx);
    const double *alpha_surf_comp_L_local = gkyl_array_cfetch(alpha_surf_comp_L, pidx);

    // Iterate in the lower direction (only conf space comparison - becuase there are no velocity space ghost cells for comparison)
    for (int dir = 0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(pdim, piter.idx, pidxl);
      pidxl[dir] = pidxl[dir]-1;
      long loc_l = gkyl_range_idx(&local, pidxl);
      const double *alpha_surf_comp_R_local = gkyl_array_cfetch(alpha_surf_comp_R, loc_l);

      // Iterate overthe number of basis on the surface
      for (int n = 0; n<surf_basis.num_basis; n++){
        TEST_CHECK(gkyl_compare_double(alpha_surf_comp_L_local[n + dir*surf_basis.num_basis], alpha_surf_comp_R_local[n + dir*surf_basis.num_basis], 1e-12));
      }
    }
  }

  // Release 
  gkyl_array_release(h_ij_inv);
  gkyl_array_release(det_h);
  gkyl_array_release(hamil);
  gkyl_array_release(alpha_surf);
  gkyl_array_release(sgn_alpha_surf);
  gkyl_array_release(const_sgn_alpha);
  gkyl_array_release(alpha_surf_comp_L);
  gkyl_array_release(alpha_surf_comp_R);
  gkyl_array_release(sgn_alpha_surf_comp_not_used);
  gkyl_eval_on_nodes_release(h_ij_inv_proj);
  gkyl_eval_on_nodes_release(det_h_proj);
  gkyl_eval_on_nodes_release(hamil_proj);
  gkyl_dg_calc_canonical_pb_vars_release(calc_vars);
}


void 
info_h_ij_inv_1x2v(double t, const double* xn, double* fout, void* ctx)
{
  // Inverse metric tensor, must be symmetric!
  // [h^{xx},h^{xy},h^{yy}]
  double R = 1.0;
  double q_theta = xn[0];
  const double q[1] = {q_theta};

  // [h^{thetatheta},h^{thetaphi},h^{phiphi}]
  fout[0] = 1.0 / pow(R, 2);
  fout[1] = 0.0;
  fout[2] = 1.0 / pow(R * sin(q[0]), 2);
}

void 
info_det_h_1x2v(double t, const double* xn, double* fout, void* ctx)
{
  // determinant of the metric tensor: J = det(h_{ij})
  double R = 1.0;
  double q_theta = xn[0];
  const double q[1] = {q_theta};
  fout[0] = pow(R, 2)*sin(q[0]);
}

void 
info_hamil_1x2v(double t, const double* xn, double* fout, void* ctx)
{
  // Canonical coordinates:
  double w0 = xn[1], w1 = xn[2];
  const double w[2] = {w0, w1};
  struct kh_2d_ctx *app = (struct kh_2d_ctx *)ctx;
  double *h_inv = malloc(3 * sizeof(double));
  info_h_ij_inv_1x2v(t, xn, h_inv, ctx); 
  fout[0] = 0.5 * h_inv[0] * w[0] * w[0] + 
            0.5 * (2.0* h_inv[1] * w[1] * w[0]) + 
            0.5 * h_inv[2] * w[1] * w[1];
  free(h_inv);
}

void 
test_1x2v(int poly_order, enum gkyl_basis_type b_type)
{
  double pi = 3.14159265359;
  double lower[] = {pi/4, -5.0, -5.0}, upper[] = {pi/2, 5.0, 5.0};
  int cells[] = {4, 4, 4};
  int vdim = 2, cdim = 1;
  int pdim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1], lower[2]}, velUpper[] = {upper[1], upper[2]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1], cells[2]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0, 0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    gkyl_cart_modal_hybrid(&basis, cdim, vdim);
  } 
  else if (b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  }
  else {
    gkyl_cart_modal_tensor(&basis, pdim, poly_order);;
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext; 
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0, 0};
  struct gkyl_range local, local_ext; 
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Allocate arrays for specified metric inverse, hamiltonian and metric determinant
  struct gkyl_array *h_ij_inv, *det_h, *hamil;
  h_ij_inv = mkarr(false, confBasis.num_basis*vdim*(vdim+1)/2, local_ext.volume);
  det_h = mkarr(false, confBasis.num_basis, local_ext.volume);
  hamil = mkarr(false, basis.num_basis, local_ext.volume);

  // Evaluate specified inverse metric function and det. at nodes to insure continuity
  struct gkyl_eval_on_nodes* h_ij_inv_proj = gkyl_eval_on_nodes_new(&confGrid, &confBasis, vdim*(vdim+1)/2, info_h_ij_inv_1x2v, 0);
  struct gkyl_eval_on_nodes* det_h_proj = gkyl_eval_on_nodes_new(&confGrid, &confBasis, 1, info_det_h_1x2v, 0);
  struct gkyl_eval_on_nodes* hamil_proj = gkyl_eval_on_nodes_new(&grid, &basis, 1, info_hamil_1x2v, 0);
  gkyl_eval_on_nodes_advance(h_ij_inv_proj, 0.0, &confLocal, h_ij_inv);
  gkyl_eval_on_nodes_advance(det_h_proj, 0.0, &confLocal, det_h);
  gkyl_eval_on_nodes_advance(hamil_proj, 0.0, &local_ext, hamil);

  // Need to figure out size of alpha_surf and sgn_alpha_surf by finding size of surface basis set 
  struct gkyl_basis surf_basis, surf_quad_basis;
  if (basis.b_type == GKYL_BASIS_MODAL_HYBRID) {
    // NOTE: If we are hybrid, allocate more memory than we need to avoid surface basis are different
    // sizes in each direction.
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, 2);
    gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, 2);
  } 
  else {
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, poly_order);
    gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, poly_order);
  }

  // always vdim + cdim
  int alpha_surf_sz = (vdim + cdim)*surf_basis.num_basis; 
  int sgn_alpha_surf_sz = (vdim + cdim)*surf_quad_basis.num_basis; // sign(alpha) is store at quadrature points

  // allocate arrays to store fields: 
  // 1. alpha_surf (surface phase space velocity)
  // 2. sgn_alpha_surf (sign(alpha_surf) at quadrature points)
  // 3. const_sgn_alpha (boolean for if sign(alpha_surf) is a constant, either +1 or -1)
  struct gkyl_array *alpha_surf = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *sgn_alpha_surf = mkarr(false, sgn_alpha_surf_sz, local_ext.volume);
  struct gkyl_array *const_sgn_alpha = mk_int_arr(false, (vdim + cdim), local_ext.volume);

  // Pre-compute alpha_surf, sgn_alpha_surf, const_sgn_alpha, and cot_vec since they are time-independent
  struct gkyl_dg_calc_canonical_pb_vars *calc_vars = gkyl_dg_calc_canonical_pb_vars_new(&grid, 
    &confBasis, &basis, false);
  gkyl_dg_calc_canonical_pb_vars_alpha_surf(calc_vars, &confLocal, &local, &local_ext, hamil,
    alpha_surf, sgn_alpha_surf, const_sgn_alpha);

  // Check continuity (Directly via Kernels, Option 2)
  double w_edge_not_used[3] = { 0.0, 0.0, 0.0 };
  double dxv[3] = {confGrid.dx[0], vel_grid.dx[0], vel_grid.dx[1] };
  struct gkyl_array *alpha_surf_comp_L = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *alpha_surf_comp_R = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *sgn_alpha_surf_comp_not_used = mkarr(false, sgn_alpha_surf_sz, local_ext.volume);

  // Loop over cells, get expansions in cells
  struct gkyl_range_iter piter;
  gkyl_range_iter_init(&piter, &local_ext);
  while(gkyl_range_iter_next(&piter)){
    long pidx = gkyl_range_idx(&local_ext, piter.idx);
    const double *hamil_local = gkyl_array_cfetch(hamil,pidx);
    double *alpha_surf_comp_L_local = gkyl_array_fetch(alpha_surf_comp_L,pidx);
    double *alpha_surf_comp_R_local = gkyl_array_fetch(alpha_surf_comp_R,pidx);
    double *sgn_alpha_surf_comp_not_used_local = gkyl_array_fetch(sgn_alpha_surf_comp_not_used,pidx);

    if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_1x2v_ser_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_1x2v_ser_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfvx_1x2v_ser_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
    } 
    else if (b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_1x2v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_1x2v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfvx_1x2v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
    }
    else if ((b_type == GKYL_BASIS_MODAL_TENSOR) && (poly_order == 1)) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_1x2v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_1x2v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfvx_1x2v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
    } 
    else if ((b_type == GKYL_BASIS_MODAL_TENSOR) && (poly_order == 2)) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_1x2v_tensor_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_1x2v_tensor_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfvx_1x2v_tensor_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
    } 
    else {
      assert(false);      
    }
  }

  // Compare suface wise element 
  // Loop over Number of surface coeff
  gkyl_range_iter_init(&piter, &local);
  int pidxl[GKYL_MAX_DIM];
  while(gkyl_range_iter_next(&piter)){
    long pidx = gkyl_range_idx(&local, piter.idx);
    const double *alpha_surf_comp_L_local = gkyl_array_cfetch(alpha_surf_comp_L, pidx);

    // Iterate in the lower direction (only conf space comparison - becuase there are no velocity space ghost cells for comparison)
    for (int dir = 0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(pdim, piter.idx, pidxl);
      pidxl[dir] = pidxl[dir]-1;
      long loc_l = gkyl_range_idx(&local, pidxl);
      const double *alpha_surf_comp_R_local = gkyl_array_cfetch(alpha_surf_comp_R, loc_l);

      // Iterate overthe number of basis on the surface
      for (int n = 0; n<surf_basis.num_basis; n++){
        TEST_CHECK(gkyl_compare_double(alpha_surf_comp_L_local[n + dir*surf_basis.num_basis], alpha_surf_comp_R_local[n + dir*surf_basis.num_basis], 1e-12));
      }
    }
  }

  // Release 
  gkyl_array_release(h_ij_inv);
  gkyl_array_release(det_h);
  gkyl_array_release(hamil);
  gkyl_array_release(alpha_surf);
  gkyl_array_release(sgn_alpha_surf);
  gkyl_array_release(const_sgn_alpha);
  gkyl_array_release(alpha_surf_comp_L);
  gkyl_array_release(alpha_surf_comp_R);
  gkyl_array_release(sgn_alpha_surf_comp_not_used);
  gkyl_eval_on_nodes_release(h_ij_inv_proj);
  gkyl_eval_on_nodes_release(det_h_proj);
  gkyl_eval_on_nodes_release(hamil_proj);
  gkyl_dg_calc_canonical_pb_vars_release(calc_vars);
}



void 
info_h_ij_inv_1x3v(double t, const double* xn, double* fout, void* ctx)
{
  // Inverse metric tensor, must be symmetric!
  double q_r = xn[0];

  // [h^{rr}, h^{rtheta},h^{rphi},h^{thetatheta},h^{thetaphi},h^{phiphi}]
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 1.0 / pow(q_r, 2);
  fout[4] = 0.0;
  fout[5] = 1.0;
}

void 
info_det_h_1x3v(double t, const double* xn, double* fout, void* ctx)
{
  // determinant of the metric tensor: J = det(h_{ij})
  double q_r = xn[0];
  fout[0] = q_r;
}

void 
info_hamil_1x3v(double t, const double* xn, double* fout, void* ctx)
{
  // Canonical coordinates:
  double q_R = xn[0], p_R_dot = xn[1], p_theta_dot = xn[2], p_phi_dot = xn[3];
  const double w[3] = {p_R_dot, p_theta_dot, p_phi_dot};
  double *h_inv = malloc(6 * sizeof(double));
  info_h_ij_inv_1x3v(t, xn, h_inv, ctx); 
  fout[0] = 0.5 * h_inv[0] * w[0] * w[0] + 
            0.5 * (2.0* h_inv[1] * w[1] * w[0]) + 
            0.5 * (2.0* h_inv[2] * w[2] * w[0]) +
            0.5 * h_inv[3] * w[1] * w[1] + 
            0.5 * (2.0* h_inv[4] * w[1] * w[2]) + 
            0.5 * h_inv[5] * w[2] * w[2];
  free(h_inv);
}

void 
test_1x3v(int poly_order, enum gkyl_basis_type b_type)
{
  double pi = 3.14159265359;
  double lower[] = { 0.5, -5.0, -5.0, -5.0 }, upper[] = { 1.5, 5.0, 5.0, 5.0 };
  int cells[] = {4, 4, 4, 4};
  int vdim = 3, cdim = 1;
  int pdim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1], lower[2], lower[3]}, velUpper[] = {upper[1], upper[2], upper[3]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1], cells[2], cells[3]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0, 0, 0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    gkyl_cart_modal_hybrid(&basis, cdim, vdim);
  } 
  else if (b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  }
  else {
    gkyl_cart_modal_tensor(&basis, pdim, poly_order);;
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext; 
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0, 0, 0};
  struct gkyl_range local, local_ext; 
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Allocate arrays for specified metric inverse, hamiltonian and metric determinant
  struct gkyl_array *h_ij_inv, *det_h, *hamil;
  h_ij_inv = mkarr(false, confBasis.num_basis*vdim*(vdim+1)/2, local_ext.volume);
  det_h = mkarr(false, confBasis.num_basis, local_ext.volume);
  hamil = mkarr(false, basis.num_basis, local_ext.volume);

  // Evaluate specified inverse metric function and det. at nodes to insure continuity
  struct gkyl_eval_on_nodes* h_ij_inv_proj = gkyl_eval_on_nodes_new(&confGrid, &confBasis, vdim*(vdim+1)/2, info_h_ij_inv_1x3v, 0);
  struct gkyl_eval_on_nodes* det_h_proj = gkyl_eval_on_nodes_new(&confGrid, &confBasis, 1, info_det_h_1x3v, 0);
  struct gkyl_eval_on_nodes* hamil_proj = gkyl_eval_on_nodes_new(&grid, &basis, 1, info_hamil_1x3v, 0);
  gkyl_eval_on_nodes_advance(h_ij_inv_proj, 0.0, &confLocal, h_ij_inv);
  gkyl_eval_on_nodes_advance(det_h_proj, 0.0, &confLocal, det_h);
  gkyl_eval_on_nodes_advance(hamil_proj, 0.0, &local_ext, hamil);

  // Need to figure out size of alpha_surf and sgn_alpha_surf by finding size of surface basis set 
  struct gkyl_basis surf_basis, surf_quad_basis;
  if (basis.b_type == GKYL_BASIS_MODAL_HYBRID) {
    // NOTE: If we are hybrid, allocate more memory than we need to avoid surface basis are different
    // sizes in each direction.
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, 2);
    gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, 2);
  } 
  else {
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, poly_order);
    gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, poly_order);
  }

  // always vdim + cdim
  int alpha_surf_sz = (vdim + cdim)*surf_basis.num_basis; 
  int sgn_alpha_surf_sz = (vdim + cdim)*surf_quad_basis.num_basis; // sign(alpha) is store at quadrature points

  // allocate arrays to store fields: 
  // 1. alpha_surf (surface phase space velocity)
  // 2. sgn_alpha_surf (sign(alpha_surf) at quadrature points)
  // 3. const_sgn_alpha (boolean for if sign(alpha_surf) is a constant, either +1 or -1)
  struct gkyl_array *alpha_surf = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *sgn_alpha_surf = mkarr(false, sgn_alpha_surf_sz, local_ext.volume);
  struct gkyl_array *const_sgn_alpha = mk_int_arr(false, (vdim + cdim), local_ext.volume);

  // Pre-compute alpha_surf, sgn_alpha_surf, const_sgn_alpha, and cot_vec since they are time-independent
  struct gkyl_dg_calc_canonical_pb_vars *calc_vars = gkyl_dg_calc_canonical_pb_vars_new(&grid, 
    &confBasis, &basis, false);
  gkyl_dg_calc_canonical_pb_vars_alpha_surf(calc_vars, &confLocal, &local, &local_ext, hamil,
    alpha_surf, sgn_alpha_surf, const_sgn_alpha);

  // Check continuity (Directly via Kernels, Option 2)
  double w_edge_not_used[4] = { 0.0, 0.0, 0.0, 0.0 };
  double dxv[4] = {confGrid.dx[0], vel_grid.dx[0], vel_grid.dx[1], vel_grid.dx[2] };
  struct gkyl_array *alpha_surf_comp_L = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *alpha_surf_comp_R = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *sgn_alpha_surf_comp_not_used = mkarr(false, sgn_alpha_surf_sz, local_ext.volume);

  // Loop over cells, get expansions in cells
  struct gkyl_range_iter piter;
  gkyl_range_iter_init(&piter, &local_ext);
  while(gkyl_range_iter_next(&piter)){
    long pidx = gkyl_range_idx(&local_ext, piter.idx);
    const double *hamil_local = gkyl_array_cfetch(hamil,pidx);
    double *alpha_surf_comp_L_local = gkyl_array_fetch(alpha_surf_comp_L,pidx);
    double *alpha_surf_comp_R_local = gkyl_array_fetch(alpha_surf_comp_R,pidx);
    double *sgn_alpha_surf_comp_not_used_local = gkyl_array_fetch(sgn_alpha_surf_comp_not_used,pidx);

    if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_1x3v_ser_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_1x3v_ser_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfvx_1x3v_ser_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
    } 
    else if (b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_1x3v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_1x3v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfvx_1x3v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
    }
    else if ((b_type == GKYL_BASIS_MODAL_TENSOR) && (poly_order == 1)) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_1x3v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_1x3v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfvx_1x3v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
    } 
    else if ((b_type == GKYL_BASIS_MODAL_TENSOR) && (poly_order == 2)) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_1x3v_tensor_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_1x3v_tensor_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfvx_1x3v_tensor_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
    } 
    else {
      assert(false);
    }
  }

  // Compare suface wise element 
  // Loop over Number of surface coeff
  gkyl_range_iter_init(&piter, &local);
  int pidxl[GKYL_MAX_DIM];
  while(gkyl_range_iter_next(&piter)){
    long pidx = gkyl_range_idx(&local, piter.idx);
    const double *alpha_surf_comp_L_local = gkyl_array_cfetch(alpha_surf_comp_L, pidx);

    // Iterate in the lower direction (only conf space comparison - becuase there are no velocity space ghost cells for comparison)
    for (int dir = 0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(pdim, piter.idx, pidxl);
      pidxl[dir] = pidxl[dir]-1;
      long loc_l = gkyl_range_idx(&local, pidxl);
      const double *alpha_surf_comp_R_local = gkyl_array_cfetch(alpha_surf_comp_R, loc_l);

      // Iterate overthe number of basis on the surface
      for (int n = 0; n<surf_basis.num_basis; n++){
        TEST_CHECK(gkyl_compare_double(alpha_surf_comp_L_local[n + dir*surf_basis.num_basis], alpha_surf_comp_R_local[n + dir*surf_basis.num_basis], 1e-12));
      }
    }
  }

  // Release 
  gkyl_array_release(h_ij_inv);
  gkyl_array_release(det_h);
  gkyl_array_release(hamil);
  gkyl_array_release(alpha_surf);
  gkyl_array_release(sgn_alpha_surf);
  gkyl_array_release(const_sgn_alpha);
  gkyl_array_release(alpha_surf_comp_L);
  gkyl_array_release(alpha_surf_comp_R);
  gkyl_array_release(sgn_alpha_surf_comp_not_used);
  gkyl_eval_on_nodes_release(h_ij_inv_proj);
  gkyl_eval_on_nodes_release(det_h_proj);
  gkyl_eval_on_nodes_release(hamil_proj);
  gkyl_dg_calc_canonical_pb_vars_release(calc_vars);
}



void 
info_h_ij_inv_2x2v(double t, const double* xn, double* fout, void* ctx)
{
  // Inverse metric tensor, must be symmetric!
  // [h^{xx},h^{xy},h^{yy}]
  double R = 1.0;
  double q_theta = xn[0], q_phi = xn[1];
  const double q[2] = {q_theta, q_phi};

  // [h^{thetatheta},h^{thetaphi},h^{phiphi}]
  fout[0] = 1.0 / pow(R, 2);
  fout[1] = 0.0;
  fout[2] = 1.0 / pow(R * sin(q[0]), 2);
}

void 
info_det_h_2x2v(double t, const double* xn, double* fout, void* ctx)
{
  // determinant of the metric tensor: J = det(h_{ij})
  double R = 1.0;
  double q_theta = xn[0], q_phi = xn[1];
  const double q[2] = {q_theta, q_phi};
  fout[0] = pow(R, 2)*sin(q[0]);
}

void 
info_hamil_2x2v(double t, const double* xn, double* fout, void* ctx)
{
  // Canonical coordinates:
  double w0 = xn[2], w1 = xn[3];
  const double w[2] = {w0, w1};
  struct kh_2d_ctx *app = (struct kh_2d_ctx *)ctx;
  double *h_inv = malloc(3 * sizeof(double));
  info_h_ij_inv_2x2v(t, xn, h_inv, ctx); 
  fout[0] = 0.5 * h_inv[0] * w[0] * w[0] + 
            0.5 * (2.0* h_inv[1] * w[1] * w[0]) + 
            0.5 * h_inv[2] * w[1] * w[1];
  free(h_inv);
}


void 
test_2x2v(int poly_order, enum gkyl_basis_type b_type)
{
  double pi = 3.14159265359;
  double lower[] = {pi/4,pi/4, -5.0, -5.0}, upper[] = {(1.01)*pi/4,(1.01)*pi/4, 5.0, 5.0};
  int cells[] = {4, 4, 4, 4};
  int vdim = 2, cdim = 2;
  int pdim = cdim + vdim;

  double confLower[] = {lower[0], lower[1]}, confUpper[] = {upper[0], upper[1]};
  double velLower[] = {lower[2], lower[3]}, velUpper[] = {upper[2], upper[3]};
  int confCells[] = {cells[0], cells[1]};
  int velCells[] = {cells[2], cells[3]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0, 0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    gkyl_cart_modal_hybrid(&basis, cdim, vdim);
  } 
  else if (b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  }
  else {
    gkyl_cart_modal_tensor(&basis, pdim, poly_order);
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = {1, 1};
  struct gkyl_range confLocal, confLocal_ext; 
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], confGhost[1], 0, 0};
  struct gkyl_range local, local_ext; 
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Allocate arrays for specified metric inverse, hamiltonian and metric determinant
  struct gkyl_array *h_ij_inv, *det_h, *hamil;
  h_ij_inv = mkarr(false, confBasis.num_basis*vdim*(vdim+1)/2, local_ext.volume);
  det_h = mkarr(false, confBasis.num_basis, local_ext.volume);
  hamil = mkarr(false, basis.num_basis, local_ext.volume);

  // Evaluate specified inverse metric function and det. at nodes to insure continuity
  struct gkyl_eval_on_nodes* h_ij_inv_proj = gkyl_eval_on_nodes_new(&confGrid, &confBasis, vdim*(vdim+1)/2, info_h_ij_inv_2x2v, 0);
  struct gkyl_eval_on_nodes* det_h_proj = gkyl_eval_on_nodes_new(&confGrid, &confBasis, 1, info_det_h_2x2v, 0);
  struct gkyl_eval_on_nodes* hamil_proj = gkyl_eval_on_nodes_new(&grid, &basis, 1, info_hamil_2x2v, 0);
  gkyl_eval_on_nodes_advance(h_ij_inv_proj, 0.0, &confLocal, h_ij_inv);
  gkyl_eval_on_nodes_advance(det_h_proj, 0.0, &confLocal, det_h);
  gkyl_eval_on_nodes_advance(hamil_proj, 0.0, &local_ext, hamil);

  // Need to figure out size of alpha_surf and sgn_alpha_surf by finding size of surface basis set 
  struct gkyl_basis surf_basis, surf_quad_basis;
  if (basis.b_type == GKYL_BASIS_MODAL_HYBRID) {
    // NOTE: If we are hybrid, allocate more memory than we need to avoid surface basis are different
    // sizes in each direction. (!) Will need to extend this based on direction
    //gkyl_cart_modal_serendip(&surf_basis, pdim-1, 2);
    gkyl_cart_modal_hybrid(&surf_basis, cdim, vdim-1);
    gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, 2);
  } 
  else {
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, poly_order);
    gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, poly_order);
  }

  // always vdim + cdim
  int alpha_surf_sz = (vdim + cdim)*surf_basis.num_basis; 
  int sgn_alpha_surf_sz = (vdim + cdim)*surf_quad_basis.num_basis; // sign(alpha) is store at quadrature points

  // allocate arrays to store fields: 
  // 1. alpha_surf (surface phase space velocity)
  // 2. sgn_alpha_surf (sign(alpha_surf) at quadrature points)
  // 3. const_sgn_alpha (boolean for if sign(alpha_surf) is a constant, either +1 or -1)
  struct gkyl_array *alpha_surf = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *sgn_alpha_surf = mkarr(false, sgn_alpha_surf_sz, local_ext.volume);
  struct gkyl_array *const_sgn_alpha = mk_int_arr(false, (vdim + cdim), local_ext.volume);

  // Pre-compute alpha_surf, sgn_alpha_surf, const_sgn_alpha, and cot_vec since they are time-independent
  struct gkyl_dg_calc_canonical_pb_vars *calc_vars = gkyl_dg_calc_canonical_pb_vars_new(&grid, 
    &confBasis, &basis, false);
  gkyl_dg_calc_canonical_pb_vars_alpha_surf(calc_vars, &confLocal, &local, &local_ext, hamil,
    alpha_surf, sgn_alpha_surf, const_sgn_alpha);

  // Check continuity (Directly via Kernels, Option 2)
  double w_edge_not_used[4] = { 0.0, 0.0, 0.0, 0.0 };
  double dxv[4] = {confGrid.dx[0], confGrid.dx[1], vel_grid.dx[0], vel_grid.dx[1]};
  struct gkyl_array *alpha_surf_comp_L = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *alpha_surf_comp_R = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *sgn_alpha_surf_comp_not_used = mkarr(false, sgn_alpha_surf_sz, local_ext.volume);

  // Loop over cells, get expansions in cells
  struct gkyl_range_iter piter;
  gkyl_range_iter_init(&piter, &local_ext);
  while(gkyl_range_iter_next(&piter)){
    long pidx = gkyl_range_idx(&local_ext, piter.idx);
    const double *hamil_local = gkyl_array_cfetch(hamil,pidx);
    double *alpha_surf_comp_L_local = gkyl_array_fetch(alpha_surf_comp_L,pidx);
    double *alpha_surf_comp_R_local = gkyl_array_fetch(alpha_surf_comp_R,pidx);
    double *sgn_alpha_surf_comp_not_used_local = gkyl_array_fetch(sgn_alpha_surf_comp_not_used,pidx);

    if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_2x2v_ser_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_2x2v_ser_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfy_2x2v_ser_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfy_2x2v_ser_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
    } 
    else if (b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_2x2v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_2x2v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfy_2x2v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfy_2x2v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
    }
    else if ((b_type == GKYL_BASIS_MODAL_TENSOR) && (poly_order == 1)) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_2x2v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_2x2v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfy_2x2v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfy_2x2v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
    } 
    else if ((b_type == GKYL_BASIS_MODAL_TENSOR) && (poly_order == 2)) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_2x2v_tensor_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_2x2v_tensor_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfy_2x2v_tensor_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfy_2x2v_tensor_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
    } 
    else {
      assert(false);
    }
  }

  // Compare suface wise element 
  // Loop over Number of surface coeff
  gkyl_range_iter_init(&piter, &local);
  int pidxl[GKYL_MAX_DIM];
  while(gkyl_range_iter_next(&piter)){
    long pidx = gkyl_range_idx(&local, piter.idx);
    const double *alpha_surf_comp_L_local = gkyl_array_cfetch(alpha_surf_comp_L, pidx);

    // Iterate in the lower direction (only conf space comparison - becuase there are no velocity space ghost cells for comparison)
    for (int dir = 0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(pdim, piter.idx, pidxl);
      pidxl[dir] = pidxl[dir]-1;
      long loc_l = gkyl_range_idx(&local, pidxl);
      const double *alpha_surf_comp_R_local = gkyl_array_cfetch(alpha_surf_comp_R, loc_l);

      // Iterate overthe number of basis on the surface
      for (int n = 0; n<surf_basis.num_basis; n++){
        TEST_CHECK(gkyl_compare_double(alpha_surf_comp_L_local[n + dir*surf_basis.num_basis], alpha_surf_comp_R_local[n + dir*surf_basis.num_basis], 1e-9));
      }
    }
  }

  // Release 
  gkyl_array_release(h_ij_inv);
  gkyl_array_release(det_h);
  gkyl_array_release(hamil);
  gkyl_array_release(alpha_surf);
  gkyl_array_release(sgn_alpha_surf);
  gkyl_array_release(const_sgn_alpha);
  gkyl_array_release(alpha_surf_comp_L);
  gkyl_array_release(alpha_surf_comp_R);
  gkyl_array_release(sgn_alpha_surf_comp_not_used);
  gkyl_eval_on_nodes_release(h_ij_inv_proj);
  gkyl_eval_on_nodes_release(det_h_proj);
  gkyl_eval_on_nodes_release(hamil_proj);
  gkyl_dg_calc_canonical_pb_vars_release(calc_vars);
}

void 
info_h_ij_inv_2x3v(double t, const double* xn, double* fout, void* ctx)
{
  // Inverse metric tensor, must be symmetric!
  double q_r = xn[0], q_theta = xn[1];

  // [h^{rr}, h^{rtheta},h^{rphi},h^{thetatheta},h^{thetaphi},h^{phiphi}]
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 1.0 / pow(q_r, 2);
  fout[4] = 0.0;
  fout[5] = 1.0 / pow((2.0 + q_r*cos(q_theta)),2);
}

void 
info_det_h_2x3v(double t, const double* xn, double* fout, void* ctx)
{
  // determinant of the metric tensor: J = det(h_{ij})
  double q_r = xn[0], q_theta = xn[1];
  fout[0] = q_r*(2.0 + q_r*cos(q_theta));
}

void 
info_hamil_2x3v(double t, const double* xn, double* fout, void* ctx)
{
  // Canonical coordinates:
  double q_R = xn[0], q_theta = xn[1], p_R_dot = xn[2], p_theta_dot = xn[3], p_phi_dot = xn[4];
  const double w[3] = {p_R_dot, p_theta_dot, p_phi_dot};
  double *h_inv = malloc(6 * sizeof(double));
  info_h_ij_inv_2x3v(t, xn, h_inv, ctx); 
  fout[0] = 0.5 * h_inv[0] * w[0] * w[0] + 
            0.5 * (2.0* h_inv[1] * w[1] * w[0]) + 
            0.5 * (2.0* h_inv[2] * w[2] * w[0]) +
            0.5 * h_inv[3] * w[1] * w[1] + 
            0.5 * (2.0* h_inv[4] * w[1] * w[2]) + 
            0.5 * h_inv[5] * w[2] * w[2];
  free(h_inv);
}


void 
test_2x3v(int poly_order, enum gkyl_basis_type b_type)
{
  double pi = 3.14159265359;
  double lower[] = { 0.5, 0.0, -5.0, -5.0, -5.0 }, upper[] = { 1.5, 2.0*pi, 5.0, 5.0, 5.0 };
  int cells[] = {4, 4, 4, 4, 4};
  int vdim = 3, cdim = 2;
  int pdim = cdim + vdim;

  double confLower[] = {lower[0], lower[1]}, confUpper[] = {upper[0], upper[1]};
  double velLower[] = {lower[2], lower[3], lower[4]}, velUpper[] = {upper[2], upper[3], upper[4]};
  int confCells[] = {cells[0], cells[1]};
  int velCells[] = {cells[2], cells[3], cells[4]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0, 0, 0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    gkyl_cart_modal_hybrid(&basis, cdim, vdim);
  } 
  else if (b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  }
  else {
    gkyl_cart_modal_tensor(&basis, pdim, poly_order);
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = {1, 1};
  struct gkyl_range confLocal, confLocal_ext; 
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], confGhost[1], 0, 0, 0};
  struct gkyl_range local, local_ext; 
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Allocate arrays for specified metric inverse, hamiltonian and metric determinant
  struct gkyl_array *h_ij_inv, *det_h, *hamil;
  h_ij_inv = mkarr(false, confBasis.num_basis*vdim*(vdim+1)/2, local_ext.volume);
  det_h = mkarr(false, confBasis.num_basis, local_ext.volume);
  hamil = mkarr(false, basis.num_basis, local_ext.volume);

  // Evaluate specified inverse metric function and det. at nodes to insure continuity
  struct gkyl_eval_on_nodes* h_ij_inv_proj = gkyl_eval_on_nodes_new(&confGrid, &confBasis, vdim*(vdim+1)/2, info_h_ij_inv_2x3v, 0);
  struct gkyl_eval_on_nodes* det_h_proj = gkyl_eval_on_nodes_new(&confGrid, &confBasis, 1, info_det_h_2x3v, 0);
  struct gkyl_eval_on_nodes* hamil_proj = gkyl_eval_on_nodes_new(&grid, &basis, 1, info_hamil_2x3v, 0);
  gkyl_eval_on_nodes_advance(h_ij_inv_proj, 0.0, &confLocal, h_ij_inv);
  gkyl_eval_on_nodes_advance(det_h_proj, 0.0, &confLocal, det_h);
  gkyl_eval_on_nodes_advance(hamil_proj, 0.0, &local_ext, hamil);

  // Need to figure out size of alpha_surf and sgn_alpha_surf by finding size of surface basis set 
  struct gkyl_basis surf_basis, surf_quad_basis;
  if (basis.b_type == GKYL_BASIS_MODAL_HYBRID) {
    // NOTE: If we are hybrid, allocate more memory than we need to avoid surface basis are different
    // sizes in each direction.
    // ** Assumes conf-face comparison in this test only ** 
    gkyl_cart_modal_hybrid(&surf_basis, cdim-1, vdim); 
    //gkyl_cart_modal_serendip(&surf_basis, pdim-1, 2);
    gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, 2);
  } 
  else {
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, poly_order);
    gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, poly_order);
  }

  // always vdim + cdim
  int alpha_surf_sz = (vdim + cdim)*surf_basis.num_basis; 
  int sgn_alpha_surf_sz = (vdim + cdim)*surf_quad_basis.num_basis; // sign(alpha) is store at quadrature points

  // allocate arrays to store fields: 
  // 1. alpha_surf (surface phase space velocity)
  // 2. sgn_alpha_surf (sign(alpha_surf) at quadrature points)
  // 3. const_sgn_alpha (boolean for if sign(alpha_surf) is a constant, either +1 or -1)
  struct gkyl_array *alpha_surf = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *sgn_alpha_surf = mkarr(false, sgn_alpha_surf_sz, local_ext.volume);
  struct gkyl_array *const_sgn_alpha = mk_int_arr(false, (vdim + cdim), local_ext.volume);

  // Pre-compute alpha_surf, sgn_alpha_surf, const_sgn_alpha, and cot_vec since they are time-independent
  struct gkyl_dg_calc_canonical_pb_vars *calc_vars = gkyl_dg_calc_canonical_pb_vars_new(&grid, 
    &confBasis, &basis, false);
  gkyl_dg_calc_canonical_pb_vars_alpha_surf(calc_vars, &confLocal, &local, &local_ext, hamil,
    alpha_surf, sgn_alpha_surf, const_sgn_alpha);

  // Check continuity (Directly via Kernels, Option 2)
  double w_edge_not_used[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
  double dxv[5] = {confGrid.dx[0], confGrid.dx[1], vel_grid.dx[0], vel_grid.dx[1], vel_grid.dx[2]};
  struct gkyl_array *alpha_surf_comp_L = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *alpha_surf_comp_R = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *sgn_alpha_surf_comp_not_used = mkarr(false, sgn_alpha_surf_sz, local_ext.volume);

 // Loop over cells, get expansions in cells
  struct gkyl_range_iter piter;
  gkyl_range_iter_init(&piter, &local_ext);
  while(gkyl_range_iter_next(&piter)){
    long pidx = gkyl_range_idx(&local_ext, piter.idx);
    const double *hamil_local = gkyl_array_cfetch(hamil,pidx);
    double *alpha_surf_comp_L_local = gkyl_array_fetch(alpha_surf_comp_L,pidx);
    double *alpha_surf_comp_R_local = gkyl_array_fetch(alpha_surf_comp_R,pidx);
    double *sgn_alpha_surf_comp_not_used_local = gkyl_array_fetch(sgn_alpha_surf_comp_not_used,pidx);

    if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_2x3v_ser_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_2x3v_ser_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfy_2x3v_ser_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfy_2x3v_ser_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
    } 
    else if (b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_2x3v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_2x3v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfy_2x3v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfy_2x3v_ser_p2(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
    }
    else if ((b_type == GKYL_BASIS_MODAL_TENSOR) && (poly_order == 1)) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_2x3v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_2x3v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfy_2x3v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfy_2x3v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
    } 
    else {
      //Tensor p2 is not availible for 2x3v
      assert(false);
    }  
  }

  // Compare suface wise element 
  // Loop over Number of surface coeff
  gkyl_range_iter_init(&piter, &local);
  int pidxl[GKYL_MAX_DIM];
  while(gkyl_range_iter_next(&piter)){
    long pidx = gkyl_range_idx(&local, piter.idx);
    const double *alpha_surf_comp_L_local = gkyl_array_cfetch(alpha_surf_comp_L, pidx);

    // Iterate in the lower direction (only conf space comparison - becuase there are no velocity space ghost cells for comparison)
    for (int dir = 0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(pdim, piter.idx, pidxl);
      pidxl[dir] = pidxl[dir]-1;
      long loc_l = gkyl_range_idx(&local, pidxl);
      const double *alpha_surf_comp_R_local = gkyl_array_cfetch(alpha_surf_comp_R, loc_l);

      // Iterate overthe number of basis on the surface
      for (int n = 0; n<surf_basis.num_basis; n++){
        TEST_CHECK(gkyl_compare_double(alpha_surf_comp_L_local[n + dir*surf_basis.num_basis], alpha_surf_comp_R_local[n + dir*surf_basis.num_basis], 1e-10));
      }
    }
  }

  // Release 
  gkyl_array_release(h_ij_inv);
  gkyl_array_release(det_h);
  gkyl_array_release(hamil);
  gkyl_array_release(alpha_surf);
  gkyl_array_release(sgn_alpha_surf);
  gkyl_array_release(const_sgn_alpha);
  gkyl_array_release(alpha_surf_comp_L);
  gkyl_array_release(alpha_surf_comp_R);
  gkyl_array_release(sgn_alpha_surf_comp_not_used);
  gkyl_eval_on_nodes_release(h_ij_inv_proj);
  gkyl_eval_on_nodes_release(det_h_proj);
  gkyl_eval_on_nodes_release(hamil_proj);
  gkyl_dg_calc_canonical_pb_vars_release(calc_vars);
}


void 
info_h_ij_inv_3x3v(double t, const double* xn, double* fout, void* ctx)
{
  // Inverse metric tensor, must be symmetric!
  double q_r = xn[0], q_theta = xn[1], q_phi = xn[2];

  // [h^{rr}, h^{rtheta},h^{rphi},h^{thetatheta},h^{thetaphi},h^{phiphi}]
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 1.0 / pow(q_r, 2);
  fout[4] = 0.0;
  fout[5] = 1.0 / pow((2.0 + q_r*cos(q_theta)),2);
}

void 
info_det_h_3x3v(double t, const double* xn, double* fout, void* ctx)
{
  // determinant of the metric tensor: J = det(h_{ij})
  double q_r = xn[0], q_theta = xn[1], q_phi = xn[2];
  fout[0] = q_r*(2.0 + q_r*cos(q_theta));
}

void 
info_hamil_3x3v(double t, const double* xn, double* fout, void* ctx)
{
  // Canonical coordinates:
  double q_R = xn[0], q_theta = xn[1], q_phi = xn[2], p_R_dot = xn[3], p_theta_dot = xn[4], p_phi_dot = xn[5];
  const double w[3] = {p_R_dot, p_theta_dot, p_phi_dot};
  double *h_inv = malloc(6 * sizeof(double));
  info_h_ij_inv_3x3v(t, xn, h_inv, ctx); 
  fout[0] = 0.5 * h_inv[0] * w[0] * w[0] + 
            0.5 * (2.0* h_inv[1] * w[1] * w[0]) + 
            0.5 * (2.0* h_inv[2] * w[2] * w[0]) +
            0.5 * h_inv[3] * w[1] * w[1] + 
            0.5 * (2.0* h_inv[4] * w[1] * w[2]) + 
            0.5 * h_inv[5] * w[2] * w[2];
  free(h_inv);
}


void 
test_3x3v(int poly_order, enum gkyl_basis_type b_type)
{
  double pi = 3.14159265359;
  double lower[] = { 0.5, 0.0, 0.0, -5.0, -5.0, -5.0 }, upper[] = { 1.5, 2.0*pi, 2.0*pi, 5.0, 5.0, 5.0 };
  int cells[] = {4, 4, 4, 4, 4, 4};
  int vdim = 3, cdim = 3;
  int pdim = cdim + vdim;

  double confLower[] = {lower[0], lower[1], lower[2]}, confUpper[] = {upper[0], upper[1], upper[2]};
  double velLower[] = {lower[3], lower[4], lower[5]}, velUpper[] = {upper[3], upper[4], upper[5]};
  int confCells[] = {cells[0], cells[1], cells[2]};
  int velCells[] = {cells[3], cells[4], cells[5]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0, 0, 0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1 && b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    gkyl_cart_modal_hybrid(&basis, cdim, vdim);
  } 
  else if (b_type == GKYL_BASIS_MODAL_SERENDIPITY) {
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  }
  else {
    gkyl_cart_modal_tensor(&basis, pdim, poly_order);
  }
  //gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  // (NOTE!) Need to change to tensor for 3x3v p1 to work:
  gkyl_cart_modal_tensor(&confBasis, cdim, poly_order);

  int confGhost[] = {1, 1, 1};
  struct gkyl_range confLocal, confLocal_ext; 
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], confGhost[1], confGhost[2], 0, 0, 0};
  struct gkyl_range local, local_ext; 
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Allocate arrays for specified metric inverse, hamiltonian and metric determinant
  struct gkyl_array *h_ij_inv, *det_h, *hamil;
  h_ij_inv = mkarr(false, confBasis.num_basis*vdim*(vdim+1)/2, local_ext.volume);
  det_h = mkarr(false, confBasis.num_basis, local_ext.volume);
  hamil = mkarr(false, basis.num_basis, local_ext.volume);

  // Evaluate specified inverse metric function and det. at nodes to insure continuity
  struct gkyl_eval_on_nodes* h_ij_inv_proj = gkyl_eval_on_nodes_new(&confGrid, &confBasis, vdim*(vdim+1)/2, info_h_ij_inv_3x3v, 0);
  struct gkyl_eval_on_nodes* det_h_proj = gkyl_eval_on_nodes_new(&confGrid, &confBasis, 1, info_det_h_3x3v, 0);
  struct gkyl_eval_on_nodes* hamil_proj = gkyl_eval_on_nodes_new(&grid, &basis, 1, info_hamil_3x3v, 0);
  gkyl_eval_on_nodes_advance(h_ij_inv_proj, 0.0, &confLocal, h_ij_inv);
  gkyl_eval_on_nodes_advance(det_h_proj, 0.0, &confLocal, det_h);
  gkyl_eval_on_nodes_advance(hamil_proj, 0.0, &local_ext, hamil);

  // Need to figure out size of alpha_surf and sgn_alpha_surf by finding size of surface basis set 
  struct gkyl_basis surf_basis, surf_quad_basis;
  if (basis.b_type == GKYL_BASIS_MODAL_HYBRID) {
    // NOTE: If we are hybrid, allocate more memory than we need to avoid surface basis are different
    // sizes in each direction.
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, 2);
    gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, 2);
  } 
  else {
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, poly_order);
    gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, poly_order);
  }

  // always vdim + cdim
  int alpha_surf_sz = (vdim + cdim)*surf_basis.num_basis; 
  int sgn_alpha_surf_sz = (vdim + cdim)*surf_quad_basis.num_basis; // sign(alpha) is store at quadrature points

  // allocate arrays to store fields: 
  // 1. alpha_surf (surface phase space velocity)
  // 2. sgn_alpha_surf (sign(alpha_surf) at quadrature points)
  // 3. const_sgn_alpha (boolean for if sign(alpha_surf) is a constant, either +1 or -1)
  struct gkyl_array *alpha_surf = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *sgn_alpha_surf = mkarr(false, sgn_alpha_surf_sz, local_ext.volume);
  struct gkyl_array *const_sgn_alpha = mk_int_arr(false, (vdim + cdim), local_ext.volume);

  // Pre-compute alpha_surf, sgn_alpha_surf, const_sgn_alpha, and cot_vec since they are time-independent
  struct gkyl_dg_calc_canonical_pb_vars *calc_vars = gkyl_dg_calc_canonical_pb_vars_new(&grid, 
    &confBasis, &basis, false);
  gkyl_dg_calc_canonical_pb_vars_alpha_surf(calc_vars, &confLocal, &local, &local_ext, hamil,
    alpha_surf, sgn_alpha_surf, const_sgn_alpha);

  // Check continuity (Directly via Kernels, Option 2)
  double w_edge_not_used[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  double dxv[6] = {confGrid.dx[0], confGrid.dx[1], confGrid.dx[2], vel_grid.dx[0], vel_grid.dx[1], vel_grid.dx[2]};
  struct gkyl_array *alpha_surf_comp_L = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *alpha_surf_comp_R = mkarr(false, alpha_surf_sz, local_ext.volume);
  struct gkyl_array *sgn_alpha_surf_comp_not_used = mkarr(false, sgn_alpha_surf_sz, local_ext.volume);

  // Loop over cells, get expansions in cells
  struct gkyl_range_iter piter;
  gkyl_range_iter_init(&piter, &local_ext);
  while(gkyl_range_iter_next(&piter)){
    long pidx = gkyl_range_idx(&local_ext, piter.idx);
    const double *hamil_local = gkyl_array_cfetch(hamil,pidx);
    double *alpha_surf_comp_L_local = gkyl_array_fetch(alpha_surf_comp_L,pidx);
    double *alpha_surf_comp_R_local = gkyl_array_fetch(alpha_surf_comp_R,pidx);
    double *sgn_alpha_surf_comp_not_used_local = gkyl_array_fetch(sgn_alpha_surf_comp_not_used,pidx);

    if ((b_type == GKYL_BASIS_MODAL_TENSOR) && (poly_order == 1)) {
      int const_sgn_alpha_surf = canonical_pb_alpha_surfx_3x3v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      int const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfx_3x3v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfy_3x3v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfy_3x3v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfz_3x3v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf_edge = canonical_pb_alpha_edge_surfz_3x3v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_R_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfvx_3x3v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local); 
      const_sgn_alpha_surf = canonical_pb_alpha_surfvy_3x3v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local);
      const_sgn_alpha_surf = canonical_pb_alpha_surfvz_3x3v_tensor_p1(w_edge_not_used, dxv, hamil_local, alpha_surf_comp_L_local, sgn_alpha_surf_comp_not_used_local);
    } 
    else {
      //Only tensor p1 is supported for 3x3v
      assert(false);
    }
  }

  // Compare suface wise element 
  // Loop over Number of surface coeff
  gkyl_range_iter_init(&piter, &local);
  int pidxl[GKYL_MAX_DIM];
  while(gkyl_range_iter_next(&piter)){
    long pidx = gkyl_range_idx(&local, piter.idx);
    const double *alpha_surf_comp_L_local = gkyl_array_cfetch(alpha_surf_comp_L, pidx);

    // Iterate in the lower direction (only conf space comparison - becuase there are no velocity space ghost cells for comparison)
    for (int dir = 0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(pdim, piter.idx, pidxl);
      pidxl[dir] = pidxl[dir]-1;
      long loc_l = gkyl_range_idx(&local, pidxl);
      const double *alpha_surf_comp_R_local = gkyl_array_cfetch(alpha_surf_comp_R, loc_l);

      // Iterate overthe number of basis on the surface
      for (int n = 0; n<surf_basis.num_basis; n++){
        TEST_CHECK(gkyl_compare_double(alpha_surf_comp_L_local[n + dir*surf_basis.num_basis], alpha_surf_comp_R_local[n + dir*surf_basis.num_basis], 1e-10));
      }
    }
  }

  // Release 
  gkyl_array_release(h_ij_inv);
  gkyl_array_release(det_h);
  gkyl_array_release(hamil);
  gkyl_array_release(alpha_surf);
  gkyl_array_release(sgn_alpha_surf);
  gkyl_array_release(const_sgn_alpha);
  gkyl_array_release(alpha_surf_comp_L);
  gkyl_array_release(alpha_surf_comp_R);
  gkyl_array_release(sgn_alpha_surf_comp_not_used);
  gkyl_eval_on_nodes_release(h_ij_inv_proj);
  gkyl_eval_on_nodes_release(det_h_proj);
  gkyl_eval_on_nodes_release(hamil_proj);
  gkyl_dg_calc_canonical_pb_vars_release(calc_vars);
}

// Check the 1x1v_p1/2 continuity
void test_1x1v_p1_continuity_tensor() { test_1x1v(1, GKYL_BASIS_MODAL_TENSOR); }
void test_1x1v_p1_continuity_ser() { test_1x1v(1, GKYL_BASIS_MODAL_SERENDIPITY); }
void test_1x1v_p2_continuity_tensor() { test_1x1v(2, GKYL_BASIS_MODAL_TENSOR); }
void test_1x1v_p2_continuity_ser() { test_1x1v(2, GKYL_BASIS_MODAL_SERENDIPITY); }
// Check the 1x2v_p1/2 continuity
void test_1x2v_p1_continuity_tensor() { test_1x2v(1, GKYL_BASIS_MODAL_TENSOR); }
void test_1x2v_p1_continuity_ser() { test_1x2v(1, GKYL_BASIS_MODAL_SERENDIPITY); }
void test_1x2v_p2_continuity_tensor() { test_1x2v(2, GKYL_BASIS_MODAL_TENSOR); }
void test_1x2v_p2_continuity_ser() { test_1x2v(2, GKYL_BASIS_MODAL_SERENDIPITY); }
// Check the 1x3v_p1/2 continuity
void test_1x3v_p1_continuity_tensor() { test_1x3v(1, GKYL_BASIS_MODAL_TENSOR); }
void test_1x3v_p1_continuity_ser() { test_1x3v(1, GKYL_BASIS_MODAL_SERENDIPITY); }
void test_1x3v_p2_continuity_tensor() { test_1x3v(2, GKYL_BASIS_MODAL_TENSOR); }
void test_1x3v_p2_continuity_ser() { test_1x3v(2, GKYL_BASIS_MODAL_SERENDIPITY); }
// Check the 2x2v_p1/2 continuity
void test_2x2v_p1_continuity_tensor() { test_2x2v(1, GKYL_BASIS_MODAL_TENSOR); }
void test_2x2v_p1_continuity_ser() { test_2x2v(1, GKYL_BASIS_MODAL_SERENDIPITY); }
void test_2x2v_p2_continuity_tensor() { test_2x2v(2, GKYL_BASIS_MODAL_TENSOR); }
void test_2x2v_p2_continuity_ser() { test_2x2v(2, GKYL_BASIS_MODAL_SERENDIPITY); }
// Check the 2x3v_p1/2 continuity
void test_2x3v_p1_continuity_tensor() { test_2x3v(1, GKYL_BASIS_MODAL_TENSOR); }
void test_2x3v_p1_continuity_ser() { test_2x3v(1, GKYL_BASIS_MODAL_SERENDIPITY); }
void test_2x3v_p2_continuity_ser() { test_2x3v(2, GKYL_BASIS_MODAL_SERENDIPITY); }
// Check the 3x3v p1 (tensor only)
void test_3x3v_p1_continuity_tensor() { test_3x3v(1, GKYL_BASIS_MODAL_TENSOR); }

TEST_LIST = {
  {"test_1x1v_p1_continuity_ten", test_1x1v_p1_continuity_tensor},
  {"test_1x1v_p1_continuity_hyb", test_1x1v_p1_continuity_ser},
  {"test_1x1v_p2_continuity_ten", test_1x1v_p2_continuity_tensor},
  {"test_1x1v_p2_continuity_ser", test_1x1v_p2_continuity_ser},
  {"test_1x2v_p1_continuity_ten", test_1x2v_p1_continuity_tensor},
  {"test_1x2v_p1_continuity_hyb", test_1x2v_p1_continuity_ser},
  {"test_1x2v_p2_continuity_ten", test_1x2v_p2_continuity_tensor},
  {"test_1x2v_p2_continuity_ser", test_1x2v_p2_continuity_ser},
  {"test_1x3v_p1_continuity_ten", test_1x3v_p1_continuity_tensor},
  {"test_1x3v_p1_continuity_hyb", test_1x3v_p1_continuity_ser},
  {"test_1x3v_p2_continuity_ten", test_1x3v_p2_continuity_tensor},
  {"test_1x3v_p2_continuity_ser", test_1x3v_p2_continuity_ser},
  {"test_2x2v_p1_continuity_ten", test_2x2v_p1_continuity_tensor},
  {"test_2x2v_p1_continuity_hyb", test_2x2v_p1_continuity_ser},
  {"test_2x2v_p2_continuity_ten", test_2x2v_p2_continuity_tensor},
  {"test_2x2v_p2_continuity_ser", test_2x2v_p2_continuity_ser},
  {"test_2x3v_p1_continuity_ten", test_2x3v_p1_continuity_tensor},
  {"test_2x3v_p1_continuity_hyb", test_2x3v_p1_continuity_ser},
  {"test_2x3v_p2_continuity_ser", test_2x3v_p2_continuity_ser},
  {"test_3x3v_p1_continuity_ten", test_3x3v_p1_continuity_tensor},
  {NULL, NULL},
};
