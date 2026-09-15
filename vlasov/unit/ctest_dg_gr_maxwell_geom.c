#include <acutest.h>

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_bc_basic.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#include <gkyl_dg_gr_maxwell_geom.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_dg_gr_maxwell_kernels.h>
#include <gkyl_dg_gr_maxwell_divide_Jc.h>
#include <gkyl_dg_gr_maxwell_divide_Jc_priv.h>
#include <gkyl_vlasov_position_map.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

static inline void
log_to_comp(int ndim, const double *eta, const double *dx, const double *xc, double *xout)
{
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d] + xc[d];
}

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr, const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;
  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d], d,
      GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d], d,
      GKYL_UPPER_EDGE, parent, ghost);
  }
}

static void
apply_expected_theta_pole_bc_cell(const struct gkyl_basis *basis, const double *skin, double *ghost)
{
  int dir = 1, nbasis = basis->num_basis;

  basis->flip_odd_sign(dir,  &skin[nbasis*0], &ghost[nbasis*0]);
  basis->flip_even_sign(dir, &skin[nbasis*1], &ghost[nbasis*1]);
  basis->flip_odd_sign(dir,  &skin[nbasis*2], &ghost[nbasis*2]);
  basis->flip_odd_sign(dir,  &skin[nbasis*3], &ghost[nbasis*3]);
  basis->flip_even_sign(dir, &skin[nbasis*4], &ghost[nbasis*4]);
  basis->flip_odd_sign(dir,  &skin[nbasis*5], &ghost[nbasis*5]);
}

void field_init_func(double t, const double *xn, double* restrict fout, void *ctx)
{

  struct gkyl_dg_gr_maxwell_geom_ctx *app = ctx;
  double r = xn[0]; 
  double theta = xn[1];
  double massBH = app->mass_bh;
  double spinBH = app->spin_bh;
  double B_0 = 1.0; // Magnitude of the fields

  // Fields
  double Dr = 0.0; // Total electric field (r-direction).
  double Dtheta = 0.0; // Total electric field (theta-direction).
  double Dphi = - 2.0 * massBH * B_0 / ( r * r * sqrt( 1.0 + 2.0 * massBH / r ) );  // Total electric field (phi-direction).

  double Br = - B_0 * cos(theta) / ( sqrt( 1.0 + 2.0 * massBH / r ) );   // Total magnetic field (r-direction).
  double Btheta = B_0 * sin(theta) / ( r * sqrt( 1.0 + 2.0 * massBH / r ) ); // Total magnetic field (theta-direction).
  double Bphi = 0.0; // Total magnetic field (phi-direction).

  // Must return conserved variables
  double rho = sqrt(r * r + spinBH * spinBH * cos(theta) * cos(theta) );
  double metric_det = rho * sqrt(2*massBH*r + rho * rho) * sin(theta);
  
  // Compute Jc * D^i and Jc * B^i
  double JDr = metric_det * Dr;
  double JDtheta = metric_det * Dtheta;
  double JDphi = metric_det * Dphi;
  double JBr = metric_det * Br;
  double JBtheta = metric_det * Btheta;
  double JBphi = metric_det * Bphi;

  // Hand off the conserved varaibles (J * Q^\xi)
  fout[0] = JDr;
  fout[1] = JDtheta;
  fout[2] = JDphi;
  fout[3] = JBr;
  fout[4] = JBtheta;
  fout[5] = JBphi;
  fout[6] = 0.0;
  fout[7] = 0.0;
}

static const double Dphi_ref[4][4] = {
  {-1.6354821637490904e+00, +3.1682899412226662e-01, +0.0000000000000000e+00, +0.0000000000000000e+00},
  {-1.6354821637490904e+00, +3.1682899412226662e-01, +0.0000000000000000e+00, +0.0000000000000000e+00},
  {-1.6354821637490904e+00, +3.1682899412226662e-01, +0.0000000000000000e+00, +0.0000000000000000e+00},
  {-1.6354821637490904e+00, +3.1682899412226662e-01, +0.0000000000000000e+00, +0.0000000000000000e+00}
};

static const double Br_ref[4][4] = {
  {-1.1141258010329544e+00, -3.9838967223318700e-02, +1.0646084038672876e-01, +3.8068321609655564e-03},
  {-4.6148601697763836e-01, -1.6501840534835825e-02, +2.5701920472327799e-01, +9.1905058326813213e-03},
  {+4.6148601697763814e-01, +1.6501840534835860e-02, +2.5701920472327799e-01, +9.1905058326813283e-03},
  {+1.1141258010329544e+00, +3.9838967223318700e-02, +1.0646084038672876e-01, +3.8068321609655564e-03}
};

static const double Btheta_ref[4][4] = {
  {+3.7263287721296423e-01, -2.9826465965693272e-02, +2.0753349447565991e-01, -1.6611499114132122e-02},
  {+8.9961534595364645e-01, -7.2007458652036282e-02, +8.5963188058500117e-02, -6.8807082244221884e-03},
  {+8.9961534595364645e-01, -7.2007458652036255e-02, -8.5963188058500034e-02, +6.8807082244221884e-03},
  {+3.7263287721296434e-01, -2.9826465965693265e-02, -2.0753349447565989e-01, +1.6611499114132122e-02}
};

static const double JDphi_ref[4][4] = {
  {-1.4915589107194009e+00, +0.0000000000000000e+00, -8.3070617727859664e-01, +0.0000000000000000e+00},
  {-3.6009417513372184e+00, -3.3306690738754696e-16, -3.4408976497590316e-01, +0.0000000000000000e+00},
  {-3.6009417513372184e+00, -2.2204460492503131e-16, +3.4408976497590282e-01, +0.0000000000000000e+00},
  {-1.4915589107194016e+00, +0.0000000000000000e+00, +8.3070617727859641e-01, -2.7755575615628914e-17}
};

static const double JBr_ref[4][4] = {
  {-1.0064417291126031e+00, -2.2936950127375039e-01, -4.9045618369336014e-01, -1.1177566171621522e-01},
  {-1.0064417291126033e+00, -2.2936950127375047e-01, +4.9045618369336008e-01, +1.1177566171621517e-01},
  {+1.0064417291126031e+00, +2.2936950127375044e-01, +4.9045618369336008e-01, +1.1177566171621517e-01},
  {+1.0064417291126035e+00, +2.2936950127375053e-01, -4.9045618369335997e-01, -1.1177566171621522e-01}
};

static const double JBtheta_ref[4][4] = {
  {+4.5544074017426023e-01, +5.2589766787906389e-02, +3.8720225028423161e-01, +4.4710264686486062e-02},
  {+2.0445592598257387e+00, +2.3608536780690687e-01, +3.8720225028423155e-01, +4.4710264686485979e-02},
  {+2.0445592598257387e+00, +2.3608536780690687e-01, -3.8720225028423133e-01, -4.4710264686486034e-02},
  {+4.5544074017426045e-01, +5.2589766787906451e-02, -3.8720225028423166e-01, -4.4710264686486104e-02}
};

static const char *em_comp_names[8] = {
  "Dr", "Dtheta", "Dphi", "Br", "Btheta", "Bphi", "phi", "psi"
};

void
test_ks_r_theta_2x_geom_p1()
{


  // GR paramters
  // populate the geometry context struct with mass and spin
  struct gkyl_dg_gr_maxwell_geom_ctx ctx = {
    .mass_bh = 1.0,
    .spin_bh = 0.0
  };
  double Mass = ctx.mass_bh;
  double Spin = ctx.spin_bh;

  // poylnomial order for this test is fixed to 1
  int poly_order = 1;

  double lower[] = {1.0, 0.0}, upper[] = {1.5, M_PI};
  int cells[] = {1, 4};
  int cdim = 2;

  double confLower[] = {lower[0], lower[1]}, confUpper[] = {upper[0], upper[1]};
  int confCells[] = {cells[0], cells[1]};

  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis confBasis;
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1, 1 };
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  // Uniform grid: default (identity) position map, so J_pos = 1 in every cell.
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = { 0 };
  struct gkyl_vlasov_position_map *pos_map = gkyl_vlasov_position_map_new(&confGrid,
    &confLocal, &confLocal_ext, &confBasis, inp_pmap, false);

  // Create the fields (JD^i, JB^i) by setting the bare minimum required fields
  // in gkyl_vlasov_app, gkyl_vm input. 
  
  // Create the geometry alpha, \beta^i, h_ij, and J_c nodally 

  // Record which configured field boundaries use the theta-pole BC.
  enum gkyl_triad_preset_geom_type triad_preset_geom_type = GKYL_TRIAD_GR_KERR_SCHILD_RTHETA;
  struct gkyl_surf_and_vol_node_arrays *lapse, *shift, *h_ij, *h_ij_inv, *det_h;

  // Evaluation of geometry at surface and volume nodal points.
  // Lapse - \alpha in the ADM split
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* lapse_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(
    &confGrid, &confBasis, 1, poly_order, gkyl_dg_gr_maxwell_preset_lapse(triad_preset_geom_type), &ctx);
  lapse = gkyl_surf_and_vol_node_arrays_new(lapse_proj, confLocal_ext.volume, false);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(lapse_proj, 0.0, &confLocal_ext, lapse);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(lapse_proj);

  // shift - \beta^i components in the ADM split (contravariant)
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* shift_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(
    &confGrid, &confBasis, 3, poly_order, gkyl_dg_gr_maxwell_preset_shift(triad_preset_geom_type), &ctx);
  shift = gkyl_surf_and_vol_node_arrays_new(shift_proj, confLocal_ext.volume, false);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(shift_proj, 0.0, &confLocal_ext, shift);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(shift_proj);

  // h_ij - Covariant components of the spatial metric (assumed to allways be a upper 
  // triangular matrix of 6 unique elements)
  // Allocate arrays for specified metric inverse
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* h_ij_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(
    &confGrid, &confBasis, 6, poly_order, gkyl_dg_gr_maxwell_preset_h_ij(triad_preset_geom_type), &ctx);
  h_ij = gkyl_surf_and_vol_node_arrays_new(h_ij_proj, confLocal_ext.volume, false);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(h_ij_proj, 0.0, &confLocal_ext, h_ij);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(h_ij_proj);

  // h_ij_inv - Contravariant components of the spatial metric inverse (assumed to allways be a upper
  // triangular matrix of 6 unique elements)
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* h_ij_inv_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(
    &confGrid, &confBasis, 6, poly_order, gkyl_dg_gr_maxwell_preset_h_ij_inv(triad_preset_geom_type), &ctx);
  h_ij_inv = gkyl_surf_and_vol_node_arrays_new(h_ij_inv_proj, confLocal_ext.volume, false);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(h_ij_inv_proj, 0.0, &confLocal_ext, h_ij_inv);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(h_ij_inv_proj);

  // Allocate arrays for the metric determinant (computed from J = sqrt(det(h_ij)))
  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes* det_h_proj = gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(
    &confGrid, &confBasis, 1, poly_order, gkyl_dg_gr_maxwell_preset_det_h(triad_preset_geom_type), &ctx);
  det_h = gkyl_surf_and_vol_node_arrays_new(det_h_proj, confLocal_ext.volume, false);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(det_h_proj, 0.0, &confLocal_ext, det_h);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(det_h_proj);

  // Number of nodes (surf and volume)
  int num_surf_nodes = (poly_order + 1);
  int num_vol_nodes = (poly_order + 1) * (poly_order + 1);
  const double *gl_nodes = gkyl_gauss_ordinates[poly_order + 1];


  // (TEST I) Geometry objects projected nodally 
  // Check the geometry has been correctly built against analytic expressions
  // Iterate over the grid, conf space, checking output
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    double xc[GKYL_MAX_DIM], xphys[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&confGrid, iter.idx, xc);

    // 1. Compare the surface-x alpha 
    const double *lapse_surf_x_d = gkyl_array_cfetch(lapse->nodal_arr_surf_x, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_surf_nodes; ++m) {
      double expected = 0;

      double eta[GKYL_MAX_DIM] = { -1.0, gl_nodes[m] };
      log_to_comp(cdim, eta, confGrid.dx, xc, xphys);

      // For each index, grab the expected value, for surf_x alpha, they happen to all be the same as 
      // alpha(r) = 1/sqrt( 1 + 2M/r ). For Schwarzschild
      double alpha_lower_x_surf = 1.0 / sqrt( 1.0 + 2.0 * Mass / xphys[0]);
      // if (iter.idx[1] == 1) expected = alpha_lower_x_surf;
      // if (iter.idx[1] == 2) expected = alpha_lower_x_surf;
      // if (iter.idx[1] == 3) expected = alpha_lower_x_surf;
      // if (iter.idx[1] == 4) expected = alpha_lower_x_surf;
      expected = alpha_lower_x_surf;
      //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,lapse_surf_x_d[m], expected);
      TEST_CHECK( gkyl_compare_double(lapse_surf_x_d[m], expected, 1e-12) );
    }

    // 2. Compare the surface-y alpha
    const double *lapse_surf_y_d = gkyl_array_cfetch(lapse->nodal_arr_surf_y, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_surf_nodes; ++m) {
      double expected = 0;

      // Compute the radius we are evaluating at
      double eta[GKYL_MAX_DIM] = { gl_nodes[m], -1.0 };
      log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
      double r_val = xphys[0];

      // For each index, grab the expected value, for surf_y alpha, they happen to all be the same as 
      // alpha(r) = 1/sqrt( 1 + 2M/r ). For Schwarzschild
      double alpha_lower_y_surf = 1.0 / sqrt( 1.0 + 2.0 * Mass / r_val);
      expected = alpha_lower_y_surf;
      //printf("(r=%1.1f), Value (Surf-y alpha) [%d]: %1.16e, Reference: %1.16e\n",r_val, m,lapse_surf_y_d[m], expected);
      TEST_CHECK( gkyl_compare_double(lapse_surf_y_d[m], expected, 1e-12) );
    }

    // 3. Compare the volume alpha
    const double *lapse_vol_d = gkyl_array_cfetch(lapse->nodal_arr_vol, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_vol_nodes; ++m) {
      double expected = 0;

      // Compute the radius we are evaluating at
      int theta_indx = m % num_surf_nodes;
      int r_indx = m / num_surf_nodes;
      double eta[GKYL_MAX_DIM] = { gl_nodes[r_indx], gl_nodes[theta_indx] };
      log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
      double r_val = xphys[0];

      // For each index, grab the expected value, for vol alpha, they happen to all be the same as 
      // alpha(r) = 1/sqrt( 1 + 2M/r ). For Schwarzschild
      double alpha_vol = 1.0 / sqrt( 1.0 + 2.0 * Mass / r_val);
      expected = alpha_vol;
      //printf("(r=%1.1f), Value (vol alpha) [%d]: %1.16e, Reference: %1.16e\n",r_val, m,lapse_vol_d[m], expected);
      TEST_CHECK( gkyl_compare_double(lapse_vol_d[m], expected, 1e-12) );
    }




    // 4. Compare the surface-x beta^i 
    const double *shift_surf_x_d = gkyl_array_cfetch(shift->nodal_arr_surf_x, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_surf_nodes; ++m) {

      // Shift has 3 components
      double NC = 3; 
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        double eta[GKYL_MAX_DIM] = { -1.0, gl_nodes[m] };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];

        // For each index, grab the expected value, for surf_x beta, they happen to all be the same as 
        // beta(r) = 1/sqrt( 1 + 2M/r ). For Schwarzschild
        double beta_lower_x_surf;
        if (n == 0) beta_lower_x_surf = ( 2.0 * Mass * r_val ) / ( r_val * r_val + 2.0 * Mass * r_val );
        if (n == 1) beta_lower_x_surf = 0.0;
        if (n == 2) beta_lower_x_surf = 0.0;
        expected = beta_lower_x_surf;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,shift_surf_x_d[m + n*num_surf_nodes], expected);
        TEST_CHECK( gkyl_compare_double(shift_surf_x_d[m + n*num_surf_nodes], expected, 1e-12) );
      }
    }


    // 5. Compare the surface-y beta^i 
    const double *shift_surf_y_d = gkyl_array_cfetch(shift->nodal_arr_surf_y, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_surf_nodes; ++m) {

      // Shift has 3 components
      double NC = 3; 
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        // Compute the radius we are evaluating at
        double eta[GKYL_MAX_DIM] = { gl_nodes[m], -1.0 };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];

        // For each index, grab the expected value, for surf_y beta, they happen to all be the same as 
        // beta(r) = 1/sqrt( 1 + 2M/r ). For Schwarzschild
        double beta_lower_y_surf;
        if (n == 0) beta_lower_y_surf = ( 2.0 * Mass * r_val ) / ( r_val * r_val + 2.0 * Mass * r_val );
        if (n == 1) beta_lower_y_surf = 0.0;
        if (n == 2) beta_lower_y_surf = 0.0;
        expected = beta_lower_y_surf;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,shift_surf_y_d[m + n*num_surf_nodes], expected);
        TEST_CHECK( gkyl_compare_double(shift_surf_y_d[m + n*num_surf_nodes], expected, 1e-12) );
      }
    }

    // 6. Compare the vol beta^i 
    const double *shift_vol_d = gkyl_array_cfetch(shift->nodal_arr_vol, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_vol_nodes; ++m) {

      // Shift has 3 components
      double NC = 3; 
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        // Compute the radius we are evaluating at
        int theta_indx = m % num_surf_nodes;
        int r_indx = m / num_surf_nodes;
        double eta[GKYL_MAX_DIM] = { gl_nodes[r_indx], gl_nodes[theta_indx] };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];

        // For each index, grab the expected value, for vol beta, they happen to all be the same as 
        // beta(r) = 1/sqrt( 1 + 2M/r ). For Schwarzschild
        double beta_vol;
        if (n == 0) beta_vol = ( 2.0 * Mass * r_val ) / ( r_val * r_val + 2.0 * Mass * r_val );
        if (n == 1) beta_vol = 0.0;
        if (n == 2) beta_vol = 0.0;
        expected = beta_vol;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,shift_vol_d[m + n*num_vol_nodes], expected);
        TEST_CHECK( gkyl_compare_double(shift_vol_d[m + n*num_vol_nodes], expected, 1e-12) );
      }
    }



    // 7. Compare the surface-x h_ij
    const double *h_ij_surf_x_d = gkyl_array_cfetch(h_ij->nodal_arr_surf_x, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_surf_nodes; ++m) {

      // h_ij has 3 components
      double NC = 6; 
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        double eta[GKYL_MAX_DIM] = { -1.0, gl_nodes[m] };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];
        double theta_val = xphys[1];

        // For each index, grab the expected value, for surf_x h_ij
        double h_ij_lower_x_surf;
        if (n == 0) h_ij_lower_x_surf = (1.0 + 2.0 * Mass / r_val);
        if (n == 1) h_ij_lower_x_surf = 0.0;
        if (n == 2) h_ij_lower_x_surf = 0.0;
        if (n == 3) h_ij_lower_x_surf = r_val * r_val;
        if (n == 4) h_ij_lower_x_surf = 0.0;
        if (n == 5) h_ij_lower_x_surf = r_val * r_val * sin(theta_val) * sin(theta_val);
        expected = h_ij_lower_x_surf;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,h_ij_surf_x_d[m + n*num_surf_nodes], expected);
        TEST_CHECK( gkyl_compare_double(h_ij_surf_x_d[m + n*num_surf_nodes], expected, 1e-12) );
      }
    }


    // 8. Compare the surface-y h_ij 
    const double *h_ij_surf_y_d = gkyl_array_cfetch(h_ij->nodal_arr_surf_y, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_surf_nodes; ++m) {

      // h_ij has 3 components
      double NC = 6; 
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        // Compute the radius we are evaluating at
        double eta[GKYL_MAX_DIM] = { gl_nodes[m], -1.0 };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];
        double theta_val = xphys[1];

        // For each index, grab the expected value, for surf_y h_ij
        double h_ij_lower_y_surf;
        if (n == 0) h_ij_lower_y_surf = (1.0 + 2.0 * Mass / r_val);
        if (n == 1) h_ij_lower_y_surf = 0.0;
        if (n == 2) h_ij_lower_y_surf = 0.0;
        if (n == 3) h_ij_lower_y_surf = r_val * r_val;
        if (n == 4) h_ij_lower_y_surf = 0.0;
        if (n == 5) h_ij_lower_y_surf = r_val * r_val * sin(theta_val) * sin(theta_val);
        expected = h_ij_lower_y_surf;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,h_ij_surf_y_d[m + n*num_surf_nodes], expected);
        TEST_CHECK( gkyl_compare_double(h_ij_surf_y_d[m + n*num_surf_nodes], expected, 1e-12) );
      }
    }

    // 9. Compare the vol h_ij
    const double *h_ij_vol_d = gkyl_array_cfetch(h_ij->nodal_arr_vol, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_vol_nodes; ++m) {

      // h_ij has 3 components
      double NC = 6; 
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        // Compute the radius we are evaluating at
        int theta_indx = m % num_surf_nodes;
        int r_indx = m / num_surf_nodes;
        double eta[GKYL_MAX_DIM] = { gl_nodes[r_indx], gl_nodes[theta_indx] };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];
        double theta_val = xphys[1];

        // For each index, grab the expected value, for vol h_ij
        double h_ij_vol;
        if (n == 0) h_ij_vol = (1.0 + 2.0 * Mass / r_val);
        if (n == 1) h_ij_vol = 0.0;
        if (n == 2) h_ij_vol = 0.0;
        if (n == 3) h_ij_vol = r_val * r_val;
        if (n == 4) h_ij_vol = 0.0;
        if (n == 5) h_ij_vol = r_val * r_val * sin(theta_val) * sin(theta_val);
        expected = h_ij_vol;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,h_ij_vol_d[m + n*num_vol_nodes], expected);
        TEST_CHECK( gkyl_compare_double(h_ij_vol_d[m + n*num_vol_nodes], expected, 1e-12) );
      }
    }


    // 10. Compare the surface-x h_ij_inv
    const double *h_ij_inv_surf_x_d = gkyl_array_cfetch(h_ij_inv->nodal_arr_surf_x, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_surf_nodes; ++m) {

      // h_ij_inv has 6 components
      double NC = 6;
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        double eta[GKYL_MAX_DIM] = { -1.0, gl_nodes[m] };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];
        double theta_val = xphys[1];

        // For each index, grab the expected value, for surf_x h_ij_inv
        double h_ij_inv_lower_x_surf;
        if (n == 0) h_ij_inv_lower_x_surf = 1.0/(1.0 + 2.0 * Mass / r_val);
        if (n == 1) h_ij_inv_lower_x_surf = 0.0;
        if (n == 2) h_ij_inv_lower_x_surf = 0.0;
        if (n == 3) h_ij_inv_lower_x_surf = 1.0/(r_val * r_val);
        if (n == 4) h_ij_inv_lower_x_surf = 0.0;
        if (n == 5) h_ij_inv_lower_x_surf = 1.0/(r_val * r_val * sin(theta_val) * sin(theta_val));
        expected = h_ij_inv_lower_x_surf;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,h_ij_inv_surf_x_d[m + n*num_surf_nodes], expected);
        TEST_CHECK( gkyl_compare_double(h_ij_inv_surf_x_d[m + n*num_surf_nodes], expected, 1e-12) );
      }
    }

    // 11. Compare the surface-y h_ij_inv
    const double *h_ij_inv_surf_y_d = gkyl_array_cfetch(h_ij_inv->nodal_arr_surf_y, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_surf_nodes; ++m) {

      // h_ij_inv has 6 components
      double NC = 6;
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        // Compute the radius we are evaluating at
        double eta[GKYL_MAX_DIM] = { gl_nodes[m], -1.0 };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];
        double theta_val = xphys[1];

        if (fabs(sin(theta_val)) < 1e-14) {
          continue;
        }

        // For each index, grab the expected value, for surf_y h_ij_inv
        double h_ij_inv_lower_y_surf;
        if (n == 0) h_ij_inv_lower_y_surf = 1.0/(1.0 + 2.0 * Mass / r_val);
        if (n == 1) h_ij_inv_lower_y_surf = 0.0;
        if (n == 2) h_ij_inv_lower_y_surf = 0.0;
        if (n == 3) h_ij_inv_lower_y_surf = 1.0/(r_val * r_val);
        if (n == 4) h_ij_inv_lower_y_surf = 0.0;
        if (n == 5) h_ij_inv_lower_y_surf = 1.0/(r_val * r_val * sin(theta_val) * sin(theta_val));
        expected = h_ij_inv_lower_y_surf;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,h_ij_inv_surf_y_d[m + n*num_surf_nodes], expected);
        TEST_CHECK( gkyl_compare_double(h_ij_inv_surf_y_d[m + n*num_surf_nodes], expected, 1e-12) );
      }
    }

    // 12. Compare the vol h_ij_inv
    const double *h_ij_inv_vol_d = gkyl_array_cfetch(h_ij_inv->nodal_arr_vol, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_vol_nodes; ++m) {

      // h_ij_inv has 6 components
      double NC = 6;
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        // Compute the radius we are evaluating at
        int theta_indx = m % num_surf_nodes;
        int r_indx = m / num_surf_nodes;
        double eta[GKYL_MAX_DIM] = { gl_nodes[r_indx], gl_nodes[theta_indx] };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];
        double theta_val = xphys[1];

        // For each index, grab the expected value, for vol h_ij_inv
        double h_ij_inv_vol;
        if (n == 0) h_ij_inv_vol = 1.0/(1.0 + 2.0 * Mass / r_val);
        if (n == 1) h_ij_inv_vol = 0.0;
        if (n == 2) h_ij_inv_vol = 0.0;
        if (n == 3) h_ij_inv_vol = 1.0/(r_val * r_val);
        if (n == 4) h_ij_inv_vol = 0.0;
        if (n == 5) h_ij_inv_vol = 1.0/(r_val * r_val * sin(theta_val) * sin(theta_val));
        expected = h_ij_inv_vol;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,h_ij_inv_vol_d[m + n*num_vol_nodes], expected);
        TEST_CHECK( gkyl_compare_double(h_ij_inv_vol_d[m + n*num_vol_nodes], expected, 1e-12) );
      }
    }

    // 13. Compare the surface-x det_h
    const double *det_h_surf_x_d = gkyl_array_cfetch(det_h->nodal_arr_surf_x, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_surf_nodes; ++m) {

      // det_h has 1 component
      double NC = 1; 
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        double eta[GKYL_MAX_DIM] = { -1.0, gl_nodes[m] };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];
        double theta_val = xphys[1];

        // For each index, grab the expected value, for surf_x det_h
        double det_h_lower_x_surf = r_val * r_val * sin(theta_val) * sqrt( 1.0 + 2.0 * Mass / r_val);
        expected = det_h_lower_x_surf;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,det_h_surf_x_d[m + n*num_surf_nodes], expected);
        TEST_CHECK( gkyl_compare_double(det_h_surf_x_d[m + n*num_surf_nodes], expected, 1e-12) );
      }
    }


    // 14. Compare the surface-y det_h
    const double *det_h_surf_y_d = gkyl_array_cfetch(det_h->nodal_arr_surf_y, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_surf_nodes; ++m) {

      // det_h has 1 component
      double NC = 1; 
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        // Compute the radius we are evaluating at
        double eta[GKYL_MAX_DIM] = { gl_nodes[m], -1.0 };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];
        double theta_val = xphys[1];

        // For each index, grab the expected value, for surf_y det_h
        double det_h_lower_y_surf = r_val * r_val * sin(theta_val) * sqrt( 1.0 + 2.0 * Mass / r_val);
        expected = det_h_lower_y_surf;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,det_h_surf_y_d[m + n*num_surf_nodes], expected);
        TEST_CHECK( gkyl_compare_double(det_h_surf_y_d[m + n*num_surf_nodes], expected, 1e-12) );
      }
    }

    // 15. Compare the vol det_h
    const double *det_h_vol_d = gkyl_array_cfetch(det_h->nodal_arr_vol, gkyl_range_idx(&confLocal, iter.idx));
    for (int m = 0; m<num_vol_nodes; ++m) {

      // det_h has 1 component
      double NC = 1; 
      for (int n = 0; n<NC; ++n) {
        double expected = 0;

        // Compute the radius we are evaluating at
        int theta_indx = m % num_surf_nodes;
        int r_indx = m / num_surf_nodes;
        double eta[GKYL_MAX_DIM] = { gl_nodes[r_indx], gl_nodes[theta_indx] };
        log_to_comp(cdim, eta, confGrid.dx, xc, xphys);
        double r_val = xphys[0];
        double theta_val = xphys[1];

        // For each index, grab the expected value, for vol det_h
        double det_h_vol = r_val * r_val * sin(theta_val) * sqrt( 1.0 + 2.0 * Mass / r_val);
        expected = det_h_vol;
        //printf("Value[%d]: %1.16e, Reference: %1.16e\n",m,det_h_vol_d[m + n*num_vol_nodes], expected);
        TEST_CHECK( gkyl_compare_double(det_h_vol_d[m + n*num_vol_nodes], expected, 1e-12) );
      }
    }
  }

  // (TEST II) division by J_c
  // Allocate the field arrays in memory
  struct gkyl_array *field_no_J_con, *field_with_J_con;
  field_no_J_con = mkarr(8*confBasis.num_basis, confLocal_ext.volume);
  field_with_J_con = mkarr(8*confBasis.num_basis, confLocal_ext.volume);

  // Project the inital conditions (Kerr-Schild Coordinates, Schwarzschild Wald solution)
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1,
     8, field_init_func, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &confLocal_ext, field_with_J_con);
  gkyl_proj_on_basis_release(proj);

  // Divide out the spatial Jacobian
  gkyl_dg_gr_maxwell_divide_Jc(&confBasis, &confLocal, det_h,
    field_with_J_con, field_no_J_con, false);

  // Test the fields before and after division by J_c (mode by mode).
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&confLocal, iter.idx);

    const double *em      = gkyl_array_cfetch(field_with_J_con, lidx);
    const double *em_no_J = gkyl_array_cfetch(field_no_J_con, lidx);

    // with J
    const double *JDr     = &em[0*confBasis.num_basis];
    const double *JDtheta = &em[1*confBasis.num_basis];
    const double *JDphi   = &em[2*confBasis.num_basis];
    const double *JBr     = &em[3*confBasis.num_basis];
    const double *JBtheta = &em[4*confBasis.num_basis];
    const double *JBphi   = &em[5*confBasis.num_basis];

    // without J
    const double *Dr      = &em_no_J[0*confBasis.num_basis];
    const double *Dtheta  = &em_no_J[1*confBasis.num_basis];
    const double *Dphi    = &em_no_J[2*confBasis.num_basis];
    const double *Br      = &em_no_J[3*confBasis.num_basis];
    const double *Btheta  = &em_no_J[4*confBasis.num_basis];
    const double *Bphi    = &em_no_J[5*confBasis.num_basis];

    int jloc = iter.idx[1] - confLocal.lower[1]; // 0..3 for Nz1 = 4

    const double *Dphi_mode_ref    = Dphi_ref[jloc];
    const double *Br_mode_ref      = Br_ref[jloc];
    const double *Btheta_mode_ref  = Btheta_ref[jloc];

    const double *JDphi_mode_ref   = JDphi_ref[jloc];
    const double *JBr_mode_ref     = JBr_ref[jloc];
    const double *JBtheta_mode_ref = JBtheta_ref[jloc];

    // Loop over all bases
    for (int k = 0; k < confBasis.num_basis; ++k) {

      // printf("\n---------\n");    
      // printf("(Dphi) Cell[%ld](basis: %d): %1.16e, Reference: %1.16e\n",lidx, k,Dphi[k],   Dphi_mode_ref[k]);
      // printf("(Br)   Cell[%ld](basis: %d): %1.16e, Reference: %1.16e\n",lidx, k,Br[k],     Br_mode_ref[k]);
      // printf("(Bphi) Cell[%ld](basis: %d): %1.16e, Reference: %1.16e\n",lidx, k,Btheta[k], Btheta_mode_ref[k]);  
      // printf("(JDphi) Cell[%ld](basis: %d): %1.16e, Reference: %1.16e\n",lidx,k,JDphi[k],   JDphi_mode_ref[k]);
      // printf("(JBr)   Cell[%ld](basis: %d): %1.16e, Reference: %1.16e\n",lidx,k,JBr[k],     JBr_mode_ref[k]);
      // printf("(JBphi) Cell[%ld](basis: %d): %1.16e, Reference: %1.16e\n",lidx,k,JBtheta[k], JBtheta_mode_ref[k]);

      // no J
      TEST_CHECK(gkyl_compare_double(Dr[k],     0.0,                1e-12));
      TEST_CHECK(gkyl_compare_double(Dtheta[k], 0.0,                1e-12));
      TEST_CHECK(gkyl_compare_double(Dphi[k],   Dphi_mode_ref[k],   1e-12));
      TEST_CHECK(gkyl_compare_double(Br[k],     Br_mode_ref[k],     1e-12));
      TEST_CHECK(gkyl_compare_double(Btheta[k], Btheta_mode_ref[k], 1e-12));
      TEST_CHECK(gkyl_compare_double(Bphi[k],   0.0,                1e-12));

      // with J
      TEST_CHECK(gkyl_compare_double(JDr[k],     0.0,                 1e-12));
      TEST_CHECK(gkyl_compare_double(JDtheta[k], 0.0,                 1e-12));
      TEST_CHECK(gkyl_compare_double(JDphi[k],   JDphi_mode_ref[k],   1e-12));
      TEST_CHECK(gkyl_compare_double(JBr[k],     JBr_mode_ref[k],     1e-12));
      TEST_CHECK(gkyl_compare_double(JBtheta[k], JBtheta_mode_ref[k], 1e-12));
      TEST_CHECK(gkyl_compare_double(JBphi[k],   0.0,                 1e-12));
    }
  }

  
  // (TEST III) Boundary Conditions
  struct skin_ghost_ranges conf_skin_ghost;
  skin_ghost_ranges_init(&conf_skin_ghost, &confLocal_ext, confGhost);

  long bc_buff_sz = 0;
  for (int d=0; d<cdim; ++d) {
    if (conf_skin_ghost.lower_skin[d].volume > bc_buff_sz) bc_buff_sz = conf_skin_ghost.lower_skin[d].volume;
    if (conf_skin_ghost.upper_skin[d].volume > bc_buff_sz) bc_buff_sz = conf_skin_ghost.upper_skin[d].volume;
  }

  struct gkyl_range radial_ghost_lo, radial_ghost_up;
  int radial_ghost_lo_lower[] = { confLocal_ext.lower[0], confLocal.lower[1] };
  int radial_ghost_lo_upper[] = { confLocal.lower[0]-1, confLocal.upper[1] };
  int radial_ghost_up_lower[] = { confLocal.upper[0]+1, confLocal.lower[1] };
  int radial_ghost_up_upper[] = { confLocal_ext.upper[0], confLocal.upper[1] };
  gkyl_sub_range_init(&radial_ghost_lo, &confLocal_ext, radial_ghost_lo_lower, radial_ghost_lo_upper);
  gkyl_sub_range_init(&radial_ghost_up, &confLocal_ext, radial_ghost_up_lower, radial_ghost_up_upper);

  struct gkyl_array *field_with_J_radial_init = mkarr(8*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *field_no_J_radial_init = mkarr(8*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_clear(field_with_J_radial_init, 0.0);
  gkyl_array_clear(field_no_J_radial_init, 0.0);

  gkyl_proj_on_basis *proj_radial_ghost = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1,
    8, field_init_func, &ctx);
  gkyl_proj_on_basis_advance(proj_radial_ghost, 0.0, &radial_ghost_lo, field_with_J_radial_init);
  gkyl_proj_on_basis_advance(proj_radial_ghost, 0.0, &radial_ghost_up, field_with_J_radial_init);
  gkyl_proj_on_basis_release(proj_radial_ghost);

  gkyl_dg_gr_maxwell_divide_Jc(&confBasis, &radial_ghost_lo, det_h,
    field_with_J_radial_init, field_no_J_radial_init, false);
  gkyl_dg_gr_maxwell_divide_Jc(&confBasis, &radial_ghost_up, det_h,
    field_with_J_radial_init, field_no_J_radial_init, false);

  struct gkyl_array *field_no_J_fixed = mkarr(8*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *field_with_J_fixed = mkarr(8*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *field_no_J_fixed_copy = mkarr(8*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *bc_buffer = mkarr(field_no_J_con->ncomp, bc_buff_sz);
  struct gkyl_array *bc_buffer_lo_fixed = mkarr(field_no_J_con->ncomp, conf_skin_ghost.lower_ghost[0].volume);
  struct gkyl_array *bc_buffer_up_fixed = mkarr(field_no_J_con->ncomp, conf_skin_ghost.upper_ghost[0].volume);
  struct gkyl_array *bc_buffer_lo_fixed_with_J = mkarr(field_no_J_con->ncomp, conf_skin_ghost.lower_ghost[0].volume);
  struct gkyl_array *bc_buffer_up_fixed_with_J = mkarr(field_no_J_con->ncomp, conf_skin_ghost.upper_ghost[0].volume);
  struct gkyl_bc_basic *bc_lo[GKYL_MAX_DIM] = { 0 }, *bc_up[GKYL_MAX_DIM] = { 0 };

  // Compute the boundary updaters. Radial boundaries are fixed from the initial field,
  // while theta boundaries use the usual theta-pole condition.
  bc_lo[0] = gkyl_bc_basic_new(0, GKYL_LOWER_EDGE, GKYL_BC_FIXED_FUNC, &confBasis,
    &conf_skin_ghost.lower_skin[0], &conf_skin_ghost.lower_ghost[0], field_no_J_con->ncomp, cdim, false);
  bc_up[0] = gkyl_bc_basic_new(0, GKYL_UPPER_EDGE, GKYL_BC_FIXED_FUNC, &confBasis,
    &conf_skin_ghost.upper_skin[0], &conf_skin_ghost.upper_ghost[0], field_no_J_con->ncomp, cdim, false);
  bc_lo[1] = gkyl_bc_basic_new(1, GKYL_LOWER_EDGE, GKYL_BC_MAXWELL_THETA_POLE, &confBasis,
    &conf_skin_ghost.lower_skin[1], &conf_skin_ghost.lower_ghost[1], field_no_J_con->ncomp, cdim, false);
  bc_up[1] = gkyl_bc_basic_new(1, GKYL_UPPER_EDGE, GKYL_BC_MAXWELL_THETA_POLE, &confBasis,
    &conf_skin_ghost.upper_skin[1], &conf_skin_ghost.upper_ghost[1], field_no_J_con->ncomp, cdim, false);

  // For further analysis, copy the field to the boundary (in primative variables)
  struct gkyl_bc_basic *bc_lo_copy[GKYL_MAX_DIM] = { 0 };
  bc_lo_copy[0] = gkyl_bc_basic_new(0, GKYL_LOWER_EDGE, GKYL_BC_COPY, &confBasis,
    &conf_skin_ghost.lower_skin[0], &conf_skin_ghost.lower_ghost[0], field_no_J_con->ncomp, cdim, false);

  gkyl_array_copy(field_no_J_fixed, field_no_J_con);
  gkyl_array_copy(field_with_J_fixed, field_with_J_con);
  gkyl_array_copy(field_no_J_fixed_copy, field_no_J_con);

  // Seed the fixed-function radial buffers from the projected radial ghost data
  // in both primitive and conservative variables.
  gkyl_array_copy_to_buffer(bc_buffer_lo_fixed->data, field_no_J_radial_init, &radial_ghost_lo);
  gkyl_array_copy_to_buffer(bc_buffer_up_fixed->data, field_no_J_radial_init, &radial_ghost_up);
  gkyl_array_copy_to_buffer(bc_buffer_lo_fixed_with_J->data, field_with_J_radial_init, &radial_ghost_lo);
  gkyl_array_copy_to_buffer(bc_buffer_up_fixed_with_J->data, field_with_J_radial_init, &radial_ghost_up);

  // Update the boundary.
  gkyl_bc_basic_advance(bc_lo[0], bc_buffer_lo_fixed, field_no_J_fixed);
  gkyl_bc_basic_advance(bc_up[0], bc_buffer_up_fixed, field_no_J_fixed);
  gkyl_bc_basic_advance(bc_lo[0], bc_buffer_lo_fixed_with_J, field_with_J_fixed);
  gkyl_bc_basic_advance(bc_up[0], bc_buffer_up_fixed_with_J, field_with_J_fixed);
  gkyl_bc_basic_advance(bc_lo[1], bc_buffer, field_no_J_fixed);
  gkyl_bc_basic_advance(bc_up[1], bc_buffer, field_no_J_fixed);
  gkyl_bc_basic_advance(bc_lo_copy[0], bc_buffer, field_no_J_fixed_copy);

  // Test whether the created BC using the gkyl_bc_basic_advance() scheme is the same
  // as directly evaluating the primative variables and dividing out Jc
  struct gkyl_range *radial_ghost_ranges[] = { &radial_ghost_lo, &radial_ghost_up };
  const char *radial_ghost_labels[] = { "lower", "upper" };
  for (int rg=0; rg<2; ++rg) {
    gkyl_range_iter_init(&iter, radial_ghost_ranges[rg]);
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&confLocal_ext, iter.idx);
      const double *expected = gkyl_array_cfetch(field_no_J_radial_init, lidx);
      const double *actual = gkyl_array_cfetch(field_no_J_fixed, lidx);

      //printf("(TEST III) Radial %s fixed ghost cell (%d,%d)\n",
      //  radial_ghost_labels[rg], iter.idx[0], iter.idx[1]);

      for (int c=0; c<8; ++c) {
        for (int m=0; m<confBasis.num_basis; ++m) {
          int k = c*confBasis.num_basis + m;
          //printf("  %s mode %d: fixed_bc=% .16e, expected_init_no_J=% .16e\n",
          //  em_comp_names[c], m, actual[k], expected[k]);
          TEST_CHECK(gkyl_compare_double(actual[k], expected[k], 1e-12));
        }
      }
    }
  }

  for (int rg=0; rg<2; ++rg) {
    gkyl_range_iter_init(&iter, radial_ghost_ranges[rg]);
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&confLocal_ext, iter.idx);
      const double *expected = gkyl_array_cfetch(field_with_J_radial_init, lidx);
      const double *actual = gkyl_array_cfetch(field_with_J_fixed, lidx);

      for (int c=0; c<8; ++c) {
        for (int m=0; m<confBasis.num_basis; ++m) {
          int k = c*confBasis.num_basis + m;
          TEST_CHECK(gkyl_compare_double(actual[k], expected[k], 1e-12));
        }
      }
    }
  }

  // Theta boundary condition. Check against direct sign flips of the inner most cell.
  for (int ir=confLocal.lower[0]; ir<=confLocal.upper[0]; ++ir) {
    int idx_skin_lo[] = { ir, confLocal.lower[1] };
    int idx_ghost_lo[] = { ir, confLocal.lower[1]-1 };
    int idx_skin_up[] = { ir, confLocal.upper[1] };
    int idx_ghost_up[] = { ir, confLocal.upper[1]+1 };
    double expected_lo[8*4] = { 0.0 }, expected_up[8*4] = { 0.0 };

    apply_expected_theta_pole_bc_cell(&confBasis,
      gkyl_array_cfetch(field_no_J_fixed, gkyl_range_idx(&confLocal_ext, idx_skin_lo)), expected_lo);
    apply_expected_theta_pole_bc_cell(&confBasis,
      gkyl_array_cfetch(field_no_J_fixed, gkyl_range_idx(&confLocal_ext, idx_skin_up)), expected_up);

    const double *actual_lo = gkyl_array_cfetch(field_no_J_fixed, gkyl_range_idx(&confLocal_ext, idx_ghost_lo));
    const double *actual_up = gkyl_array_cfetch(field_no_J_fixed, gkyl_range_idx(&confLocal_ext, idx_ghost_up));

    for (int k=0; k<6*confBasis.num_basis; ++k) {
      TEST_CHECK(gkyl_compare_double(actual_lo[k], expected_lo[k], 1e-12));
      TEST_CHECK(gkyl_compare_double(actual_up[k], expected_up[k], 1e-12));
    }
  }

  gkyl_array_release(field_with_J_radial_init);
  gkyl_array_release(field_no_J_radial_init);
  gkyl_array_release(bc_buffer_lo_fixed);
  gkyl_array_release(bc_buffer_up_fixed);
  gkyl_array_release(bc_buffer_lo_fixed_with_J);
  gkyl_array_release(bc_buffer_up_fixed_with_J);
  gkyl_array_release(field_with_J_fixed);

  gkyl_bc_basic_release(bc_lo_copy[0]);
  for (int d=0; d<cdim; ++d) {
    gkyl_bc_basic_release(bc_lo[d]);
    gkyl_bc_basic_release(bc_up[d]);
  }
  gkyl_array_release(bc_buffer);



  // (TEST IV) Configuration space flux and maximum eigenvalue calculation
  // Compute the x-surface
  int theta_pole_lo[3] = {0,1,0};
  struct gkyl_array *conf_flux_surf = mkarr(cdim*8*num_surf_nodes, confLocal_ext.volume);
  struct gkyl_array *cflrate = mkarr(1, confLocal_ext.volume);
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM]; 
  gkyl_range_iter_init(&iter, &confLocal);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(cdim, iter.idx, idx);
    long cidx = gkyl_range_idx(&confLocal, idx);

    // Grab the cell center location for NC bracket calculation 
    double xcC[GKYL_MAX_DIM], xcR[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&confGrid, idx, xcC);

    const double *field_no_J_con_c = gkyl_array_cfetch(field_no_J_fixed, cidx); 
    const double *field_no_J_con_copy_c = gkyl_array_cfetch(field_no_J_fixed_copy, cidx); 
    const double *field_con_c = gkyl_array_cfetch(field_with_J_con, cidx); 
    double *cflrate_d = gkyl_array_fetch(cflrate, cidx);
    double *flux = gkyl_array_fetch(conf_flux_surf, cidx); 

    // Each cell owns *lower* fluxes in each configuration-space direction. 
    // So we need the distribution function in our current cell, and the cell
    // one index lower in each direction. If we are at the lower configuration-space
    // edge, we call ghost cells
    for (int dir = 0; dir<cdim; ++dir) {

      // Select the geometry elements based on direction
      const double *lapse_d = 0;
      const double *shift_d = 0;
      const double *h_ij_d = 0;
      const double *h_ij_inv_d = 0;
      const double *det_h_d = 0;

      if (dir == 0) {
        lapse_d = gkyl_array_cfetch(lapse->nodal_arr_surf_x, cidx);
        shift_d = gkyl_array_cfetch(shift->nodal_arr_surf_x, cidx);
        h_ij_d = gkyl_array_cfetch(h_ij->nodal_arr_surf_x, cidx);
        h_ij_inv_d = gkyl_array_cfetch(h_ij_inv->nodal_arr_surf_x, cidx);
        det_h_d = gkyl_array_cfetch(det_h->nodal_arr_surf_x, cidx);
      } else if (dir == 1) {
        lapse_d = gkyl_array_cfetch(lapse->nodal_arr_surf_y, cidx);
        shift_d = gkyl_array_cfetch(shift->nodal_arr_surf_y, cidx);
        h_ij_d = gkyl_array_cfetch(h_ij->nodal_arr_surf_y, cidx);
        h_ij_inv_d = gkyl_array_cfetch(h_ij_inv->nodal_arr_surf_y, cidx);
        det_h_d = gkyl_array_cfetch(det_h->nodal_arr_surf_y, cidx);
      } else {
        lapse_d = gkyl_array_cfetch(lapse->nodal_arr_surf_z, cidx);
        shift_d = gkyl_array_cfetch(shift->nodal_arr_surf_z, cidx);
        h_ij_d = gkyl_array_cfetch(h_ij->nodal_arr_surf_z, cidx);
        h_ij_inv_d = gkyl_array_cfetch(h_ij_inv->nodal_arr_surf_z, cidx);
        det_h_d = gkyl_array_cfetch(det_h->nodal_arr_surf_z, cidx);
      }

      // Create an index for the left cell (which may be a ghost cell) 
      gkyl_copy_int_arr(cdim, iter.idx, idx_l);
      idx_l[dir] = idx_l[dir]-1;
      long cidx_l = gkyl_range_idx(&confLocal_ext, idx_l);
      const double *jacob_pos_l = gkyl_array_cfetch(pos_map->jacob_pos, cidx_l);
      const double *jacob_pos_c = gkyl_array_cfetch(pos_map->jacob_pos, cidx);
      const double *field_no_J_con_l = gkyl_array_cfetch(field_no_J_fixed, cidx_l);
      const double *field_no_J_con_copy_l = gkyl_array_cfetch(field_no_J_fixed_copy, cidx_l);
      const double *field_con_l = gkyl_array_cfetch(field_with_J_con, cidx_l); 

      // For Points not along the domain-edge in theta, compute the left hand surface 
      // conf-flux.
      int theta_pole = 0;

      // If at the left edge in the theta dir
      if(idx[dir] == confLocal.lower[dir] && theta_pole_lo[dir]) {
        theta_pole = 1;
      }
      

      // 2 quadrature nodes on a 2x surface, 6 field components.
      double alpha_quad_x[3] = {0.0}; 
      double flux_l_x[18] = {0.0};
      double flux_r_x[18] = {0.0};
      double alpha_quad_y[3] = {0.0}; 
      double flux_l_y[18] = {0.0};
      double flux_r_y[18] = {0.0};

      gkyl_dg_gr_maxwell_inp meq_struct;
      const gkyl_dg_gr_maxwell_inp *meq = &meq_struct;
      meq_struct.chi = 1.0;
      meq_struct.gamma = 1.0;

      // Compute the fluxes in the first two directions
      if (dir == 0) {
        dg_gr_maxwell_alpha_quad_x_2x_ser_p1(meq, xcC, confGrid.dx, theta_pole, jacob_pos_l, jacob_pos_c,
          lapse_d, shift_d, h_ij_d, h_ij_inv_d, det_h_d, field_con_l, field_con_c, field_no_J_con_l, field_no_J_con_c, flux_l_x, flux_r_x, alpha_quad_x);

        double cflrate = lax_flux_x_2x_ser_p1(confGrid.dx, theta_pole, jacob_pos_l, jacob_pos_c, det_h_d, flux_l_x, flux_r_x, alpha_quad_x,
          field_con_l, field_con_c, field_no_J_con_l, field_no_J_con_c, flux);
      }
      else if (dir == 1) {
        dg_gr_maxwell_alpha_quad_y_2x_ser_p1(meq, xcC, confGrid.dx, theta_pole, jacob_pos_l, jacob_pos_c,
          lapse_d, shift_d, h_ij_d, h_ij_inv_d, det_h_d, field_con_l, field_con_c, field_no_J_con_l, field_no_J_con_c, flux_l_y, flux_r_y, alpha_quad_y);

        double cflrate = lax_flux_y_2x_ser_p1(confGrid.dx, theta_pole, jacob_pos_l, jacob_pos_c, det_h_d, flux_l_y, flux_r_y, alpha_quad_y,
          field_con_l, field_con_c, field_no_J_con_l, field_no_J_con_c, flux);
      }

      // Dir = 1 comparison, for Jumps in quantities at the theta edges 
      if (dir == 1) { 
        // Compute and check |alpha^i| J_c \Delta U ~ |alpha_i| \Delta Q
        double Dx_con_l_n[2] = {0.0};
        double Dy_con_l_n[2] = {0.0};
        double Dz_con_l_n[2] = {0.0};
        double Bx_con_l_n[2] = {0.0};
        double By_con_l_n[2] = {0.0};
        double Bz_con_l_n[2] = {0.0};

        double Dx_con_r_n[2] = {0.0};
        double Dy_con_r_n[2] = {0.0};
        double Dz_con_r_n[2] = {0.0};
        double Bx_con_r_n[2] = {0.0};
        double By_con_r_n[2] = {0.0};
        double Bz_con_r_n[2] = {0.0};

        double JDx_con_l_n[2] = {0.0};
        double JDy_con_l_n[2] = {0.0};
        double JDz_con_l_n[2] = {0.0};
        double JBx_con_l_n[2] = {0.0};
        double JBy_con_l_n[2] = {0.0};
        double JBz_con_l_n[2] = {0.0};

        double JDx_con_r_n[2] = {0.0};
        double JDy_con_r_n[2] = {0.0};
        double JDz_con_r_n[2] = {0.0};
        double JBx_con_r_n[2] = {0.0};
        double JBy_con_r_n[2] = {0.0};
        double JBz_con_r_n[2] = {0.0};

        double Dx_con_copy_l_n[2] = {0.0};
        double Dy_con_copy_l_n[2] = {0.0};
        double Dz_con_copy_l_n[2] = {0.0};
        double Bx_con_copy_l_n[2] = {0.0};
        double By_con_copy_l_n[2] = {0.0};
        double Bz_con_copy_l_n[2] = {0.0};

        double Dx_con_copy_r_n[2] = {0.0};
        double Dy_con_copy_r_n[2] = {0.0};
        double Dz_con_copy_r_n[2] = {0.0};
        double Bx_con_copy_r_n[2] = {0.0};
        double By_con_copy_r_n[2] = {0.0};
        double Bz_con_copy_r_n[2] = {0.0};

        const double *Dx_l = &field_no_J_con_l[0];
        const double *Dy_l = &field_no_J_con_l[4];
        const double *Dz_l = &field_no_J_con_l[8];
        const double *Bx_l = &field_no_J_con_l[12];
        const double *By_l = &field_no_J_con_l[16];
        const double *Bz_l = &field_no_J_con_l[20];

        Dx_con_l_n[0] = -(0.8660254037844386*Dx_l[3])+0.8660254037844386*Dx_l[2]-0.5*Dx_l[1]+0.5*Dx_l[0];
        Dy_con_l_n[0] = -(0.8660254037844386*Dy_l[3])+0.8660254037844386*Dy_l[2]-0.5*Dy_l[1]+0.5*Dy_l[0];
        Dz_con_l_n[0] = -(0.8660254037844386*Dz_l[3])+0.8660254037844386*Dz_l[2]-0.5*Dz_l[1]+0.5*Dz_l[0];
        Bx_con_l_n[0] = -(0.8660254037844386*Bx_l[3])+0.8660254037844386*Bx_l[2]-0.5*Bx_l[1]+0.5*Bx_l[0];
        By_con_l_n[0] = -(0.8660254037844386*By_l[3])+0.8660254037844386*By_l[2]-0.5*By_l[1]+0.5*By_l[0];
        Bz_con_l_n[0] = -(0.8660254037844386*Bz_l[3])+0.8660254037844386*Bz_l[2]-0.5*Bz_l[1]+0.5*Bz_l[0];

        Dx_con_l_n[1] = 0.8660254037844386*Dx_l[3]+0.8660254037844386*Dx_l[2]+0.5*Dx_l[1]+0.5*Dx_l[0];
        Dy_con_l_n[1] = 0.8660254037844386*Dy_l[3]+0.8660254037844386*Dy_l[2]+0.5*Dy_l[1]+0.5*Dy_l[0];
        Dz_con_l_n[1] = 0.8660254037844386*Dz_l[3]+0.8660254037844386*Dz_l[2]+0.5*Dz_l[1]+0.5*Dz_l[0];
        Bx_con_l_n[1] = 0.8660254037844386*Bx_l[3]+0.8660254037844386*Bx_l[2]+0.5*Bx_l[1]+0.5*Bx_l[0];
        By_con_l_n[1] = 0.8660254037844386*By_l[3]+0.8660254037844386*By_l[2]+0.5*By_l[1]+0.5*By_l[0];
        Bz_con_l_n[1] = 0.8660254037844386*Bz_l[3]+0.8660254037844386*Bz_l[2]+0.5*Bz_l[1]+0.5*Bz_l[0];

        const double *Dx_r = &field_no_J_con_c[0];
        const double *Dy_r = &field_no_J_con_c[4];
        const double *Dz_r = &field_no_J_con_c[8];
        const double *Bx_r = &field_no_J_con_c[12];
        const double *By_r = &field_no_J_con_c[16];
        const double *Bz_r = &field_no_J_con_c[20];

        Dx_con_r_n[0] = 0.8660254037844386*Dx_r[3]-0.8660254037844386*Dx_r[2]-0.5*Dx_r[1]+0.5*Dx_r[0];
        Dy_con_r_n[0] = 0.8660254037844386*Dy_r[3]-0.8660254037844386*Dy_r[2]-0.5*Dy_r[1]+0.5*Dy_r[0];
        Dz_con_r_n[0] = 0.8660254037844386*Dz_r[3]-0.8660254037844386*Dz_r[2]-0.5*Dz_r[1]+0.5*Dz_r[0];
        Bx_con_r_n[0] = 0.8660254037844386*Bx_r[3]-0.8660254037844386*Bx_r[2]-0.5*Bx_r[1]+0.5*Bx_r[0];
        By_con_r_n[0] = 0.8660254037844386*By_r[3]-0.8660254037844386*By_r[2]-0.5*By_r[1]+0.5*By_r[0];
        Bz_con_r_n[0] = 0.8660254037844386*Bz_r[3]-0.8660254037844386*Bz_r[2]-0.5*Bz_r[1]+0.5*Bz_r[0];

        Dx_con_r_n[1] = -(0.8660254037844386*Dx_r[3])-0.8660254037844386*Dx_r[2]+0.5*Dx_r[1]+0.5*Dx_r[0];
        Dy_con_r_n[1] = -(0.8660254037844386*Dy_r[3])-0.8660254037844386*Dy_r[2]+0.5*Dy_r[1]+0.5*Dy_r[0];
        Dz_con_r_n[1] = -(0.8660254037844386*Dz_r[3])-0.8660254037844386*Dz_r[2]+0.5*Dz_r[1]+0.5*Dz_r[0];
        Bx_con_r_n[1] = -(0.8660254037844386*Bx_r[3])-0.8660254037844386*Bx_r[2]+0.5*Bx_r[1]+0.5*Bx_r[0];
        By_con_r_n[1] = -(0.8660254037844386*By_r[3])-0.8660254037844386*By_r[2]+0.5*By_r[1]+0.5*By_r[0];
        Bz_con_r_n[1] = -(0.8660254037844386*Bz_r[3])-0.8660254037844386*Bz_r[2]+0.5*Bz_r[1]+0.5*Bz_r[0];

        // The theta-pole BC applies parity in the ghost cell, while the
        // flux kernel additionally enforces regularity by zeroing the
        // theta components at the face.
        if (theta_pole) {
          Dy_con_l_n[0] = 0.0;
          Dy_con_l_n[1] = 0.0;
          Dy_con_r_n[0] = 0.0;
          Dy_con_r_n[1] = 0.0;
          By_con_l_n[0] = 0.0;
          By_con_l_n[1] = 0.0;
          By_con_r_n[0] = 0.0;
          By_con_r_n[1] = 0.0;
        }

        // Acquire new pointer locations.
        Dx_l = &field_no_J_con_copy_l[0];
        Dy_l = &field_no_J_con_copy_l[4];
        Dz_l = &field_no_J_con_copy_l[8];
        Bx_l = &field_no_J_con_copy_l[12];
        By_l = &field_no_J_con_copy_l[16];
        Bz_l = &field_no_J_con_copy_l[20];

        Dx_con_copy_l_n[0] = -(0.8660254037844386*Dx_l[3])+0.8660254037844386*Dx_l[2]-0.5*Dx_l[1]+0.5*Dx_l[0];
        Dy_con_copy_l_n[0] = -(0.8660254037844386*Dy_l[3])+0.8660254037844386*Dy_l[2]-0.5*Dy_l[1]+0.5*Dy_l[0];
        Dz_con_copy_l_n[0] = -(0.8660254037844386*Dz_l[3])+0.8660254037844386*Dz_l[2]-0.5*Dz_l[1]+0.5*Dz_l[0];
        Bx_con_copy_l_n[0] = -(0.8660254037844386*Bx_l[3])+0.8660254037844386*Bx_l[2]-0.5*Bx_l[1]+0.5*Bx_l[0];
        By_con_copy_l_n[0] = -(0.8660254037844386*By_l[3])+0.8660254037844386*By_l[2]-0.5*By_l[1]+0.5*By_l[0];
        Bz_con_copy_l_n[0] = -(0.8660254037844386*Bz_l[3])+0.8660254037844386*Bz_l[2]-0.5*Bz_l[1]+0.5*Bz_l[0];

        Dx_con_copy_l_n[1] = 0.8660254037844386*Dx_l[3]+0.8660254037844386*Dx_l[2]+0.5*Dx_l[1]+0.5*Dx_l[0];
        Dy_con_copy_l_n[1] = 0.8660254037844386*Dy_l[3]+0.8660254037844386*Dy_l[2]+0.5*Dy_l[1]+0.5*Dy_l[0];
        Dz_con_copy_l_n[1] = 0.8660254037844386*Dz_l[3]+0.8660254037844386*Dz_l[2]+0.5*Dz_l[1]+0.5*Dz_l[0];
        Bx_con_copy_l_n[1] = 0.8660254037844386*Bx_l[3]+0.8660254037844386*Bx_l[2]+0.5*Bx_l[1]+0.5*Bx_l[0];
        By_con_copy_l_n[1] = 0.8660254037844386*By_l[3]+0.8660254037844386*By_l[2]+0.5*By_l[1]+0.5*By_l[0];
        Bz_con_copy_l_n[1] = 0.8660254037844386*Bz_l[3]+0.8660254037844386*Bz_l[2]+0.5*Bz_l[1]+0.5*Bz_l[0];

        // Acquire new pointer locations.
        Dx_r = &field_no_J_con_copy_c[0];
        Dy_r = &field_no_J_con_copy_c[4];
        Dz_r = &field_no_J_con_copy_c[8];
        Bx_r = &field_no_J_con_copy_c[12];
        By_r = &field_no_J_con_copy_c[16];
        Bz_r = &field_no_J_con_copy_c[20];

        Dx_con_copy_r_n[0] = 0.8660254037844386*Dx_r[3]-0.8660254037844386*Dx_r[2]-0.5*Dx_r[1]+0.5*Dx_r[0];
        Dy_con_copy_r_n[0] = 0.8660254037844386*Dy_r[3]-0.8660254037844386*Dy_r[2]-0.5*Dy_r[1]+0.5*Dy_r[0];
        Dz_con_copy_r_n[0] = 0.8660254037844386*Dz_r[3]-0.8660254037844386*Dz_r[2]-0.5*Dz_r[1]+0.5*Dz_r[0];
        Bx_con_copy_r_n[0] = 0.8660254037844386*Bx_r[3]-0.8660254037844386*Bx_r[2]-0.5*Bx_r[1]+0.5*Bx_r[0];
        By_con_copy_r_n[0] = 0.8660254037844386*By_r[3]-0.8660254037844386*By_r[2]-0.5*By_r[1]+0.5*By_r[0];
        Bz_con_copy_r_n[0] = 0.8660254037844386*Bz_r[3]-0.8660254037844386*Bz_r[2]-0.5*Bz_r[1]+0.5*Bz_r[0];

        Dx_con_copy_r_n[1] = -(0.8660254037844386*Dx_r[3])-0.8660254037844386*Dx_r[2]+0.5*Dx_r[1]+0.5*Dx_r[0];
        Dy_con_copy_r_n[1] = -(0.8660254037844386*Dy_r[3])-0.8660254037844386*Dy_r[2]+0.5*Dy_r[1]+0.5*Dy_r[0];
        Dz_con_copy_r_n[1] = -(0.8660254037844386*Dz_r[3])-0.8660254037844386*Dz_r[2]+0.5*Dz_r[1]+0.5*Dz_r[0];
        Bx_con_copy_r_n[1] = -(0.8660254037844386*Bx_r[3])-0.8660254037844386*Bx_r[2]+0.5*Bx_r[1]+0.5*Bx_r[0];
        By_con_copy_r_n[1] = -(0.8660254037844386*By_r[3])-0.8660254037844386*By_r[2]+0.5*By_r[1]+0.5*By_r[0];
        Bz_con_copy_r_n[1] = -(0.8660254037844386*Bz_r[3])-0.8660254037844386*Bz_r[2]+0.5*Bz_r[1]+0.5*Bz_r[0];

        const double *JDx_l = &field_con_l[0];
        const double *JDy_l = &field_con_l[4];
        const double *JDz_l = &field_con_l[8];
        const double *JBx_l = &field_con_l[12];
        const double *JBy_l = &field_con_l[16];
        const double *JBz_l = &field_con_l[20];

        JDx_con_l_n[0] = -(0.8660254037844386*JDx_l[3])+0.8660254037844386*JDx_l[2]-0.5*JDx_l[1]+0.5*JDx_l[0];
        JDy_con_l_n[0] = -(0.8660254037844386*JDy_l[3])+0.8660254037844386*JDy_l[2]-0.5*JDy_l[1]+0.5*JDy_l[0];
        JDz_con_l_n[0] = -(0.8660254037844386*JDz_l[3])+0.8660254037844386*JDz_l[2]-0.5*JDz_l[1]+0.5*JDz_l[0];
        JBx_con_l_n[0] = -(0.8660254037844386*JBx_l[3])+0.8660254037844386*JBx_l[2]-0.5*JBx_l[1]+0.5*JBx_l[0];
        JBy_con_l_n[0] = -(0.8660254037844386*JBy_l[3])+0.8660254037844386*JBy_l[2]-0.5*JBy_l[1]+0.5*JBy_l[0];
        JBz_con_l_n[0] = -(0.8660254037844386*JBz_l[3])+0.8660254037844386*JBz_l[2]-0.5*JBz_l[1]+0.5*JBz_l[0];

        JDx_con_l_n[1] = 0.8660254037844386*JDx_l[3]+0.8660254037844386*JDx_l[2]+0.5*JDx_l[1]+0.5*JDx_l[0];
        JDy_con_l_n[1] = 0.8660254037844386*JDy_l[3]+0.8660254037844386*JDy_l[2]+0.5*JDy_l[1]+0.5*JDy_l[0];
        JDz_con_l_n[1] = 0.8660254037844386*JDz_l[3]+0.8660254037844386*JDz_l[2]+0.5*JDz_l[1]+0.5*JDz_l[0];
        JBx_con_l_n[1] = 0.8660254037844386*JBx_l[3]+0.8660254037844386*JBx_l[2]+0.5*JBx_l[1]+0.5*JBx_l[0];
        JBy_con_l_n[1] = 0.8660254037844386*JBy_l[3]+0.8660254037844386*JBy_l[2]+0.5*JBy_l[1]+0.5*JBy_l[0];
        JBz_con_l_n[1] = 0.8660254037844386*JBz_l[3]+0.8660254037844386*JBz_l[2]+0.5*JBz_l[1]+0.5*JBz_l[0];

        const double *JDx_r = &field_con_c[0];
        const double *JDy_r = &field_con_c[4];
        const double *JDz_r = &field_con_c[8];
        const double *JBx_r = &field_con_c[12];
        const double *JBy_r = &field_con_c[16];
        const double *JBz_r = &field_con_c[20];

        JDx_con_r_n[0] = 0.8660254037844386*JDx_r[3]-0.8660254037844386*JDx_r[2]-0.5*JDx_r[1]+0.5*JDx_r[0];
        JDy_con_r_n[0] = 0.8660254037844386*JDy_r[3]-0.8660254037844386*JDy_r[2]-0.5*JDy_r[1]+0.5*JDy_r[0];
        JDz_con_r_n[0] = 0.8660254037844386*JDz_r[3]-0.8660254037844386*JDz_r[2]-0.5*JDz_r[1]+0.5*JDz_r[0];
        JBx_con_r_n[0] = 0.8660254037844386*JBx_r[3]-0.8660254037844386*JBx_r[2]-0.5*JBx_r[1]+0.5*JBx_r[0];
        JBy_con_r_n[0] = 0.8660254037844386*JBy_r[3]-0.8660254037844386*JBy_r[2]-0.5*JBy_r[1]+0.5*JBy_r[0];
        JBz_con_r_n[0] = 0.8660254037844386*JBz_r[3]-0.8660254037844386*JBz_r[2]-0.5*JBz_r[1]+0.5*JBz_r[0];

        JDx_con_r_n[1] = -(0.8660254037844386*JDx_r[3])-0.8660254037844386*JDx_r[2]+0.5*JDx_r[1]+0.5*JDx_r[0];
        JDy_con_r_n[1] = -(0.8660254037844386*JDy_r[3])-0.8660254037844386*JDy_r[2]+0.5*JDy_r[1]+0.5*JDy_r[0];
        JDz_con_r_n[1] = -(0.8660254037844386*JDz_r[3])-0.8660254037844386*JDz_r[2]+0.5*JDz_r[1]+0.5*JDz_r[0];
        JBx_con_r_n[1] = -(0.8660254037844386*JBx_r[3])-0.8660254037844386*JBx_r[2]+0.5*JBx_r[1]+0.5*JBx_r[0];
        JBy_con_r_n[1] = -(0.8660254037844386*JBy_r[3])-0.8660254037844386*JBy_r[2]+0.5*JBy_r[1]+0.5*JBy_r[0];
        JBz_con_r_n[1] = -(0.8660254037844386*JBz_r[3])-0.8660254037844386*JBz_r[2]+0.5*JBz_r[1]+0.5*JBz_r[0];

        // Compute the two different Jumps.
        double a_Jc_dU[12];
        double dU[12];
        double a_dQ[12];
        double Q_mag[12];
        double Ul[12];
        double Ur[12];
        for (int j = 0; j<2; ++j) {
          a_Jc_dU[0 + j] = alpha_quad_y[j] * det_h_d[j] * ( Dx_con_r_n[j] - Dx_con_l_n[j] );
          a_Jc_dU[2 + j] = alpha_quad_y[j] * det_h_d[j] * ( Dy_con_r_n[j] - Dy_con_l_n[j] );
          a_Jc_dU[4 + j] = alpha_quad_y[j] * det_h_d[j] * ( Dz_con_r_n[j] - Dz_con_l_n[j] );
          a_Jc_dU[6 + j] = alpha_quad_y[j] * det_h_d[j] * ( Bx_con_r_n[j] - Bx_con_l_n[j] );
          a_Jc_dU[8 + j] = alpha_quad_y[j] * det_h_d[j] * ( By_con_r_n[j] - By_con_l_n[j] );
          a_Jc_dU[10 + j] = alpha_quad_y[j] * det_h_d[j] * ( Bz_con_r_n[j] - Bz_con_l_n[j] );
          dU[0 + j] = ( Dx_con_r_n[j] - Dx_con_l_n[j] );
          dU[2 + j] = ( Dy_con_r_n[j] - Dy_con_l_n[j] );
          dU[4 + j] = ( Dz_con_r_n[j] - Dz_con_l_n[j] );
          dU[6 + j] = ( Bx_con_r_n[j] - Bx_con_l_n[j] );
          dU[8 + j] = ( By_con_r_n[j] - By_con_l_n[j] );
          dU[10 + j] = ( Bz_con_r_n[j] - Bz_con_l_n[j] );
          a_dQ[0 + j] = alpha_quad_y[j] * ( JDx_con_r_n[j] - JDx_con_l_n[j] );
          a_dQ[2 + j] = alpha_quad_y[j] * ( JDy_con_r_n[j] - JDy_con_l_n[j] );
          a_dQ[4 + j] = alpha_quad_y[j] * ( JDz_con_r_n[j] - JDz_con_l_n[j] );
          a_dQ[6 + j] = alpha_quad_y[j] * ( JBx_con_r_n[j] - JBx_con_l_n[j] );
          a_dQ[8 + j] = alpha_quad_y[j] * ( JBy_con_r_n[j] - JBy_con_l_n[j] );
          a_dQ[10 + j] = alpha_quad_y[j] * ( JBz_con_r_n[j] - JBz_con_l_n[j] );
          Q_mag[0 + j] = ( JDx_con_r_n[j] );
          Q_mag[2 + j] = ( JDy_con_r_n[j] );
          Q_mag[4 + j] = ( JDz_con_r_n[j] );
          Q_mag[6 + j] = ( JBx_con_r_n[j] );
          Q_mag[8 + j] = ( JBy_con_r_n[j] );
          Q_mag[10 + j] = ( JBz_con_r_n[j] );
          Ur[0 + j] = ( Dx_con_r_n[j] );
          Ur[2 + j] = ( Dy_con_r_n[j] );
          Ur[4 + j] = ( Dz_con_r_n[j] );
          Ur[6 + j] = ( Bx_con_r_n[j] );
          Ur[8 + j] = ( By_con_r_n[j] );
          Ur[10 + j] = ( Bz_con_r_n[j] );
          Ul[0 + j] = ( Dx_con_l_n[j] );
          Ul[2 + j] = ( Dy_con_l_n[j] );
          Ul[4 + j] = ( Dz_con_l_n[j] );
          Ul[6 + j] = ( Bx_con_l_n[j] );
          Ul[8 + j] = ( By_con_l_n[j] );
          Ul[10 + j] = ( Bz_con_l_n[j] );
        }

        // Specifically verify the theta-pole face states are continuous
        // after applying the same regularity condition as the flux kernel.
        if (theta_pole) {
          for (int j = 0; j<12; ++j) {
            //printf("(theta-pole): dU[%d] %1.16e, Ur: %1.16e, Ul: %1.16e\n", j, dU[j], Ur[j], Ul[j]);
            TEST_CHECK(gkyl_compare_double(dU[j], 0.0, 1e-12));
          }
        }

        // Check the flux-Jump conditions (It is critical this passes!)
        for (int j = 0; j<12; ++j) {
          //printf("(Index: %d) A_plus_dQ_y: %1.16e, A_minus_dQ_y: %1.16e, flux_l_y: %1.16e, flux_r_y: %1.16e\n",j,A_plus_dQ_y[j], A_minus_dQ_y[j], flux_l_y[j], flux_r_y[j]);
        }

        for (int j = 0; j<12; ++j) {
          //printf("(Index: %d) A_plus_dQ_y + A_minus_dQ_y: %1.16e, flux_r_y - flux_l_y: %1.16e\n",j,A_plus_dQ_y[j] + A_minus_dQ_y[j], flux_r_y[j] - flux_l_y[j]);

          // (6/22/2026 - GJ - Temporarily disabling the flux-Jump check until the Roe solver is updated for new Eigensystem including the scalar field corrections.)
          if (theta_pole == 0) {
            //TEST_CHECK(gkyl_compare_double(A_plus_dQ_y[j] + A_minus_dQ_y[j], flux_r_y[j] - flux_l_y[j], 1e-12));
          }
        }

        // print the comparison in the Jumps
        for (int j = 0; j<12; ++j) {
          //printf("(Q_ref[%d] = %1.16e) a_Jc_dU[%d]: %1.16e, a_dQ[%d]: %1.16e\n",j,Q_mag[j],j,a_Jc_dU[j], j, a_dQ[j]);
        }
      }

      // Dir = 0 comparison, for Jumps in quantities at the radial boundary
      if (dir == 0) { 
        // Compute and check |alpha^i| J_c \Delta U ~ |alpha_i| \Delta Q
        double Dx_con_l_n[2] = {0.0};
        double Dy_con_l_n[2] = {0.0};
        double Dz_con_l_n[2] = {0.0};
        double Bx_con_l_n[2] = {0.0};
        double By_con_l_n[2] = {0.0};
        double Bz_con_l_n[2] = {0.0};

        double Dx_con_r_n[2] = {0.0};
        double Dy_con_r_n[2] = {0.0};
        double Dz_con_r_n[2] = {0.0};
        double Bx_con_r_n[2] = {0.0};
        double By_con_r_n[2] = {0.0};
        double Bz_con_r_n[2] = {0.0};

        double JDx_con_l_n[2] = {0.0};
        double JDy_con_l_n[2] = {0.0};
        double JDz_con_l_n[2] = {0.0};
        double JBx_con_l_n[2] = {0.0};
        double JBy_con_l_n[2] = {0.0};
        double JBz_con_l_n[2] = {0.0};

        double JDx_con_r_n[2] = {0.0};
        double JDy_con_r_n[2] = {0.0};
        double JDz_con_r_n[2] = {0.0};
        double JBx_con_r_n[2] = {0.0};
        double JBy_con_r_n[2] = {0.0};
        double JBz_con_r_n[2] = {0.0};

        double Dx_con_copy_l_n[2] = {0.0};
        double Dy_con_copy_l_n[2] = {0.0};
        double Dz_con_copy_l_n[2] = {0.0};
        double Bx_con_copy_l_n[2] = {0.0};
        double By_con_copy_l_n[2] = {0.0};
        double Bz_con_copy_l_n[2] = {0.0};

        double Dx_con_copy_r_n[2] = {0.0};
        double Dy_con_copy_r_n[2] = {0.0};
        double Dz_con_copy_r_n[2] = {0.0};
        double Bx_con_copy_r_n[2] = {0.0};
        double By_con_copy_r_n[2] = {0.0};
        double Bz_con_copy_r_n[2] = {0.0};
        
        const double *Dx_l = &field_no_J_con_l[0]; 
        const double *Dy_l = &field_no_J_con_l[4]; 
        const double *Dz_l = &field_no_J_con_l[8]; 
        const double *Bx_l = &field_no_J_con_l[12]; 
        const double *By_l = &field_no_J_con_l[16]; 
        const double *Bz_l = &field_no_J_con_l[20]; 
        
        Dx_con_l_n[0] = -(0.8660254037844386*Dx_l[3])-0.5*Dx_l[2]+0.8660254037844386*Dx_l[1]+0.5*Dx_l[0];
        Dy_con_l_n[0] = -(0.8660254037844386*Dy_l[3])-0.5*Dy_l[2]+0.8660254037844386*Dy_l[1]+0.5*Dy_l[0];
        Dz_con_l_n[0] = -(0.8660254037844386*Dz_l[3])-0.5*Dz_l[2]+0.8660254037844386*Dz_l[1]+0.5*Dz_l[0];
        Bx_con_l_n[0] = -(0.8660254037844386*Bx_l[3])-0.5*Bx_l[2]+0.8660254037844386*Bx_l[1]+0.5*Bx_l[0];
        By_con_l_n[0] = -(0.8660254037844386*By_l[3])-0.5*By_l[2]+0.8660254037844386*By_l[1]+0.5*By_l[0];
        Bz_con_l_n[0] = -(0.8660254037844386*Bz_l[3])-0.5*Bz_l[2]+0.8660254037844386*Bz_l[1]+0.5*Bz_l[0];
        
        Dx_con_l_n[1] = 0.8660254037844386*Dx_l[3]+0.5*Dx_l[2]+0.8660254037844386*Dx_l[1]+0.5*Dx_l[0];
        Dy_con_l_n[1] = 0.8660254037844386*Dy_l[3]+0.5*Dy_l[2]+0.8660254037844386*Dy_l[1]+0.5*Dy_l[0];
        Dz_con_l_n[1] = 0.8660254037844386*Dz_l[3]+0.5*Dz_l[2]+0.8660254037844386*Dz_l[1]+0.5*Dz_l[0];
        Bx_con_l_n[1] = 0.8660254037844386*Bx_l[3]+0.5*Bx_l[2]+0.8660254037844386*Bx_l[1]+0.5*Bx_l[0];
        By_con_l_n[1] = 0.8660254037844386*By_l[3]+0.5*By_l[2]+0.8660254037844386*By_l[1]+0.5*By_l[0];
        Bz_con_l_n[1] = 0.8660254037844386*Bz_l[3]+0.5*Bz_l[2]+0.8660254037844386*Bz_l[1]+0.5*Bz_l[0];

        const double *Dx_r = &field_no_J_con_c[0]; 
        const double *Dy_r = &field_no_J_con_c[4]; 
        const double *Dz_r = &field_no_J_con_c[8]; 
        const double *Bx_r = &field_no_J_con_c[12]; 
        const double *By_r = &field_no_J_con_c[16]; 
        const double *Bz_r = &field_no_J_con_c[20]; 
        
        Dx_con_r_n[0] = 0.8660254037844386*Dx_r[3]-0.5*Dx_r[2]-0.8660254037844386*Dx_r[1]+0.5*Dx_r[0];
        Dy_con_r_n[0] = 0.8660254037844386*Dy_r[3]-0.5*Dy_r[2]-0.8660254037844386*Dy_r[1]+0.5*Dy_r[0];
        Dz_con_r_n[0] = 0.8660254037844386*Dz_r[3]-0.5*Dz_r[2]-0.8660254037844386*Dz_r[1]+0.5*Dz_r[0];
        Bx_con_r_n[0] = 0.8660254037844386*Bx_r[3]-0.5*Bx_r[2]-0.8660254037844386*Bx_r[1]+0.5*Bx_r[0];
        By_con_r_n[0] = 0.8660254037844386*By_r[3]-0.5*By_r[2]-0.8660254037844386*By_r[1]+0.5*By_r[0];
        Bz_con_r_n[0] = 0.8660254037844386*Bz_r[3]-0.5*Bz_r[2]-0.8660254037844386*Bz_r[1]+0.5*Bz_r[0];
        
        Dx_con_r_n[1] = -(0.8660254037844386*Dx_r[3])+0.5*Dx_r[2]-0.8660254037844386*Dx_r[1]+0.5*Dx_r[0];
        Dy_con_r_n[1] = -(0.8660254037844386*Dy_r[3])+0.5*Dy_r[2]-0.8660254037844386*Dy_r[1]+0.5*Dy_r[0];
        Dz_con_r_n[1] = -(0.8660254037844386*Dz_r[3])+0.5*Dz_r[2]-0.8660254037844386*Dz_r[1]+0.5*Dz_r[0];
        Bx_con_r_n[1] = -(0.8660254037844386*Bx_r[3])+0.5*Bx_r[2]-0.8660254037844386*Bx_r[1]+0.5*Bx_r[0];
        By_con_r_n[1] = -(0.8660254037844386*By_r[3])+0.5*By_r[2]-0.8660254037844386*By_r[1]+0.5*By_r[0];
        Bz_con_r_n[1] = -(0.8660254037844386*Bz_r[3])+0.5*Bz_r[2]-0.8660254037844386*Bz_r[1]+0.5*Bz_r[0];

        // Acqiure new pointer locations
        Dx_l = &field_no_J_con_copy_l[0]; 
        Dy_l = &field_no_J_con_copy_l[4]; 
        Dz_l = &field_no_J_con_copy_l[8]; 
        Bx_l = &field_no_J_con_copy_l[12]; 
        By_l = &field_no_J_con_copy_l[16]; 
        Bz_l = &field_no_J_con_copy_l[20]; 
        
        Dx_con_copy_l_n[0] = -(0.8660254037844386*Dx_l[3])-0.5*Dx_l[2]+0.8660254037844386*Dx_l[1]+0.5*Dx_l[0];
        Dy_con_copy_l_n[0] = -(0.8660254037844386*Dy_l[3])-0.5*Dy_l[2]+0.8660254037844386*Dy_l[1]+0.5*Dy_l[0];
        Dz_con_copy_l_n[0] = -(0.8660254037844386*Dz_l[3])-0.5*Dz_l[2]+0.8660254037844386*Dz_l[1]+0.5*Dz_l[0];
        Bx_con_copy_l_n[0] = -(0.8660254037844386*Bx_l[3])-0.5*Bx_l[2]+0.8660254037844386*Bx_l[1]+0.5*Bx_l[0];
        By_con_copy_l_n[0] = -(0.8660254037844386*By_l[3])-0.5*By_l[2]+0.8660254037844386*By_l[1]+0.5*By_l[0];
        Bz_con_copy_l_n[0] = -(0.8660254037844386*Bz_l[3])-0.5*Bz_l[2]+0.8660254037844386*Bz_l[1]+0.5*Bz_l[0];
        
        Dx_con_copy_l_n[1] = 0.8660254037844386*Dx_l[3]+0.5*Dx_l[2]+0.8660254037844386*Dx_l[1]+0.5*Dx_l[0];
        Dy_con_copy_l_n[1] = 0.8660254037844386*Dy_l[3]+0.5*Dy_l[2]+0.8660254037844386*Dy_l[1]+0.5*Dy_l[0];
        Dz_con_copy_l_n[1] = 0.8660254037844386*Dz_l[3]+0.5*Dz_l[2]+0.8660254037844386*Dz_l[1]+0.5*Dz_l[0];
        Bx_con_copy_l_n[1] = 0.8660254037844386*Bx_l[3]+0.5*Bx_l[2]+0.8660254037844386*Bx_l[1]+0.5*Bx_l[0];
        By_con_copy_l_n[1] = 0.8660254037844386*By_l[3]+0.5*By_l[2]+0.8660254037844386*By_l[1]+0.5*By_l[0];
        Bz_con_copy_l_n[1] = 0.8660254037844386*Bz_l[3]+0.5*Bz_l[2]+0.8660254037844386*Bz_l[1]+0.5*Bz_l[0];

        // Acqiure new pointer locations
        Dx_r = &field_no_J_con_copy_c[0]; 
        Dy_r = &field_no_J_con_copy_c[4]; 
        Dz_r = &field_no_J_con_copy_c[8]; 
        Bx_r = &field_no_J_con_copy_c[12]; 
        By_r = &field_no_J_con_copy_c[16]; 
        Bz_r = &field_no_J_con_copy_c[20]; 
        
        Dx_con_copy_r_n[0] = 0.8660254037844386*Dx_r[3]-0.5*Dx_r[2]-0.8660254037844386*Dx_r[1]+0.5*Dx_r[0];
        Dy_con_copy_r_n[0] = 0.8660254037844386*Dy_r[3]-0.5*Dy_r[2]-0.8660254037844386*Dy_r[1]+0.5*Dy_r[0];
        Dz_con_copy_r_n[0] = 0.8660254037844386*Dz_r[3]-0.5*Dz_r[2]-0.8660254037844386*Dz_r[1]+0.5*Dz_r[0];
        Bx_con_copy_r_n[0] = 0.8660254037844386*Bx_r[3]-0.5*Bx_r[2]-0.8660254037844386*Bx_r[1]+0.5*Bx_r[0];
        By_con_copy_r_n[0] = 0.8660254037844386*By_r[3]-0.5*By_r[2]-0.8660254037844386*By_r[1]+0.5*By_r[0];
        Bz_con_copy_r_n[0] = 0.8660254037844386*Bz_r[3]-0.5*Bz_r[2]-0.8660254037844386*Bz_r[1]+0.5*Bz_r[0];
        
        Dx_con_copy_r_n[1] = -(0.8660254037844386*Dx_r[3])+0.5*Dx_r[2]-0.8660254037844386*Dx_r[1]+0.5*Dx_r[0];
        Dy_con_copy_r_n[1] = -(0.8660254037844386*Dy_r[3])+0.5*Dy_r[2]-0.8660254037844386*Dy_r[1]+0.5*Dy_r[0];
        Dz_con_copy_r_n[1] = -(0.8660254037844386*Dz_r[3])+0.5*Dz_r[2]-0.8660254037844386*Dz_r[1]+0.5*Dz_r[0];
        Bx_con_copy_r_n[1] = -(0.8660254037844386*Bx_r[3])+0.5*Bx_r[2]-0.8660254037844386*Bx_r[1]+0.5*Bx_r[0];
        By_con_copy_r_n[1] = -(0.8660254037844386*By_r[3])+0.5*By_r[2]-0.8660254037844386*By_r[1]+0.5*By_r[0];
        Bz_con_copy_r_n[1] = -(0.8660254037844386*Bz_r[3])+0.5*Bz_r[2]-0.8660254037844386*Bz_r[1]+0.5*Bz_r[0];
        
        const double *JDx_l = &field_con_l[0]; 
        const double *JDy_l = &field_con_l[4]; 
        const double *JDz_l = &field_con_l[8]; 
        const double *JBx_l = &field_con_l[12]; 
        const double *JBy_l = &field_con_l[16]; 
        const double *JBz_l = &field_con_l[20]; 
        
        JDx_con_l_n[0] = -(0.8660254037844386*JDx_l[3])-0.5*JDx_l[2]+0.8660254037844386*JDx_l[1]+0.5*JDx_l[0];
        JDy_con_l_n[0] = -(0.8660254037844386*JDy_l[3])-0.5*JDy_l[2]+0.8660254037844386*JDy_l[1]+0.5*JDy_l[0];
        JDz_con_l_n[0] = -(0.8660254037844386*JDz_l[3])-0.5*JDz_l[2]+0.8660254037844386*JDz_l[1]+0.5*JDz_l[0];
        JBx_con_l_n[0] = -(0.8660254037844386*JBx_l[3])-0.5*JBx_l[2]+0.8660254037844386*JBx_l[1]+0.5*JBx_l[0];
        JBy_con_l_n[0] = -(0.8660254037844386*JBy_l[3])-0.5*JBy_l[2]+0.8660254037844386*JBy_l[1]+0.5*JBy_l[0];
        JBz_con_l_n[0] = -(0.8660254037844386*JBz_l[3])-0.5*JBz_l[2]+0.8660254037844386*JBz_l[1]+0.5*JBz_l[0];
        
        JDx_con_l_n[1] = 0.8660254037844386*JDx_l[3]+0.5*JDx_l[2]+0.8660254037844386*JDx_l[1]+0.5*JDx_l[0];
        JDy_con_l_n[1] = 0.8660254037844386*JDy_l[3]+0.5*JDy_l[2]+0.8660254037844386*JDy_l[1]+0.5*JDy_l[0];
        JDz_con_l_n[1] = 0.8660254037844386*JDz_l[3]+0.5*JDz_l[2]+0.8660254037844386*JDz_l[1]+0.5*JDz_l[0];
        JBx_con_l_n[1] = 0.8660254037844386*JBx_l[3]+0.5*JBx_l[2]+0.8660254037844386*JBx_l[1]+0.5*JBx_l[0];
        JBy_con_l_n[1] = 0.8660254037844386*JBy_l[3]+0.5*JBy_l[2]+0.8660254037844386*JBy_l[1]+0.5*JBy_l[0];
        JBz_con_l_n[1] = 0.8660254037844386*JBz_l[3]+0.5*JBz_l[2]+0.8660254037844386*JBz_l[1]+0.5*JBz_l[0];

        const double *JDx_r = &field_con_c[0]; 
        const double *JDy_r = &field_con_c[4]; 
        const double *JDz_r = &field_con_c[8]; 
        const double *JBx_r = &field_con_c[12]; 
        const double *JBy_r = &field_con_c[16]; 
        const double *JBz_r = &field_con_c[20]; 
        
        JDx_con_r_n[0] = 0.8660254037844386*JDx_r[3]-0.5*JDx_r[2]-0.8660254037844386*JDx_r[1]+0.5*JDx_r[0];
        JDy_con_r_n[0] = 0.8660254037844386*JDy_r[3]-0.5*JDy_r[2]-0.8660254037844386*JDy_r[1]+0.5*JDy_r[0];
        JDz_con_r_n[0] = 0.8660254037844386*JDz_r[3]-0.5*JDz_r[2]-0.8660254037844386*JDz_r[1]+0.5*JDz_r[0];
        JBx_con_r_n[0] = 0.8660254037844386*JBx_r[3]-0.5*JBx_r[2]-0.8660254037844386*JBx_r[1]+0.5*JBx_r[0];
        JBy_con_r_n[0] = 0.8660254037844386*JBy_r[3]-0.5*JBy_r[2]-0.8660254037844386*JBy_r[1]+0.5*JBy_r[0];
        JBz_con_r_n[0] = 0.8660254037844386*JBz_r[3]-0.5*JBz_r[2]-0.8660254037844386*JBz_r[1]+0.5*JBz_r[0];
        
        JDx_con_r_n[1] = -(0.8660254037844386*JDx_r[3])+0.5*JDx_r[2]-0.8660254037844386*JDx_r[1]+0.5*JDx_r[0];
        JDy_con_r_n[1] = -(0.8660254037844386*JDy_r[3])+0.5*JDy_r[2]-0.8660254037844386*JDy_r[1]+0.5*JDy_r[0];
        JDz_con_r_n[1] = -(0.8660254037844386*JDz_r[3])+0.5*JDz_r[2]-0.8660254037844386*JDz_r[1]+0.5*JDz_r[0];
        JBx_con_r_n[1] = -(0.8660254037844386*JBx_r[3])+0.5*JBx_r[2]-0.8660254037844386*JBx_r[1]+0.5*JBx_r[0];
        JBy_con_r_n[1] = -(0.8660254037844386*JBy_r[3])+0.5*JBy_r[2]-0.8660254037844386*JBy_r[1]+0.5*JBy_r[0];
        JBz_con_r_n[1] = -(0.8660254037844386*JBz_r[3])+0.5*JBz_r[2]-0.8660254037844386*JBz_r[1]+0.5*JBz_r[0];

        // Compute the two different Jumps
        // Iterate over nodes
        double a_Jc_dU[12]; 
        double dU_copy[12]; 
        double a_dQ[12]; 
        double Q_mag[12]; 
        double U_mag[12]; 
        for (int j = 0; j<2; ++j) {
          a_Jc_dU[0 + j] = alpha_quad_x[j] * det_h_d[j] * ( Dx_con_r_n[j] - Dx_con_l_n[j] ); 
          a_Jc_dU[2 + j] = alpha_quad_x[j] * det_h_d[j] * ( Dy_con_r_n[j] - Dy_con_l_n[j] ); 
          a_Jc_dU[4 + j] = alpha_quad_x[j] * det_h_d[j] * ( Dz_con_r_n[j] - Dz_con_l_n[j] ); 
          a_Jc_dU[6 + j] = alpha_quad_x[j] * det_h_d[j] * ( Bx_con_r_n[j] - Bx_con_l_n[j] ); 
          a_Jc_dU[8 + j] = alpha_quad_x[j] * det_h_d[j] * ( By_con_r_n[j] - By_con_l_n[j] ); 
          a_Jc_dU[10 + j] = alpha_quad_x[j] * det_h_d[j] * ( Bz_con_r_n[j] - Bz_con_l_n[j] ); 
          dU_copy[0 + j] = ( Dx_con_copy_r_n[j] - Dx_con_copy_l_n[j] ); 
          dU_copy[2 + j] = ( Dy_con_copy_r_n[j] - Dy_con_copy_l_n[j] ); 
          dU_copy[4 + j] = ( Dz_con_copy_r_n[j] - Dz_con_copy_l_n[j] ); 
          dU_copy[6 + j] = ( Bx_con_copy_r_n[j] - Bx_con_copy_l_n[j] ); 
          dU_copy[8 + j] = ( By_con_copy_r_n[j] - By_con_copy_l_n[j] ); 
          dU_copy[10 + j] = ( Bz_con_copy_r_n[j] - Bz_con_copy_l_n[j] ); 
          a_dQ[0 + j] = alpha_quad_x[j] * ( JDx_con_r_n[j] - JDx_con_l_n[j] ); 
          a_dQ[2 + j] = alpha_quad_x[j] * ( JDy_con_r_n[j] - JDy_con_l_n[j] ); 
          a_dQ[4 + j] = alpha_quad_x[j] * ( JDz_con_r_n[j] - JDz_con_l_n[j] ); 
          a_dQ[6 + j] = alpha_quad_x[j] * ( JBx_con_r_n[j] - JBx_con_l_n[j] ); 
          a_dQ[8 + j] = alpha_quad_x[j] * ( JBy_con_r_n[j] - JBy_con_l_n[j] ); 
          a_dQ[10 + j] = alpha_quad_x[j] * ( JBz_con_r_n[j] - JBz_con_l_n[j] ); 
          Q_mag[0 + j] = ( JDx_con_r_n[j] ); 
          Q_mag[2 + j] = ( JDy_con_r_n[j] ); 
          Q_mag[4 + j] = ( JDz_con_r_n[j] ); 
          Q_mag[6 + j] = ( JBx_con_r_n[j] ); 
          Q_mag[8 + j] = ( JBy_con_r_n[j] ); 
          Q_mag[10 + j] = ( JBz_con_r_n[j] ); 
          U_mag[0 + j] = ( Dx_con_r_n[j] ); 
          U_mag[2 + j] = ( Dy_con_r_n[j] ); 
          U_mag[4 + j] = ( Dz_con_r_n[j] ); 
          U_mag[6 + j] = ( Bx_con_r_n[j] ); 
          U_mag[8 + j] = ( By_con_r_n[j] ); 
          U_mag[10 + j] = ( Bz_con_r_n[j] );
        }

        // Check the flux-Jump conditions (It is critical this passes!)
        for (int j = 0; j<12; ++j) {
          // printf("(Index: %d) A_plus_dQ_x: %1.16e, A_minus_dQ_x: %1.16e, flux_l_x: %1.16e, flux_r_x: %1.16e\n",j,A_plus_dQ_x[j], A_minus_dQ_x[j], flux_l_x[j], flux_r_x[j]);
          //printf("(Index: %d) A_plus_dQ_x + A_minus_dQ_x: %1.16e, flux_r_x - flux_l_x: %1.16e\n",j,A_plus_dQ_x[j] + A_minus_dQ_x[j], flux_r_x[j] - flux_l_x[j]);
          //TEST_CHECK(gkyl_compare_double(A_plus_dQ_x[j] + A_minus_dQ_x[j], flux_r_x[j] - flux_l_x[j], 1e-12));
        }

        // print the comparison in the Jumps
        for (int j = 0; j<12; ++j) {
          //printf("(Q_ref[%d] = %1.16e) a_Jc_dU[%d]: %1.16e, a_dQ[%d]: %1.16e\n",j,Q_mag[j],j,a_Jc_dU[j], j, a_dQ[j]);
        }
        // print the comparison in the Jumps
        for (int j = 0; j<12; ++j) {
          //printf("(U_ref[%d] = %1.16e) dU[%d]: %1.16e\n",j,U_mag[j],j,dU_copy[j]);
        }
      }

      // Check the charateristic speeds at quadratures
      if (dir == 0) {
        for (int j = 0; j<2; ++j) {
          
          // Eigenvalues at the inner radial boundary are [0.0, -1.0, and -0.33], so |alpha| is 1:
          //printf("[X] alpha(Node: %d): alpha: %1.16e \n",j,alpha_quad_x[j]);
          TEST_CHECK(gkyl_compare_double(alpha_quad_x[j], 1.0, 1e-12));

          // Check the field components are properly computed.
          for (int i = 0; i<6; ++i) { 
            //printf("[X] Field(Comp: %d, Node: %d): Field_l: %1.16e, Field_r: %1.16e \n",i,j, field_no_J_con_l[j + i*2], field_no_J_con_c[j + i*2]);
            //printf("[X] Field(Comp: %d, Node: %d): Field_l: %1.16e, Field_r: %1.16e \n",i,j, field_no_J_con_l[j + i*2], field_no_J_con_c[j + i*2]);
            //printf("[X] Flux(Comp: %d, Node: %d): flux_l: %1.16e, flux_r: %1.16e \n",i,j, flux_l_x[j + i*2], flux_r_x[j + i*2]);
          }
        }
      }
    }
  }

  // Release vmg data
  gkyl_vlasov_position_map_release(pos_map);
  gkyl_surf_and_vol_node_arrays_release(lapse);
  gkyl_surf_and_vol_node_arrays_release(shift);
  gkyl_surf_and_vol_node_arrays_release(h_ij);
  gkyl_surf_and_vol_node_arrays_release(h_ij_inv);
  gkyl_surf_and_vol_node_arrays_release(det_h);
  
  
  // Release field arrays
  gkyl_array_release(field_with_J_con);
  gkyl_array_release(field_no_J_con);
  gkyl_array_release(conf_flux_surf);
  gkyl_array_release(cflrate);

  // Release fields
  gkyl_array_release(field_no_J_fixed);
  gkyl_array_release(field_no_J_fixed_copy);

}

void test_ks_r_theta_2x_p1() { test_ks_r_theta_2x_geom_p1(); }

TEST_LIST = {
  { "test_ks_r_theta_2x_p1", test_ks_r_theta_2x_p1},
  {NULL, NULL}
};
