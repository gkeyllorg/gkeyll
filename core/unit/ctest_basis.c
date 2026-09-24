#include <math.h>

#include <acutest.h>
#include <gkyl_basis.h>
#include <gkyl_util.h>
#include <gkyl_basis_gkhyb_1x2v_p1_surfx1_eval_quad.h>
#include <gkyl_basis_gkhyb_1x2v_p1_surfx2_eval_quad.h>
#include <gkyl_basis_gkhyb_1x2v_p1_surfx3_eval_quad.h>
#include <gkyl_basis_gkhyb_1x2v_p1_upwind_quad_to_modal.h>

void
test_ser_1d_members(struct gkyl_basis basis1)
{
  TEST_CHECK( basis1.ndim == 1 );
  TEST_CHECK( basis1.poly_order == 1 );
  TEST_CHECK( basis1.num_basis == 2 );
  TEST_CHECK( strcmp(basis1.id, "serendipity") == 0 );
  TEST_CHECK( basis1.b_type == GKYL_BASIS_MODAL_SERENDIPITY );

  double z[basis1.num_basis], b[basis1.num_basis];

  z[0] = 0.0;
  basis1.eval(z, b);

  TEST_CHECK( gkyl_compare(b[0], 1/sqrt(2.0), 1e-15) );
  TEST_CHECK( b[1] == 0.0 );

  z[0] = 0.5; basis1.eval(z, b);
  
  TEST_CHECK( gkyl_compare(b[0], 1/sqrt(2.0), 1e-15) );
  TEST_CHECK( b[1] == sqrt(3.0/2.0)*0.5 );

  double nodes[basis1.ndim*basis1.num_basis];
  basis1.node_list(nodes);

  TEST_CHECK( nodes[0] == -1 );
  TEST_CHECK( nodes[1] == 1 );

  double f[basis1.num_basis];
  for (int i=0; i<basis1.num_basis; ++i) f[i] = 1.0;

  z[0] = 0.5;
  TEST_CHECK ( gkyl_compare(1.319479216882342, basis1.eval_expand(z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(1.224744871391589, basis1.eval_grad_expand(0, z, f), 1e-15) );  

  z[0] = 0.25;
  TEST_CHECK ( gkyl_compare(1.013292999034445, basis1.eval_expand(z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(1.224744871391589, basis1.eval_grad_expand(0, z, f), 1e-15) );
}

void
test_ser_1d()
{
  struct gkyl_basis basis1;
  gkyl_cart_modal_serendip(&basis1, 1, 1);
  test_ser_1d_members(basis1);

  struct gkyl_basis *basis2 = gkyl_cart_modal_serendip_new(1, 1);
  TEST_CHECK(basis1.num_basis == gkyl_cart_modal_basis_get_num_basis(basis2));
  test_ser_1d_members(*basis2);
  gkyl_cart_modal_basis_release(basis2);
}

void
test_ser_2d_members(struct gkyl_basis basis)
{
  TEST_CHECK( basis.ndim == 2 );
  TEST_CHECK( basis.poly_order == 2 );
  TEST_CHECK( basis.num_basis == 8 );
  TEST_CHECK( strcmp(basis.id, "serendipity") == 0 );
  TEST_CHECK( basis.b_type == GKYL_BASIS_MODAL_SERENDIPITY );

  double z[basis.ndim], b[basis.num_basis];

  z[0] = 0.0; z[1] = 0.0;
  basis.eval(z, b);

  TEST_CHECK( gkyl_compare(0.5, b[0], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[1], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[2], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[3], 1e-15) );
  TEST_CHECK( gkyl_compare(-sqrt(5.0)/4, b[4], 1e-15) );
  TEST_CHECK( gkyl_compare(-sqrt(5.0)/4, b[5], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[6], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[7], 1e-15) );

  double fin[basis.num_basis], fout[basis.num_basis];
  for (int i=0; i<basis.num_basis; ++i) {
    fin[i] = 1.0;
    fout[i] = 0.0;
  }

  basis.flip_odd_sign(0, fin, fout);
  TEST_CHECK( fin[0] == fout[0] );
  TEST_CHECK( -fin[1] == fout[1] );
  TEST_CHECK( fin[2] == fout[2] );  
  TEST_CHECK( -fin[3] == fout[3] );
  TEST_CHECK( fin[4] == fout[4] );
  TEST_CHECK( fin[5] == fout[5] );
  TEST_CHECK( fin[6] == fout[6] );
  TEST_CHECK( -fin[7] == fout[7] );

  basis.flip_odd_sign(1, fin, fout);
  TEST_CHECK( fin[0] == fout[0] );
  TEST_CHECK( fin[1] == fout[1] );
  TEST_CHECK( -fin[2] == fout[2] );  
  TEST_CHECK( -fin[3] == fout[3] );
  TEST_CHECK( fin[4] == fout[4] );
  TEST_CHECK( fin[5] == fout[5] );
  TEST_CHECK( -fin[6] == fout[6] );
  TEST_CHECK( fin[7] == fout[7] );

  basis.flip_even_sign(0, fin, fout);
  TEST_CHECK( -fin[0] == fout[0] );
  TEST_CHECK( fin[1] == fout[1] );
  TEST_CHECK( -fin[2] == fout[2] );  
  TEST_CHECK( fin[3] == fout[3] );
  TEST_CHECK( -fin[4] == fout[4] );
  TEST_CHECK( -fin[5] == fout[5] );
  TEST_CHECK( -fin[6] == fout[6] );
  TEST_CHECK( fin[7] == fout[7] );

  basis.flip_even_sign(1, fin, fout);
  TEST_CHECK( -fin[0] == fout[0] );
  TEST_CHECK( -fin[1] == fout[1] );
  TEST_CHECK( fin[2] == fout[2] );  
  TEST_CHECK( fin[3] == fout[3] );
  TEST_CHECK( -fin[4] == fout[4] );
  TEST_CHECK( -fin[5] == fout[5] );
  TEST_CHECK( fin[6] == fout[6] );
  TEST_CHECK( -fin[7] == fout[7] );  

  double nodes[basis.ndim*basis.num_basis];
  basis.node_list(nodes);

  TEST_CHECK( nodes[0] == -1 );
  TEST_CHECK( nodes[1] == -1 );

  TEST_CHECK( nodes[2] == 0 );
  TEST_CHECK( nodes[3] == -1 );

  TEST_CHECK( nodes[4] == 1 );
  TEST_CHECK( nodes[5] == -1 );

  double f[basis.num_basis];
  for (int i=0; i<basis.num_basis; ++i) f[i] = 1.0;

  z[0] = 0.5; z[1] = 0.5;
  TEST_CHECK ( gkyl_compare(1.219455447459001, basis.eval_expand(z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(4.503383682599098, basis.eval_grad_expand(0, z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(4.503383682599098, basis.eval_grad_expand(1, z, f), 1e-15) );

  z[0] = 0.25; z[1] = -0.75;
  TEST_CHECK ( gkyl_compare(0.4723022336170485, basis.eval_expand(z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(0.1559433418554234, basis.eval_grad_expand(0, z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(-3.150527379222043, basis.eval_grad_expand(1, z, f), 1e-15) );
}

void
test_ser_2d()
{
  struct gkyl_basis basis1;
  gkyl_cart_modal_serendip(&basis1, 2, 2);
  test_ser_2d_members(basis1);

  struct gkyl_basis *basis2 = gkyl_cart_modal_serendip_new(2, 2);
  test_ser_2d_members(*basis2);
  gkyl_cart_modal_basis_release(basis2);
}

void
test_ten_2d_members(struct gkyl_basis basis)
{
  TEST_CHECK( basis.ndim == 2 );
  TEST_CHECK( basis.poly_order == 2 );
  TEST_CHECK( basis.num_basis == 9 );
  TEST_CHECK( strcmp(basis.id, "tensor") == 0 );
  TEST_CHECK( basis.b_type == GKYL_BASIS_MODAL_TENSOR );

  double z[basis.ndim], b[basis.num_basis];

  z[0] = 0.0; z[1] = 0.0;
  basis.eval(z, b);

  TEST_CHECK( gkyl_compare(0.5, b[0], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[1], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[2], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[3], 1e-15) );
  TEST_CHECK( gkyl_compare(-sqrt(5.0)/4, b[4], 1e-15) );
  TEST_CHECK( gkyl_compare(-sqrt(5.0)/4, b[5], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[6], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[7], 1e-15) );
  TEST_CHECK( gkyl_compare(5.0/8.0, b[8], 1e-15) );

  double f[basis.num_basis];
  for (int i=0; i<basis.num_basis; ++i) f[i] = 1.0;

  z[0] = 0.5; z[1] = 0.5;
  TEST_CHECK ( gkyl_compare(1.258517947459001, basis.eval_expand(z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(4.034633682599098, basis.eval_grad_expand(0, z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(4.034633682599098, basis.eval_grad_expand(1, z, f), 1e-15) );  

  z[0] = 0.25; z[1] = -0.75;
  TEST_CHECK ( gkyl_compare(0.1231811398670484, basis.eval_expand(z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(0.8004745918554234, basis.eval_grad_expand(0, z, f), 1e-15) );
  TEST_CHECK ( gkyl_compare(-0.8653711292220426, basis.eval_grad_expand(1, z, f), 1e-15) );
}

void
test_ten_2d_members_p3(struct gkyl_basis basis)
{
  TEST_CHECK( basis.ndim == 2 );
  TEST_CHECK( basis.poly_order == 3 );
  TEST_CHECK( basis.num_basis == 16 );
  TEST_CHECK( strcmp(basis.id, "tensor") == 0 );
  TEST_CHECK( basis.b_type == GKYL_BASIS_MODAL_TENSOR );

  double z[basis.ndim], b[basis.num_basis];

  z[0] = 0.0; z[1] = 0.0;
  basis.eval(z, b);

  TEST_CHECK( gkyl_compare(0.5, b[0], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[1], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[2], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[3], 1e-15) );
  TEST_CHECK( gkyl_compare(-sqrt(5.0)/4, b[4], 1e-15) );
  TEST_CHECK( gkyl_compare(-sqrt(5.0)/4, b[5], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[6], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[7], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[8], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[9], 1e-15) );
  TEST_CHECK( gkyl_compare(5.0/8.0, b[10], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[11], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[12], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[13], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[14], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[15], 1e-15) );

  double f[basis.num_basis];
  for (int i=0; i<basis.num_basis; ++i) f[i] = 1.0;

  z[0] = 0.5; z[1] = 0.5;
  TEST_CHECK ( gkyl_compare(0.09202080373491306, basis.eval_expand(z, f), 1e-14) );
  TEST_CHECK ( gkyl_compare(1.303799542808271, basis.eval_grad_expand(0, z, f), 1e-14) );
  TEST_CHECK ( gkyl_compare(1.303799542808271, basis.eval_grad_expand(1, z, f), 1e-14) );

  z[0] = 0.25; z[1] = -0.75;
  TEST_CHECK ( gkyl_compare(-0.1193909952935715, basis.eval_expand(z, f), 1e-14) );
  TEST_CHECK ( gkyl_compare(0.2231373667479314, basis.eval_grad_expand(0, z, f), 1e-14) );
  TEST_CHECK ( gkyl_compare(-0.7090977834887839, basis.eval_grad_expand(1, z, f), 1e-14) );
}

void
test_ten_2d()
{
  struct gkyl_basis basis1;
  gkyl_cart_modal_tensor(&basis1, 2, 2);
  test_ten_2d_members(basis1);

  struct gkyl_basis *basis2 = gkyl_cart_modal_tensor_new(2, 2);
  test_ten_2d_members(*basis2);
  gkyl_cart_modal_basis_release(basis2);

  struct gkyl_basis basis3;
  gkyl_cart_modal_tensor(&basis3, 2, 3);
  test_ten_2d_members_p3(basis3);
}

void
test_hyb_members(struct gkyl_basis basis)
{
  TEST_CHECK( basis.ndim == 2 );
  TEST_CHECK( basis.poly_order == 1 );
  TEST_CHECK( basis.num_basis == 6 );
  TEST_CHECK( strcmp(basis.id, "hybrid") == 0 );
  TEST_CHECK( basis.b_type == GKYL_BASIS_MODAL_HYBRID );

  double z[basis.ndim], b[basis.num_basis];

  z[0] = 0.0; z[1] = 0.0;
  basis.eval(z, b);

  TEST_CHECK( gkyl_compare(0.5, b[0], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[1], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[2], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[3], 1e-15) );
  TEST_CHECK( gkyl_compare(-sqrt(5.0)/4, b[4], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[6], 1e-15) );

  double fin[basis.num_basis], fout[basis.num_basis];
  for (int i=0; i<basis.num_basis; ++i) {
    fin[i] = 1.0;
    fout[i] = 0.0;
  }

  basis.flip_odd_sign(0, fin, fout);
  TEST_CHECK( fin[0] == fout[0] );
  TEST_CHECK( -fin[1] == fout[1] );
  TEST_CHECK( fin[2] == fout[2] );  
  TEST_CHECK( -fin[3] == fout[3] );
  TEST_CHECK( fin[4] == fout[4] );
  TEST_CHECK( -fin[5] == fout[5] );

}

void
test_hyb()
{
  struct gkyl_basis basis1;
  gkyl_cart_modal_hybrid(&basis1, 1, 1);
  test_hyb_members(basis1);

  struct gkyl_basis *basis2 = gkyl_cart_modal_hybrid_new(1, 1);
  test_hyb_members(*basis2);
  gkyl_cart_modal_basis_release(basis2);
}

void
test_gkhyb_members(struct gkyl_basis basis)
{
  TEST_CHECK( basis.ndim == 3 );
  TEST_CHECK( basis.poly_order == 1 );
  TEST_CHECK( basis.num_basis == 12 );
  TEST_CHECK( strcmp(basis.id, "gkhybrid") == 0 );
  TEST_CHECK( basis.b_type == GKYL_BASIS_MODAL_GKHYBRID );

  double z[basis.ndim], b[basis.num_basis];

  z[0] = 0.0; z[1] = 0.0; z[2] = 0.0;
  basis.eval(z, b);

  TEST_CHECK( gkyl_compare(1./sqrt(pow(2.,3)), b[0], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[1], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[2], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[3], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[4], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[5], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[6], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[7], 1e-15) );
  TEST_CHECK( gkyl_compare(-sqrt(5.0)/sqrt(pow(2.,5)), b[8], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[9], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[10], 1e-15) );
  TEST_CHECK( gkyl_compare(0.0, b[11], 1e-15) );

  double fin[basis.num_basis], fout[basis.num_basis];
  for (int i=0; i<basis.num_basis; ++i) {
    fin[i] = 1.0;
    fout[i] = 0.0;
  }

  basis.flip_odd_sign(0, fin, fout);
  TEST_CHECK( fin[0] == fout[0] );
  TEST_CHECK( -fin[1] == fout[1] );
  TEST_CHECK( fin[2] == fout[2] );  
  TEST_CHECK( fin[3] == fout[3] );  
  TEST_CHECK( -fin[4] == fout[4] );
  TEST_CHECK( -fin[5] == fout[5] );
  TEST_CHECK( fin[6] == fout[6] );
  TEST_CHECK( -fin[7] == fout[7] );
  TEST_CHECK( fin[8] == fout[8] );
  TEST_CHECK( -fin[9] == fout[9] );
  TEST_CHECK( fin[10] == fout[10] );
  TEST_CHECK( -fin[11] == fout[11] );
}

void
test_gkhyb_1x2v_surf_eval_nod(struct gkyl_basis basis)
{
  // Accepted results here were generated with
  // ms-basis_surf_quad_upwind_accepted_results.mac
  double fin[basis.num_basis];
  for (int i=0; i<basis.num_basis; ++i)
    fin[i] = i+1.0;

  typedef double (*doubleFunc_t)(const double* GKYL_RESTRICT f);
  int numnod;

  // Evaluate at left nodes on surf perp to dir 1.
  doubleFunc_t funcs1l[] = {gkhyb_1x2v_p1_surfx1_eval_quad_node_0_l,
                            gkhyb_1x2v_p1_surfx1_eval_quad_node_1_l,
                            gkhyb_1x2v_p1_surfx1_eval_quad_node_2_l,
                            gkhyb_1x2v_p1_surfx1_eval_quad_node_3_l,
                            gkhyb_1x2v_p1_surfx1_eval_quad_node_4_l,
                            gkhyb_1x2v_p1_surfx1_eval_quad_node_5_l};
  numnod = sizeof(funcs1l)/sizeof(doubleFunc_t);
  double fout1l_a[] = {
    0.5477225575051661*fin[11]-0.3162277660168379*fin[10]-0.5477225575051661*fin[9]+0.3162277660168379*fin[8]-
    0.8215838362577489*fin[7]+0.4743416490252568*fin[6]+0.6123724356957944*fin[5]+0.8215838362577489*fin[4]-0.3535533905932737*fin[3]-
    0.4743416490252568*fin[2]-0.6123724356957944*fin[1]+0.3535533905932737*fin[0],-0.6846531968814573*fin[11]+0.3952847075210473*fin[10]+
    0.6846531968814574*fin[9]-0.3952847075210473*fin[8]+0.6123724356957944*fin[5]-0.3535533905932737*fin[3]-0.6123724356957944*fin[1]+
    0.3535533905932737*fin[0],0.5477225575051661*fin[11]-0.3162277660168379*fin[10]-0.5477225575051661*fin[9]+0.3162277660168379*fin[8]+
    0.8215838362577489*fin[7]-0.4743416490252568*fin[6]+0.6123724356957944*fin[5]-0.8215838362577489*fin[4]-0.3535533905932737*fin[3]+
    0.4743416490252568*fin[2]-0.6123724356957944*fin[1]+0.3535533905932737*fin[0],-0.5477225575051661*fin[11]+0.3162277660168379*fin[10]-
    0.5477225575051661*fin[9]+0.3162277660168379*fin[8]+0.8215838362577489*fin[7]-0.4743416490252568*fin[6]-0.6123724356957944*fin[5]+
    0.8215838362577489*fin[4]+0.3535533905932737*fin[3]-0.4743416490252568*fin[2]-0.6123724356957944*fin[1]+0.3535533905932737*fin[0],
    0.6846531968814573*fin[11]-0.3952847075210473*fin[10]+0.6846531968814574*fin[9]-0.3952847075210473*fin[8]-0.6123724356957944*fin[5]+
    0.3535533905932737*fin[3]-0.6123724356957944*fin[1]+0.3535533905932737*fin[0],-0.5477225575051661*fin[11]+0.3162277660168379*fin[10]-
    0.5477225575051661*fin[9]+0.3162277660168379*fin[8]-0.8215838362577489*fin[7]+0.4743416490252568*fin[6]-0.6123724356957944*fin[5]-
    0.8215838362577489*fin[4]+0.3535533905932737*fin[3]+0.4743416490252568*fin[2]-0.6123724356957944*fin[1]+0.3535533905932737*fin[0]
  };
  for (int i=0; i<numnod; i++) {
    doubleFunc_t sfunc = funcs1l[i];
    double fsurf = sfunc(fin);
    TEST_CHECK( gkyl_compare(fout1l_a[i], fsurf, 1e-12) );
    TEST_MSG("Expected: %.13e in i=%d", fout1l_a[i], i);
    TEST_MSG("Produced: %.13e", fsurf);
  }

  // Evaluate at right nodes on surf perp to dir 1.
  doubleFunc_t funcs1r[] = {gkhyb_1x2v_p1_surfx1_eval_quad_node_0_r,
                            gkhyb_1x2v_p1_surfx1_eval_quad_node_1_r,
                            gkhyb_1x2v_p1_surfx1_eval_quad_node_2_r,
                            gkhyb_1x2v_p1_surfx1_eval_quad_node_3_r,
                            gkhyb_1x2v_p1_surfx1_eval_quad_node_4_r,
                            gkhyb_1x2v_p1_surfx1_eval_quad_node_5_r};
  double fout1r_a[] = {
    -0.5477225575051661*fin[11]-0.3162277660168379*fin[10]+0.5477225575051661*fin[9]+0.3162277660168379*fin[8]+
    0.8215838362577489*fin[7]+0.4743416490252568*fin[6]-0.6123724356957944*fin[5]-0.8215838362577489*fin[4]-0.3535533905932737*fin[3]-
    0.4743416490252568*fin[2]+0.6123724356957944*fin[1]+0.3535533905932737*fin[0],0.6846531968814573*fin[11]+0.3952847075210473*fin[10]-
    0.6846531968814574*fin[9]-0.3952847075210473*fin[8]-0.6123724356957944*fin[5]-0.3535533905932737*fin[3]+0.6123724356957944*fin[1]+
    0.3535533905932737*fin[0],-0.5477225575051661*fin[11]-0.3162277660168379*fin[10]+0.5477225575051661*fin[9]+0.3162277660168379*fin[8]-
    0.8215838362577489*fin[7]-0.4743416490252568*fin[6]-0.6123724356957944*fin[5]+0.8215838362577489*fin[4]-0.3535533905932737*fin[3]+
    0.4743416490252568*fin[2]+0.6123724356957944*fin[1]+0.3535533905932737*fin[0],0.5477225575051661*fin[11]+0.3162277660168379*fin[10]+
    0.5477225575051661*fin[9]+0.3162277660168379*fin[8]-0.8215838362577489*fin[7]-0.4743416490252568*fin[6]+0.6123724356957944*fin[5]-
    0.8215838362577489*fin[4]+0.3535533905932737*fin[3]-0.4743416490252568*fin[2]+0.6123724356957944*fin[1]+0.3535533905932737*fin[0],-
    0.6846531968814573*fin[11]-0.3952847075210473*fin[10]-0.6846531968814574*fin[9]-0.3952847075210473*fin[8]+0.6123724356957944*fin[5]+
    0.3535533905932737*fin[3]+0.6123724356957944*fin[1]+0.3535533905932737*fin[0],0.5477225575051661*fin[11]+0.3162277660168379*fin[10]+
    0.5477225575051661*fin[9]+0.3162277660168379*fin[8]+0.8215838362577489*fin[7]+0.4743416490252568*fin[6]+0.6123724356957944*fin[5]+
    0.8215838362577489*fin[4]+0.3535533905932737*fin[3]+0.4743416490252568*fin[2]+0.6123724356957944*fin[1]+0.3535533905932737*fin[0]
  };
  for (int i=0; i<numnod; i++) {
    doubleFunc_t sfunc = funcs1r[i];
    double fsurf = sfunc(fin);
    TEST_CHECK( gkyl_compare(fout1r_a[i], fsurf, 1e-12) );
    TEST_MSG("Expected: %.13e in i=%d", fout1r_a[i], i);
    TEST_MSG("Produced: %.13e", fsurf);
  }

  // Evaluate at left nodes on surf perp to dir 2.
  doubleFunc_t funcs2l[] = {gkhyb_1x2v_p1_surfx2_eval_quad_node_0_l,
                            gkhyb_1x2v_p1_surfx2_eval_quad_node_1_l,
                            gkhyb_1x2v_p1_surfx2_eval_quad_node_2_l,
                            gkhyb_1x2v_p1_surfx2_eval_quad_node_3_l};
  numnod = sizeof(funcs2l)/sizeof(doubleFunc_t);
  double fout2l_a[] = {
    0.7905694150420947*fin[11]-0.7905694150420948*(fin[10]+fin[9])+0.7905694150420947*fin[8]-0.6123724356957944*fin[7]+
    0.6123724356957944*fin[6]+0.3535533905932737*fin[5]+0.6123724356957944*fin[4]-0.3535533905932737*fin[3]-0.6123724356957944*fin[2]-
    0.3535533905932737*fin[1]+0.3535533905932737*fin[0],-0.7905694150420947*fin[11]+0.7905694150420948*fin[10]-0.7905694150420948*fin[9]+
    0.7905694150420947*fin[8]+0.6123724356957944*fin[7]-0.6123724356957944*fin[6]-0.3535533905932737*fin[5]+0.6123724356957944*fin[4]+
    0.3535533905932737*fin[3]-0.6123724356957944*fin[2]-0.3535533905932737*fin[1]+0.3535533905932737*fin[0],-0.7905694150420947*fin[11]-
    0.7905694150420948*fin[10]+0.7905694150420948*fin[9]+0.7905694150420947*fin[8]+0.6123724356957944*(fin[7]+fin[6])-0.3535533905932737*
    fin[5]-0.6123724356957944*fin[4]-0.3535533905932737*fin[3]-0.6123724356957944*fin[2]+0.3535533905932737*(fin[1]+fin[0]),
    0.7905694150420947*fin[11]+0.7905694150420948*(fin[10]+fin[9])+0.7905694150420947*fin[8]-0.6123724356957944*(fin[7]+fin[6])+
    0.3535533905932737*fin[5]-0.6123724356957944*fin[4]+0.3535533905932737*fin[3]-0.6123724356957944*fin[2]+0.3535533905932737*
    (fin[1]+fin[0])
  };
  for (int i=0; i<numnod; i++) {
    doubleFunc_t sfunc = funcs2l[i];
    double fsurf = sfunc(fin);
    TEST_CHECK( gkyl_compare(fout2l_a[i], fsurf, 1e-12) );
    TEST_MSG("Expected: %.13e in i=%d", fout2l_a[i], i);
    TEST_MSG("Produced: %.13e", fsurf);
  }

  // Evaluate at right nodes on surf perp to dir 2.
  doubleFunc_t funcs2r[] = {gkhyb_1x2v_p1_surfx2_eval_quad_node_0_r,
                            gkhyb_1x2v_p1_surfx2_eval_quad_node_1_r,
                            gkhyb_1x2v_p1_surfx2_eval_quad_node_2_r,
                            gkhyb_1x2v_p1_surfx2_eval_quad_node_3_r};
  double fout2r_a[] = {
    0.7905694150420947*fin[11]-0.7905694150420948*(fin[10]+fin[9])+0.7905694150420947*fin[8]+0.6123724356957944*fin[7]-
    0.6123724356957944*fin[6]+0.3535533905932737*fin[5]-0.6123724356957944*fin[4]-0.3535533905932737*fin[3]+0.6123724356957944*fin[2]-
    0.3535533905932737*fin[1]+0.3535533905932737*fin[0],-0.7905694150420947*fin[11]+0.7905694150420948*fin[10]-0.7905694150420948*fin[9]+
    0.7905694150420947*fin[8]-0.6123724356957944*fin[7]+0.6123724356957944*fin[6]-0.3535533905932737*fin[5]-0.6123724356957944*fin[4]+
    0.3535533905932737*fin[3]+0.6123724356957944*fin[2]-0.3535533905932737*fin[1]+0.3535533905932737*fin[0],-0.7905694150420947*fin[11]-
    0.7905694150420948*fin[10]+0.7905694150420948*fin[9]+0.7905694150420947*fin[8]-0.6123724356957944*(fin[7]+fin[6])-0.3535533905932737*
    fin[5]+0.6123724356957944*fin[4]-0.3535533905932737*fin[3]+0.6123724356957944*fin[2]+0.3535533905932737*(fin[1]+fin[0]),
    0.7905694150420947*fin[11]+0.7905694150420948*(fin[10]+fin[9])+0.7905694150420947*fin[8]+0.6123724356957944*(fin[7]+fin[6])+
    0.3535533905932737*fin[5]+0.6123724356957944*fin[4]+0.3535533905932737*fin[3]+0.6123724356957944*fin[2]+0.3535533905932737*
    (fin[1]+fin[0])
  };
  for (int i=0; i<numnod; i++) {
    doubleFunc_t sfunc = funcs2r[i];
    double fsurf = sfunc(fin);
    TEST_CHECK( gkyl_compare(fout2r_a[i], fsurf, 1e-12) );
    TEST_MSG("Expected: %.13e in i=%d", fout2r_a[i], i);
    TEST_MSG("Produced: %.13e", fsurf);
  }

  // Evaluate at left nodes on surf perp to dir 3.
  doubleFunc_t funcs3l[] = {gkhyb_1x2v_p1_surfx3_eval_quad_node_0_l,
                            gkhyb_1x2v_p1_surfx3_eval_quad_node_1_l,
                            gkhyb_1x2v_p1_surfx3_eval_quad_node_2_l,
                            gkhyb_1x2v_p1_surfx3_eval_quad_node_3_l,
                            gkhyb_1x2v_p1_surfx3_eval_quad_node_4_l,
                            gkhyb_1x2v_p1_surfx3_eval_quad_node_5_l};
  numnod = sizeof(funcs3l)/sizeof(doubleFunc_t);
  double fout3l_a[] = {
    0.5477225575051661*fin[11]-0.5477225575051661*fin[10]-0.3162277660168379*fin[9]+0.3162277660168379*fin[8]-
    0.8215838362577489*fin[7]+0.8215838362577489*fin[6]+0.6123724356957944*fin[5]+0.4743416490252568*fin[4]-0.6123724356957944*fin[3]-
    0.4743416490252568*fin[2]-0.3535533905932737*fin[1]+0.3535533905932737*fin[0],-0.6846531968814573*fin[11]+0.6846531968814574*fin[10]+
    0.3952847075210473*fin[9]-0.3952847075210473*fin[8]+0.6123724356957944*fin[5]-0.6123724356957944*fin[3]-0.3535533905932737*fin[1]+
    0.3535533905932737*fin[0],0.5477225575051661*fin[11]-0.5477225575051661*fin[10]-0.3162277660168379*fin[9]+0.3162277660168379*fin[8]+
    0.8215838362577489*fin[7]-0.8215838362577489*fin[6]+0.6123724356957944*fin[5]-0.4743416490252568*fin[4]-0.6123724356957944*fin[3]+
    0.4743416490252568*fin[2]-0.3535533905932737*fin[1]+0.3535533905932737*fin[0],-0.5477225575051661*(fin[11]+fin[10])+
    0.3162277660168379*fin[9]+0.3162277660168379*fin[8]+0.8215838362577489*(fin[7]+fin[6])-0.6123724356957944*fin[5]-0.4743416490252568*
    fin[4]-0.6123724356957944*fin[3]-0.4743416490252568*fin[2]+0.3535533905932737*(fin[1]+fin[0]),0.6846531968814573*fin[11]+
    0.6846531968814574*fin[10]-0.3952847075210473*fin[9]-0.3952847075210473*fin[8]-0.6123724356957944*(fin[5]+fin[3])+0.3535533905932737*
    (fin[1]+fin[0]),-0.5477225575051661*(fin[11]+fin[10])+0.3162277660168379*fin[9]+0.3162277660168379*fin[8]-0.8215838362577489*
    (fin[7]+fin[6])-0.6123724356957944*fin[5]+0.4743416490252568*fin[4]-0.6123724356957944*fin[3]+0.4743416490252568*fin[2]+
    0.3535533905932737*(fin[1]+fin[0])
  };
  for (int i=0; i<numnod; i++) {
    doubleFunc_t sfunc = funcs3l[i];
    double fsurf = sfunc(fin);
    TEST_CHECK( gkyl_compare(fout3l_a[i], fsurf, 1e-12) );
    TEST_MSG("Expected: %.13e in i=%d", fout3l_a[i], i);
    TEST_MSG("Produced: %.13e", fsurf);
  }

  // Evaluate at right nodes on surf perp to dir 3.
  doubleFunc_t funcs3r[] = {gkhyb_1x2v_p1_surfx3_eval_quad_node_0_r,
                            gkhyb_1x2v_p1_surfx3_eval_quad_node_1_r,
                            gkhyb_1x2v_p1_surfx3_eval_quad_node_2_r,
                            gkhyb_1x2v_p1_surfx3_eval_quad_node_3_r,
                            gkhyb_1x2v_p1_surfx3_eval_quad_node_4_r,
                            gkhyb_1x2v_p1_surfx3_eval_quad_node_5_r};
  double fout3r_a[] = {
    -0.5477225575051661*fin[11]+0.5477225575051661*fin[10]-0.3162277660168379*fin[9]+0.3162277660168379*fin[8]+
    0.8215838362577489*fin[7]-0.8215838362577489*fin[6]-0.6123724356957944*fin[5]+0.4743416490252568*fin[4]+0.6123724356957944*fin[3]-
    0.4743416490252568*fin[2]-0.3535533905932737*fin[1]+0.3535533905932737*fin[0],0.6846531968814573*fin[11]-0.6846531968814574*fin[10]+
    0.3952847075210473*fin[9]-0.3952847075210473*fin[8]-0.6123724356957944*fin[5]+0.6123724356957944*fin[3]-0.3535533905932737*fin[1]+
    0.3535533905932737*fin[0],-0.5477225575051661*fin[11]+0.5477225575051661*fin[10]-0.3162277660168379*fin[9]+0.3162277660168379*fin[8]-
    0.8215838362577489*fin[7]+0.8215838362577489*fin[6]-0.6123724356957944*fin[5]-0.4743416490252568*fin[4]+0.6123724356957944*fin[3]+
    0.4743416490252568*fin[2]-0.3535533905932737*fin[1]+0.3535533905932737*fin[0],0.5477225575051661*(fin[11]+fin[10])+0.3162277660168379*
    fin[9]+0.3162277660168379*fin[8]-0.8215838362577489*(fin[7]+fin[6])+0.6123724356957944*fin[5]-0.4743416490252568*fin[4]+
    0.6123724356957944*fin[3]-0.4743416490252568*fin[2]+0.3535533905932737*(fin[1]+fin[0]),-0.6846531968814573*fin[11]-0.6846531968814574*
    fin[10]-0.3952847075210473*fin[9]-0.3952847075210473*fin[8]+0.6123724356957944*(fin[5]+fin[3])+0.3535533905932737*(fin[1]+fin[0]),
    0.5477225575051661*(fin[11]+fin[10])+0.3162277660168379*fin[9]+0.3162277660168379*fin[8]+0.8215838362577489*(fin[7]+fin[6])+
    0.6123724356957944*fin[5]+0.4743416490252568*fin[4]+0.6123724356957944*fin[3]+0.4743416490252568*fin[2]+0.3535533905932737*
    (fin[1]+fin[0])
  };
  for (int i=0; i<numnod; i++) {
    doubleFunc_t sfunc = funcs3r[i];
    double fsurf = sfunc(fin);
    TEST_CHECK( gkyl_compare(fout3r_a[i], fsurf, 1e-12) );
    TEST_MSG("Expected: %.13e in i=%d", fout3r_a[i], i);
    TEST_MSG("Produced: %.13e", fsurf);
  }
}

void
test_gkhyb_1x2v_upwind_quad_to_modal(struct gkyl_basis basis)
{
  // Accepted results here were generated with
  // ms-basis_surf_quad_upwind_accepted_results.mac
  double fin[basis.num_basis];
  for (int i=0; i<basis.num_basis; ++i)
    fin[i] = i+1.0;

  typedef void (*voidFunc_t)(const double* fUpwindQuad, double* GKYL_RESTRICT fUpwind);
  voidFunc_t funcs[] = {gkhyb_1x2v_p1_xdir_upwind_quad_to_modal,
                        gkhyb_1x2v_p1_vpardir_upwind_quad_to_modal,
                        gkhyb_1x2v_p1_mudir_upwind_quad_to_modal};
  int numdirs = sizeof(funcs)/sizeof(voidFunc_t);

  typedef double (*doubleFunc_t)(const double* GKYL_RESTRICT f);
  int dir, numnod;

  // Evaluate at left nodes on surf perp to dir 1.
  dir = 0;
  doubleFunc_t funcs1l[] = {gkhyb_1x2v_p1_surfx1_eval_quad_node_0_l,
                            gkhyb_1x2v_p1_surfx1_eval_quad_node_1_l,
                            gkhyb_1x2v_p1_surfx1_eval_quad_node_2_l,
                            gkhyb_1x2v_p1_surfx1_eval_quad_node_3_l,
                            gkhyb_1x2v_p1_surfx1_eval_quad_node_4_l,
                            gkhyb_1x2v_p1_surfx1_eval_quad_node_5_l};
  numnod = sizeof(funcs1l)/sizeof(doubleFunc_t);
  double f1l_n[numnod], f1l_c[numnod];
  for (int i=0; i<numnod; i++) {
    doubleFunc_t sfunc = funcs1l[i];
    f1l_n[i] = sfunc(fin);
  }
  // Transform nodal coefficients to modal coefficients.
  voidFunc_t n2mfunc1l = funcs[dir];
  n2mfunc1l(f1l_n, f1l_c);
  // Check results.
  double fout1l_a[] = {
    1.863864650467316e-16*fin[9]+1.643351536295413e-16*fin[8]-1.224744871391586*fin[1]+0.7071067811865468*fin[0],
    0.7071067811865475*fin[2]-1.224744871391588*fin[4],1.863864650467316e-16*fin[11]-1.224744871391586*fin[5]+0.7071067811865468*fin[3],
    0.7071067811865475*fin[6]-1.224744871391588*fin[7],0.7071067811865472*fin[8]-1.22474487139159*fin[9],0.7071067811865472*fin[10]-
    1.22474487139159*fin[11]
  };
  for (int i=0; i<numnod; i++) {
    TEST_CHECK( gkyl_compare(fout1l_a[i], f1l_c[i], 1e-12) );
    TEST_MSG("Expected: %.13e in i=%d", fout1l_a[i], i);
    TEST_MSG("Produced: %.13e", f1l_c[i]);
  }

  // Evaluate at left nodes on surf perp to dir 2.
  dir = 1;
  doubleFunc_t funcs2l[] = {gkhyb_1x2v_p1_surfx2_eval_quad_node_0_l,
                            gkhyb_1x2v_p1_surfx2_eval_quad_node_1_l,
                            gkhyb_1x2v_p1_surfx2_eval_quad_node_2_l,
                            gkhyb_1x2v_p1_surfx2_eval_quad_node_3_l};
  numnod = sizeof(funcs2l)/sizeof(doubleFunc_t);
  double f2l_n[numnod], f2l_c[numnod];
  for (int i=0; i<numnod; i++) {
    doubleFunc_t sfunc = funcs2l[i];
    f2l_n[i] = sfunc(fin);
  }
  // Transform nodal coefficients to modal coefficients.
  voidFunc_t n2mfunc2l = funcs[dir];
  n2mfunc2l(f2l_n, f2l_c);
  // Check results.
  double fout2l_a[] = {
    1.581138830084189*fin[8]-1.224744871391586*fin[2]+0.7071067811865468*fin[0],1.581138830084189*fin[9]-1.224744871391586*fin[4]
    +0.7071067811865468*fin[1],1.581138830084189*fin[10]-1.224744871391586*fin[6]+0.7071067811865468*fin[3],1.581138830084189*fin[11]-
    1.224744871391586*fin[7]+0.7071067811865468*fin[5]
  };
  for (int i=0; i<numnod; i++) {
    TEST_CHECK( gkyl_compare(fout2l_a[i], f2l_c[i], 1e-12) );
    TEST_MSG("Expected: %.13e in i=%d", fout2l_a[i], i);
    TEST_MSG("Produced: %.13e", f2l_c[i]);
  }

  // Evaluate at left nodes on surf perp to dir 3.
  dir = 2;
  doubleFunc_t funcs3l[] = {gkhyb_1x2v_p1_surfx3_eval_quad_node_0_l,
                            gkhyb_1x2v_p1_surfx3_eval_quad_node_1_l,
                            gkhyb_1x2v_p1_surfx3_eval_quad_node_2_l,
                            gkhyb_1x2v_p1_surfx3_eval_quad_node_3_l,
                            gkhyb_1x2v_p1_surfx3_eval_quad_node_4_l,
                            gkhyb_1x2v_p1_surfx3_eval_quad_node_5_l};
  numnod = sizeof(funcs3l)/sizeof(doubleFunc_t);
  double f3l_n[numnod], f3l_c[numnod];
  for (int i=0; i<numnod; i++) {
    doubleFunc_t sfunc = funcs3l[i];
    f3l_n[i] = sfunc(fin);
  }
  // Transform nodal coefficients to modal coefficients.
  voidFunc_t n2mfunc3l = funcs[dir];
  n2mfunc3l(f3l_n, f3l_c);
  // Check results.
  double fout3l_a[] = {
    1.863864650467316e-16*fin[10]+1.643351536295413e-16*fin[8]-1.224744871391586*fin[3]+0.7071067811865468*fin[0],
    1.863864650467316e-16*fin[11]-1.224744871391586*fin[5]+0.7071067811865468*fin[1],0.7071067811865475*fin[2]-1.224744871391588*fin[6],
    0.7071067811865475*fin[4]-1.224744871391588*fin[7],0.7071067811865472*fin[8]-1.22474487139159*fin[10],0.7071067811865472*fin[9]-
    1.22474487139159*fin[11]
  };
  for (int i=0; i<numnod; i++) {
    TEST_CHECK( gkyl_compare(fout3l_a[i], f3l_c[i], 1e-12) );
    TEST_MSG("Expected: %.13e in i=%d", fout3l_a[i], i);
    TEST_MSG("Produced: %.13e", f3l_c[i]);
  }
}

void
test_gkhyb()
{
  struct gkyl_basis basis1;
  gkyl_cart_modal_gkhybrid(&basis1, 1, 2);
  test_gkhyb_members(basis1);

  struct gkyl_basis *basis2 = gkyl_cart_modal_gkhybrid_new(1, 2);
  test_gkhyb_members(*basis2);
  // Test surf quad node evaluation and modal-to-nodal transform.
  test_gkhyb_1x2v_surf_eval_nod(*basis2);
  test_gkhyb_1x2v_upwind_quad_to_modal(*basis2);
  gkyl_cart_modal_basis_release(basis2);
}

#ifdef GKYL_HAVE_CUDA

int dev_cu_ser_2d(struct gkyl_basis *basis);

void
test_cu_ser_2d_members(struct gkyl_basis *basis)
{
  int nfail = dev_cu_ser_2d(basis);

  TEST_CHECK(nfail == 0);
}

void
test_cu_ser_2d()
{
  struct gkyl_basis *basis1 = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  gkyl_cart_modal_serendip_cu_dev(basis1, 2, 2);
  test_cu_ser_2d_members(basis1);
  gkyl_cart_modal_basis_release_cu(basis1);

  struct gkyl_basis *basis2 = gkyl_cart_modal_serendip_cu_dev_new(2, 2);
  test_cu_ser_2d_members(basis2);
  gkyl_cart_modal_basis_release_cu(basis2);
}
#endif

TEST_LIST = {
  { "ser_1d", test_ser_1d },
  { "ser_2d", test_ser_2d },
  { "ten_2d", test_ten_2d },
  { "hyb", test_hyb },
  { "gkhyb", test_gkhyb },
#ifdef GKYL_HAVE_CUDA
  { "cu_ser_2d", test_cu_ser_2d },
#endif    
  { NULL, NULL },
};
