// Additional unit tests for gkyl_basis (modal serendipity/tensor bases)
#include <acutest.h>
#include <gkyl_basis.h>
#include <gkyl_util.h>
#include <math.h>

static int ipow(int b, int e) { int r=1; for (int i=0;i<e;++i) r*=b; return r; }

void
test_basis_num_basis_serendip()
{
  struct gkyl_basis b;

  // 1D serendipity: num_basis = poly_order+1
  gkyl_cart_modal_serendip(&b, 1, 1); TEST_CHECK( b.num_basis == 2 );
  gkyl_cart_modal_serendip(&b, 1, 2); TEST_CHECK( b.num_basis == 3 );
  gkyl_cart_modal_serendip(&b, 1, 3); TEST_CHECK( b.num_basis == 4 );

  // 2D serendipity: p1=4, p2=8, p3=12
  gkyl_cart_modal_serendip(&b, 2, 1); TEST_CHECK( b.num_basis == 4 );
  gkyl_cart_modal_serendip(&b, 2, 2); TEST_CHECK( b.num_basis == 8 );
  gkyl_cart_modal_serendip(&b, 2, 3); TEST_CHECK( b.num_basis == 12 );

  // 3D serendipity: p1=8, p2=20, p3=32
  gkyl_cart_modal_serendip(&b, 3, 1); TEST_CHECK( b.num_basis == 8 );
  gkyl_cart_modal_serendip(&b, 3, 2); TEST_CHECK( b.num_basis == 20 );

  // metadata sanity
  gkyl_cart_modal_serendip(&b, 2, 2);
  TEST_CHECK( b.ndim == 2 );
  TEST_CHECK( b.poly_order == 2 );
  TEST_CHECK( b.b_type == GKYL_BASIS_MODAL_SERENDIPITY );
}

void
test_basis_num_basis_tensor()
{
  struct gkyl_basis b;
  // tensor basis: num_basis = (poly_order+1)^ndim
  for (int ndim=1; ndim<=3; ++ndim) {
    for (int p=1; p<=2; ++p) {
      gkyl_cart_modal_tensor(&b, ndim, p);
      TEST_CHECK( b.num_basis == (unsigned) ipow(p+1, ndim) );
      TEST_CHECK( b.b_type == GKYL_BASIS_MODAL_TENSOR );
    }
  }
}

void
test_basis_eval_expand_consistency()
{
  // eval_expand(z,f) must equal sum_i f[i]*b_i(z), where b = eval(z).
  struct gkyl_basis b;
  gkyl_cart_modal_serendip(&b, 2, 2);

  double z[2] = { 0.3, -0.7 };
  double bz[8];
  b.eval(z, bz);

  double f[8] = { 1.0, 2.0, -1.5, 0.5, 3.0, -2.0, 0.25, -0.75 };
  double ev = b.eval_expand(z, f);

  double sum = 0.0;
  for (int i=0; i<8; ++i) sum += f[i]*bz[i];
  TEST_CHECK( gkyl_compare_double(ev, sum, 1e-13) );
}

void
test_basis_const_mode()
{
  // The first (0th) basis function is constant: same value at any z, and >0.
  struct gkyl_basis b;
  gkyl_cart_modal_serendip(&b, 1, 2);

  double z1[1] = { -0.9 }, z2[1] = { 0.4 };
  double b1[3], b2[3];
  b.eval(z1, b1);
  b.eval(z2, b2);
  TEST_CHECK( gkyl_compare_double(b1[0], b2[0], 1e-14) );
  TEST_CHECK( b1[0] > 0.0 );

  // Pure constant expansion evaluates to the same value everywhere.
  double f[3] = { 1.0, 0.0, 0.0 };
  TEST_CHECK( gkyl_compare_double(b.eval_expand(z1, f), b.eval_expand(z2, f), 1e-14) );
}

void
test_basis_flip_odd_involution()
{
  // Applying flip_odd_sign twice in the same direction recovers the input.
  struct gkyl_basis b;
  gkyl_cart_modal_serendip(&b, 2, 2);

  double f[8] = { 1.0, 2.0, -1.5, 0.5, 3.0, -2.0, 0.25, -0.75 };
  double f1[8], f2[8];
  b.flip_odd_sign(0, f, f1);
  b.flip_odd_sign(0, f1, f2);
  for (int i=0; i<8; ++i)
    TEST_CHECK( gkyl_compare_double(f2[i], f[i], 1e-14) );
}

void
test_basis_flip_odd_reflection()
{
  // flip_odd_sign(dir) on coefficients corresponds to reflecting z->-z in dir:
  // expand(flip(f))(z) == expand(f)(z with z[dir] negated).
  struct gkyl_basis b;
  gkyl_cart_modal_serendip(&b, 1, 2);

  double f[3] = { 1.5, -0.5, 0.8 };
  double ff[3];
  b.flip_odd_sign(0, f, ff);

  double z[1] = { 0.6 }, zm[1] = { -0.6 };
  TEST_CHECK( gkyl_compare_double(b.eval_expand(z, ff), b.eval_expand(zm, f), 1e-13) );
}

void
test_basis_new_release()
{
  struct gkyl_basis *b = gkyl_cart_modal_serendip_new(3, 1);
  TEST_CHECK( b->num_basis == 8 );
  TEST_CHECK( b->ndim == 3 );
  gkyl_cart_modal_basis_release(b);

  struct gkyl_basis *t = gkyl_cart_modal_tensor_new(2, 2);
  TEST_CHECK( t->num_basis == 9 );
  gkyl_cart_modal_basis_release(t);
}

TEST_LIST = {
  { "num_basis_serendip", test_basis_num_basis_serendip },
  { "num_basis_tensor", test_basis_num_basis_tensor },
  { "eval_expand_consistency", test_basis_eval_expand_consistency },
  { "const_mode", test_basis_const_mode },
  { "flip_odd_involution", test_basis_flip_odd_involution },
  { "flip_odd_reflection", test_basis_flip_odd_reflection },
  { "new_release", test_basis_new_release },
  { NULL, NULL },
};
