#include <acutest.h>

#include <gkyl_mat_triples.h>
#include <gkyl_superlu.h>
#include <gkyl_superlu_ops.h>
#include <gkyl_util.h>

#include <stdbool.h>

void test_cusolver_qr_dev();
void test_cusolver_rf_dev();
void test_cusolver_ops_dev();
void test_cusolver_ops_multiple_rhs_dev();
void test_cusolver_ops_multiple_prob_dev();
void test_cudss_simple_dev();
void test_cudss_ops_dev();
void test_cudss_ops_update_amat_dev();
void test_cudss_ops_multiple_rhs_dev();

void test_slu_example_ho()
{
/*  
 * This is the small 5x5 example used in the Sections 2 and 3 of the 
 * Users' Guide to illustrate how to call a SuperLU routine, and the
 * matrix data structures used by SuperLU.
 *
 */
  SuperMatrix A, L, U, B;
  double   *a, *rhs;
  double   s, u, p, e, r, l;
  int      *asub, *xa;
  int      *perm_r; /* row permutations from partial pivoting */
  int      *perm_c; /* column permutation vector */
  int      nrhs, info, i, m, n, nnz, permc_spec;
  superlu_options_t options;
  SuperLUStat_t stat;

  /* Initialize matrix A. */
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  m = n = 5;
  nnz = 12;
  if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a[].");
  if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub[].");
  if ( !(xa = intMalloc(n+1)) ) ABORT("Malloc fails for xa[].");
  
  s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  a[0] = s; a[1] = l; a[2] = l; a[3] = u; a[4] = l; a[5] = l;
  a[6] = u; a[7] = p; a[8] = u; a[9] = e; a[10]= u; a[11]= r;
  asub[0] = 0; asub[1] = 1; asub[2] = 4; asub[3] = 1;
  asub[4] = 2; asub[5] = 4; asub[6] = 0; asub[7] = 2;
  asub[8] = 0; asub[9] = 3; asub[10]= 3; asub[11]= 4;
  xa[0] = 0; xa[1] = 3; xa[2] = 6; xa[3] = 8; xa[4] = 10; xa[5] = nnz;

  /* Create matrix A in the format expected by SuperLU. */
  dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    
  /* Create right-hand side matrix B. */
  nrhs = 1;
  if ( !(rhs = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");
  /* B : transpose([1,1,1,1,1]);*/
  for (i = 0; i < m; ++i) rhs[i] = 1.0;
  dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);

  if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
  if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");

  /* Set the default input options. */
  set_default_options(&options);
  options.ColPerm = NATURAL;

  /* Initialize the statistics variables. */
  StatInit(&stat);

  /* Solve linear system */
  dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

  /* Solution is: [-1/32, 11/168, 3/224, 1/16, 11/336] */
  TEST_CHECK( gkyl_compare(-1.0/32.0, rhs[0], 1e-14) );
  TEST_CHECK( gkyl_compare( 11.0/168.0, rhs[1], 1e-14) );
  TEST_CHECK( gkyl_compare( 3.0/224.0, rhs[2], 1e-14) );
  TEST_CHECK( gkyl_compare( 1.0/16.0, rhs[3], 1e-14) );
  TEST_CHECK( gkyl_compare( 11.0/336.0, rhs[4], 1e-14) );

  /* De-allocate storage */
  SUPERLU_FREE (rhs);
  SUPERLU_FREE (perm_r);
  SUPERLU_FREE (perm_c);
  Destroy_CompCol_Matrix(&A);
  Destroy_SuperMatrix_Store(&B);
  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);
  StatFree(&stat);
}

void test_superlu_ops(const bool separateLUdecomp)
{
  /*  
   * Like test_slu_example_ho but using superlu_ops.
   *
   */
  double s, u, p, e, r, l;
  int    nrhs, m, n, nnz;

  /* Initialize matrix A. */
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  m = n = 5;
  nrhs = 1;
  
  s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  struct gkyl_mat_triples **tri_arr = gkyl_malloc(sizeof(struct gkyl_mat_triples *));
  tri_arr[0] = gkyl_mat_triples_new(m, n);
  struct gkyl_mat_triples *tri = tri_arr[0];
  // row 0
  gkyl_mat_triples_insert(tri, 0, 0, s);
  gkyl_mat_triples_insert(tri, 0, 2, u);
  gkyl_mat_triples_insert(tri, 0, 3, u);
  // row 1
  gkyl_mat_triples_insert(tri, 1, 0, l);
  gkyl_mat_triples_insert(tri, 1, 1, u);
  // row 2
  gkyl_mat_triples_insert(tri, 2, 1, l);
  gkyl_mat_triples_insert(tri, 2, 2, p);
  // row 3
  gkyl_mat_triples_insert(tri, 3, 3, e);
  gkyl_mat_triples_insert(tri, 3, 4, u);
  // row 4
  gkyl_mat_triples_insert(tri, 4, 0, l);
  gkyl_mat_triples_insert(tri, 4, 1, l);
  gkyl_mat_triples_insert(tri, 4, 4, r);

  // Create the SuperLU linear problem setup.
  gkyl_superlu_prob *sluprob = gkyl_superlu_prob_new(1, m, n, nrhs);

  // Allocate the A matrix from triples.
  gkyl_superlu_amat_from_triples(sluprob, tri_arr);
  gkyl_mat_triples_release(tri_arr[0]);
  gkyl_free(tri_arr);

  // Create right-hand side matrix B = transpose([1,1,1,1,1]).
  gkyl_mat_triples *triRHS = gkyl_mat_triples_new(m, nrhs);
  gkyl_mat_triples_insert(triRHS, 0, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 1, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 2, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 3, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 4, 0, 1.0);
  gkyl_superlu_brhs_from_triples(sluprob, triRHS);
  gkyl_mat_triples_release(triRHS);

  // Solve linear system.
  if (separateLUdecomp) {
    gkyl_superlu_ludecomp(sluprob);
    gkyl_superlu_solve(sluprob);
  } else {
    gkyl_superlu_solve(sluprob);
  }

  // Solution is: [-1/32, 11/168, 3/224, 1/16, 11/336].
  TEST_CHECK( gkyl_compare(-1.0/32.0,   gkyl_superlu_get_rhs_lin(sluprob,0), 1e-14) );
  TEST_CHECK( gkyl_compare( 11.0/168.0, gkyl_superlu_get_rhs_lin(sluprob,1), 1e-14) );
  TEST_CHECK( gkyl_compare( 3.0/224.0,  gkyl_superlu_get_rhs_lin(sluprob,2), 1e-14) );
  TEST_CHECK( gkyl_compare( 1.0/16.0,   gkyl_superlu_get_rhs_lin(sluprob,3), 1e-14) );
  TEST_CHECK( gkyl_compare( 11.0/336.0, gkyl_superlu_get_rhs_lin(sluprob,4), 1e-14) );

  gkyl_superlu_prob_release(sluprob);

}

void test_superlu_ops_basic_ho() {
  test_superlu_ops(false);
}

void test_superlu_ops_separateLU_ho() {
  test_superlu_ops(true);
}

void test_superlu_ops_basic_update_amat_ho() {
  /*  
   * Like test_slu_example_ho but using superlu_ops.
   *
   */
  double s, u, p, e, r, l;
  int    nrhs, m, n, nnz;

  /* Initialize matrix A. */
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  m = n = 5;
  nrhs = 1;
  
  s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  struct gkyl_mat_triples **tri_arr = gkyl_malloc(sizeof(struct gkyl_mat_triples *));
  tri_arr[0] = gkyl_mat_triples_new(m, n);
  struct gkyl_mat_triples *tri = tri_arr[0];
  // row 0
  gkyl_mat_triples_insert(tri, 0, 0, s);
  gkyl_mat_triples_insert(tri, 0, 2, u);
  gkyl_mat_triples_insert(tri, 0, 3, u);
  // row 1
  gkyl_mat_triples_insert(tri, 1, 0, l);
  gkyl_mat_triples_insert(tri, 1, 1, u);
  // row 2
  gkyl_mat_triples_insert(tri, 2, 1, l);
  gkyl_mat_triples_insert(tri, 2, 2, p);
  // row 3
  gkyl_mat_triples_insert(tri, 3, 3, e);
  gkyl_mat_triples_insert(tri, 3, 4, u);
  // row 4
  gkyl_mat_triples_insert(tri, 4, 0, l);
  gkyl_mat_triples_insert(tri, 4, 1, l);
  gkyl_mat_triples_insert(tri, 4, 4, r);

  // Create the SuperLU linear problem setup.
  gkyl_superlu_prob *sluprob = gkyl_superlu_prob_new(1, m, n, nrhs);

  // Allocate the A matrix from triples.
  gkyl_superlu_amat_from_triples(sluprob, tri_arr);

  // Create right-hand side matrix B = transpose([1,1,1,1,1]).
  gkyl_mat_triples *triRHS = gkyl_mat_triples_new(m, nrhs);
  gkyl_mat_triples_insert(triRHS, 0, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 1, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 2, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 3, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 4, 0, 1.0);
  gkyl_superlu_brhs_from_triples(sluprob, triRHS);

  // Solve linear system.
  gkyl_superlu_solve(sluprob);

  // Solution is: [-1/32, 11/168, 3/224, 1/16, 11/336].
  TEST_CHECK( gkyl_compare(-1.0/32.0  , gkyl_superlu_get_rhs_lin(sluprob,0), 1e-14) );
  TEST_CHECK( gkyl_compare( 11.0/168.0, gkyl_superlu_get_rhs_lin(sluprob,1), 1e-14) );
  TEST_CHECK( gkyl_compare( 3.0/224.0 , gkyl_superlu_get_rhs_lin(sluprob,2), 1e-14) );
  TEST_CHECK( gkyl_compare( 1.0/16.0  , gkyl_superlu_get_rhs_lin(sluprob,3), 1e-14) );
  TEST_CHECK( gkyl_compare( 11.0/336.0, gkyl_superlu_get_rhs_lin(sluprob,4), 1e-14) );

  // Now update the LHS matrix. Multiply it by a constant so the solution should be the same as before but divided by that constant.
  double prob_fac = 1.3;

  // row 0
  gkyl_mat_triples_insert(tri, 0, 0, s*prob_fac);
  gkyl_mat_triples_insert(tri, 0, 2, u*prob_fac);
  gkyl_mat_triples_insert(tri, 0, 3, u*prob_fac);
  // row 1
  gkyl_mat_triples_insert(tri, 1, 0, l*prob_fac);
  gkyl_mat_triples_insert(tri, 1, 1, u*prob_fac);
  // row 2
  gkyl_mat_triples_insert(tri, 2, 1, l*prob_fac);
  gkyl_mat_triples_insert(tri, 2, 2, p*prob_fac);
  // row 3
  gkyl_mat_triples_insert(tri, 3, 3, e*prob_fac);
  gkyl_mat_triples_insert(tri, 3, 4, u*prob_fac);
  // row 4
  gkyl_mat_triples_insert(tri, 4, 0, l*prob_fac);
  gkyl_mat_triples_insert(tri, 4, 1, l*prob_fac);
  gkyl_mat_triples_insert(tri, 4, 4, r*prob_fac);

  gkyl_superlu_amat_update_from_triples(sluprob, tri_arr);

  // Reset RHS.
  gkyl_superlu_brhs_from_triples(sluprob, triRHS);

  // Solve linear system.
  gkyl_superlu_solve(sluprob);

  // Solution is: (1/prob_fac)*[-1/32, 11/168, 3/224, 1/16, 11/336].
  TEST_CHECK( gkyl_compare((1.0/prob_fac)*(-1.0/32.0  ), gkyl_superlu_get_rhs_lin(sluprob,0), 1e-14) );
  TEST_CHECK( gkyl_compare((1.0/prob_fac)*( 11.0/168.0), gkyl_superlu_get_rhs_lin(sluprob,1), 1e-14) );
  TEST_CHECK( gkyl_compare((1.0/prob_fac)*( 3.0/224.0 ), gkyl_superlu_get_rhs_lin(sluprob,2), 1e-14) );
  TEST_CHECK( gkyl_compare((1.0/prob_fac)*( 1.0/16.0  ), gkyl_superlu_get_rhs_lin(sluprob,3), 1e-14) );
  TEST_CHECK( gkyl_compare((1.0/prob_fac)*( 11.0/336.0), gkyl_superlu_get_rhs_lin(sluprob,4), 1e-14) );

  gkyl_mat_triples_release(tri_arr[0]);
  gkyl_free(tri_arr);
  gkyl_mat_triples_release(triRHS);
  gkyl_superlu_prob_release(sluprob);

}

double superlu_test_answer(double s, double u, double p, double e, double r, double l, int idx) {
  // Solution is: [-1/32, 11/168, 3/224, 1/16, 11/336].
  // for a unit RHS vector and
  //  s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
  double sol;
  switch (idx) {
    case 0:
      sol = (e*l*r + e*p*r - l*p*u - e*r*u - p*r*u + p*pow(u,2))/
        (e*pow(l,2)*r + e*p*r*s - pow(l,2)*p*u + l*p*pow(u,2));
      break;
    case 1:
      sol = (r*(-(e*l*p) + e*p*s + e*l*u + l*p*u))/
        (u*(e*pow(l,2)*r + e*p*r*s - pow(l,2)*p*u + l*p*pow(u,2)));
      break;
    case 2:
      sol = -((-(e*pow(l,2)*r) + e*l*r*s + pow(l,2)*r*u - e*r*s*u + pow(l,2)*pow(u,2) -
        l*pow(u,3))/(u*(e*pow(l,2)*r + e*p*r*s - pow(l,2)*p*u + l*p*pow(u,2))));
      break;
    case 3:
      sol = (-(pow(l,2)*p) + pow(l,2)*r + l*p*s + p*r*s + pow(l,2)*u + l*p*u - p*s*u -
      l*pow(u,2))/(e*pow(l,2)*r + e*p*r*s - pow(l,2)*p*u + l*p*pow(u,2));
      break;
    case 4:
      sol = (e*pow(l,2)*p - e*l*p*s - e*pow(l,2)*u - e*l*p*u - pow(l,2)*p*u + e*p*s*u +
        e*l*pow(u,2) + l*p*pow(u,2))/
        (u*(e*pow(l,2)*r + e*p*r*s - pow(l,2)*p*u + l*p*pow(u,2)));
      break;
  }
  return sol;
};

void test_superlu_ops_multiple_prob_ho()
{
  double s, u, p, e, r, l;
  int    nprob, m, n;

  /* Initialize matrix A. */
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  m = n = 5;
  nprob = 7;

  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  struct gkyl_mat_triples **tri_arr = (struct gkyl_mat_triples **) gkyl_malloc(nprob*sizeof(struct gkyl_mat_triples *));
  for (size_t k=0; k<nprob; k++) {
    tri_arr[k] = gkyl_mat_triples_new(m, n);
    struct gkyl_mat_triples *tri = tri_arr[k];

    s = 19.0*(k+1)/nprob; u = 21.0*(k+1)/nprob; p = 16.0*(k+1)/nprob; e = 5.0*(k+1)/nprob; r = 18.0*(k+1)/nprob; l = 12.0*(k+1)/nprob;

    // row 0
    gkyl_mat_triples_insert(tri, 0, 0, s);
    gkyl_mat_triples_insert(tri, 0, 2, u);
    gkyl_mat_triples_insert(tri, 0, 3, u);
    // row 1
    gkyl_mat_triples_insert(tri, 1, 0, l);
    gkyl_mat_triples_insert(tri, 1, 1, u);
    // row 2
    gkyl_mat_triples_insert(tri, 2, 1, l);
    gkyl_mat_triples_insert(tri, 2, 2, p);
    // row 3
    gkyl_mat_triples_insert(tri, 3, 3, e);
    gkyl_mat_triples_insert(tri, 3, 4, u);
    // row 4
    gkyl_mat_triples_insert(tri, 4, 0, l);
    gkyl_mat_triples_insert(tri, 4, 1, l);
    gkyl_mat_triples_insert(tri, 4, 4, r);
  }

  // Create the SuperLU linear problem setup.
  gkyl_superlu_prob *prob = gkyl_superlu_prob_new(nprob, m, n, 1);

  // Allocate the A matrix from triples.
  gkyl_superlu_amat_from_triples(prob, tri_arr);
  for (size_t k=0; k<nprob; k++)
    gkyl_mat_triples_release(tri_arr[k]);
  gkyl_free(tri_arr);

  // Create right-hand side matrix B = transpose([1,1,1,1,1]).
  gkyl_mat_triples *triRHS = gkyl_mat_triples_new(m, nprob);
  for (int k=0; k<nprob; k++) {
    gkyl_mat_triples_insert(triRHS, 0, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 1, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 2, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 3, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 4, k, 1.0);
  }
  gkyl_superlu_brhs_from_triples(prob, triRHS);
  gkyl_mat_triples_release(triRHS);

  gkyl_superlu_solve(prob);

  for (int k=0; k<nprob; k++) {
    s = 19.0*(k+1)/nprob; u = 21.0*(k+1)/nprob; p = 16.0*(k+1)/nprob; e = 5.0*(k+1)/nprob; r = 18.0*(k+1)/nprob; l = 12.0*(k+1)/nprob;
    for (int i=0; i<m; i++)
      TEST_CHECK( gkyl_compare_double( superlu_test_answer(s,u,p,e,r,l,i), gkyl_superlu_get_rhs_lin(prob,k*m+i), 1e-10) );
  }

  gkyl_superlu_prob_release(prob);
}

void test_superlu_ops_multiple_prob_update_amat_ho()
{
  double s, u, p, e, r, l;
  int    nprob, m, n;

  /* Initialize matrix A. */
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  m = n = 5;
  nprob = 7;

  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  struct gkyl_mat_triples **tri_arr = (struct gkyl_mat_triples **) gkyl_malloc(nprob*sizeof(struct gkyl_mat_triples *));
  for (size_t k=0; k<nprob; k++) {
    tri_arr[k] = gkyl_mat_triples_new(m, n);
    struct gkyl_mat_triples *tri = tri_arr[k];

    s = 19.0*(k+1)/nprob; u = 21.0*(k+1)/nprob; p = 16.0*(k+1)/nprob; e = 5.0*(k+1)/nprob; r = 18.0*(k+1)/nprob; l = 12.0*(k+1)/nprob;

    // row 0
    gkyl_mat_triples_insert(tri, 0, 0, s);
    gkyl_mat_triples_insert(tri, 0, 2, u);
    gkyl_mat_triples_insert(tri, 0, 3, u);
    // row 1
    gkyl_mat_triples_insert(tri, 1, 0, l);
    gkyl_mat_triples_insert(tri, 1, 1, u);
    // row 2
    gkyl_mat_triples_insert(tri, 2, 1, l);
    gkyl_mat_triples_insert(tri, 2, 2, p);
    // row 3
    gkyl_mat_triples_insert(tri, 3, 3, e);
    gkyl_mat_triples_insert(tri, 3, 4, u);
    // row 4
    gkyl_mat_triples_insert(tri, 4, 0, l);
    gkyl_mat_triples_insert(tri, 4, 1, l);
    gkyl_mat_triples_insert(tri, 4, 4, r);
  }

  // Create the SuperLU linear problem setup.
  gkyl_superlu_prob *prob = gkyl_superlu_prob_new(nprob, m, n, 1);

  // Allocate the A matrix from triples.
  gkyl_superlu_amat_from_triples(prob, tri_arr);

  // Create right-hand side matrix B = transpose([1,1,1,1,1]).
  gkyl_mat_triples *triRHS = gkyl_mat_triples_new(m, nprob);
  for (int k=0; k<nprob; k++) {
    gkyl_mat_triples_insert(triRHS, 0, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 1, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 2, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 3, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 4, k, 1.0);
  }
  gkyl_superlu_brhs_from_triples(prob, triRHS);

  gkyl_superlu_solve(prob);

  for (int k=0; k<nprob; k++) {
    s = 19.0*(k+1)/nprob; u = 21.0*(k+1)/nprob; p = 16.0*(k+1)/nprob; e = 5.0*(k+1)/nprob; r = 18.0*(k+1)/nprob; l = 12.0*(k+1)/nprob;
    for (int i=0; i<m; i++)
      TEST_CHECK( gkyl_compare_double( superlu_test_answer(s,u,p,e,r,l,i), gkyl_superlu_get_rhs_lin(prob,k*m+i), 1e-10) );
  }

  // Now update the LHS matrix. Multiply it by a constant, and multiply the RHS
  // by 2X that constant, so the solution should be the same but divided by 2.
  double prob_fac = 1.3;
  double *prob_fac_per_rhs = (double *) gkyl_malloc(nprob*sizeof(double));
  for (int k=0; k<nprob; k++)
    prob_fac_per_rhs[k] = (k+1)*prob_fac;

  for (size_t k=0; k<nprob; k++) {
    struct gkyl_mat_triples *tri = tri_arr[k];

    s = 19.0*(k+1)/nprob; u = 21.0*(k+1)/nprob; p = 16.0*(k+1)/nprob; e = 5.0*(k+1)/nprob; r = 18.0*(k+1)/nprob; l = 12.0*(k+1)/nprob;

    // row 0
    gkyl_mat_triples_insert(tri, 0, 0, s*prob_fac_per_rhs[k]);
    gkyl_mat_triples_insert(tri, 0, 2, u*prob_fac_per_rhs[k]);
    gkyl_mat_triples_insert(tri, 0, 3, u*prob_fac_per_rhs[k]);
    // row 1
    gkyl_mat_triples_insert(tri, 1, 0, l*prob_fac_per_rhs[k]);
    gkyl_mat_triples_insert(tri, 1, 1, u*prob_fac_per_rhs[k]);
    // row 2
    gkyl_mat_triples_insert(tri, 2, 1, l*prob_fac_per_rhs[k]);
    gkyl_mat_triples_insert(tri, 2, 2, p*prob_fac_per_rhs[k]);
    // row 3
    gkyl_mat_triples_insert(tri, 3, 3, e*prob_fac_per_rhs[k]);
    gkyl_mat_triples_insert(tri, 3, 4, u*prob_fac_per_rhs[k]);
    // row 4
    gkyl_mat_triples_insert(tri, 4, 0, l*prob_fac_per_rhs[k]);
    gkyl_mat_triples_insert(tri, 4, 1, l*prob_fac_per_rhs[k]);
    gkyl_mat_triples_insert(tri, 4, 4, r*prob_fac_per_rhs[k]);
  }

  // Create the cuSolver linear problem setup.
  gkyl_superlu_amat_update_from_triples(prob, tri_arr);

  // Reset RHS.
  for (int k=0; k<nprob; k++) {
    gkyl_mat_triples_insert(triRHS, 0, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 1, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 2, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 3, k, 1.0);
    gkyl_mat_triples_insert(triRHS, 4, k, 1.0);
  }
  gkyl_superlu_brhs_from_triples(prob, triRHS);

  gkyl_superlu_solve(prob);

  for (int k=0; k<nprob; k++) {
    double fac = prob_fac_per_rhs[k];
    s = fac*19.0*(k+1)/nprob;
    u = fac*21.0*(k+1)/nprob;
    p = fac*16.0*(k+1)/nprob;
    e = fac*5.0 *(k+1)/nprob;
    r = fac*18.0*(k+1)/nprob;
    l = fac*12.0*(k+1)/nprob;
    for (int i=0; i<m; i++)
      TEST_CHECK( gkyl_compare_double( superlu_test_answer(s,u,p,e,r,l,i), gkyl_superlu_get_rhs_lin(prob,k*m+i), 1e-10) );
  }

  gkyl_free(prob_fac_per_rhs);
  for (size_t k=0; k<nprob; k++)
    gkyl_mat_triples_release(tri_arr[k]);

  gkyl_free(tri_arr);
  gkyl_mat_triples_release(triRHS);
  gkyl_superlu_prob_release(prob);
}

TEST_LIST = {
  { "slu_example_ho", test_slu_example_ho },
  { "superlu_ops_basic_ho", test_superlu_ops_basic_ho },
  { "superlu_ops_basic_update_amat_ho", test_superlu_ops_basic_update_amat_ho },
  { "superlu_ops_separateLU_ho", test_superlu_ops_separateLU_ho },
  { "superlu_ops_multiple_prob_ho", test_superlu_ops_multiple_prob_ho },
  { "superlu_ops_multiple_prob_update_amat_ho", test_superlu_ops_multiple_prob_update_amat_ho },
#ifdef GKYL_HAVE_CUDA
#ifdef GKYL_HAVE_CUDSS
  { "cudss_simple_dev", test_cudss_simple_dev },
  { "cudss_ops_dev", test_cudss_ops_dev },
  { "cudss_ops_update_amat_dev", test_cudss_ops_update_amat_dev },
  { "cudss_ops_multiple_rhs_dev", test_cudss_ops_multiple_rhs_dev },
#else
  { "cusolver_qr_dev", test_cusolver_qr_dev },
  { "cusolver_rf_dev", test_cusolver_rf_dev },
  { "cusolver_ops_dev", test_cusolver_ops_dev },
  { "cusolver_ops_multiple_rhs_dev", test_cusolver_ops_multiple_rhs_dev },
  { "cusolver_ops_multiple_prob_dev", test_cusolver_ops_multiple_prob_dev },
#endif
#endif
  { NULL, NULL }
};
