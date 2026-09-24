#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_hyper_dg.h>

static struct gkyl_array*
mkarr1(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (use_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

//#define mkarr1(use_gpu, nc, size) (fprintf(stderr, "mkarr1: %d\n", __LINE__), mkarr1_(use_gpu, nc, size))

int hyper_dg_kernel_test(const gkyl_hyper_dg *slvr);

void
test_vlasov_1x2v_p2_(bool use_gpu)
{
  // initialize grid and ranges
  int cdim = 1, vdim = 2;
  int pdim = cdim+vdim;

  int cells[] = {24, 12, 12};
  int ghost[] = {1, 0, 0};
  double lower[] = {0., -1., -1.};
  double upper[] = {1., 1., 1.};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  // initialize basis
  int poly_order = 2;
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in velocity space). */
    gkyl_cart_modal_hybrid(&basis, cdim, vdim);
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // initialize eqn
  struct gkyl_dg_eqn *eqn;
  enum gkyl_field_id field_id = GKYL_FIELD_E_B;
  enum gkyl_model_id model_id = GKYL_MODEL_DEFAULT;
  eqn = gkyl_dg_vlasov_new(&confBasis, &basis, &confRange, &phaseRange, model_id, field_id, use_gpu);

  // initialize hyper_dg slvr
  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2};
  int zero_flux_flags[2*GKYL_MAX_DIM] = {0, 1, 1, 0, 1, 1};

  gkyl_hyper_dg *slvr;
  slvr = gkyl_hyper_dg_new(&phaseGrid, &basis, eqn, pdim, up_dirs, zero_flux_flags, 1, use_gpu);

  // initialize arrays
  struct gkyl_array *fin, *rhs, *cflrate, *qmem;
  struct gkyl_array *fin_h, *qmem_h, *rhs_h;
  
  fin = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  rhs = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  cflrate = mkarr1(use_gpu, 1, phaseRange_ext.volume);
  qmem = mkarr1(use_gpu, 8*confBasis.num_basis, confRange_ext.volume);

  double *cfl_ptr;
  if (use_gpu)
    cfl_ptr = gkyl_cu_malloc(sizeof(double));
  else
    cfl_ptr = gkyl_malloc(sizeof(double));

  // set initial condition
  int nf = phaseRange_ext.volume*basis.num_basis;
  double *fin_d;
  if (use_gpu) {
    fin_h = mkarr1(false, basis.num_basis, phaseRange_ext.volume);
    fin_d = fin_h->data;
  } else {
    fin_d = fin->data;
  }
  for(int i=0; i< nf; i++) {
    fin_d[i] = (double)(2*i+11 % nf) / nf  * ((i%2 == 0) ? 1 : -1);
  }
  if (use_gpu) gkyl_array_copy(fin, fin_h);

  int nem = confRange_ext.volume*confBasis.num_basis;
  double *qmem_d;
  if (use_gpu) {
    qmem_h = mkarr1(false, 8*confBasis.num_basis, confRange_ext.volume);
    qmem_d = qmem_h->data;
  } else {
    qmem_d = qmem->data;
  }
  for(int i=0; i< nem; i++) {
    qmem_d[i] = (double)(-i+27 % nem) / nem  * ((i%2 == 0) ? 1 : -1);
  }
  if (use_gpu) gkyl_array_copy(qmem, qmem_h);

  // run hyper_dg_advance
  int nrep = 10;
  for(int n=0; n<nrep; n++) {
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_vlasov_set_auxfields(eqn,
      (struct gkyl_dg_vlasov_auxfields) {.field = qmem, .cot_vec = 0, 
      .alpha_surf = 0, .sgn_alpha_surf = 0, .const_sgn_alpha = 0 }); // Must set EM fields to use.

    gkyl_hyper_dg_advance(slvr, &phaseRange, fin, cflrate, rhs);

    gkyl_array_reduce(cfl_ptr, cflrate, GKYL_MAX);
  }

  double cfl_ptr_h[1];
  if (use_gpu)
    gkyl_cu_memcpy(cfl_ptr_h, cfl_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    cfl_ptr_h[0] = cfl_ptr[0];
  TEST_CHECK( gkyl_compare_double(cfl_ptr_h[0], 1.2589437866921e+02, 1e-12) );

  // get linear index of first non-ghost cell
  // 1-indexed for interfacing with G2 Lua layer
  int idx[] = {1, 1, 1, 1, 1};
  int linl = gkyl_range_idx(&phaseRange, idx);

  rhs_h = mkarr1(false, basis.num_basis, phaseRange_ext.volume);
  gkyl_array_copy(rhs_h, rhs);

  // check that ghost cells are empty
  double val = 0;
  double *rhs_d;
  int i = 0;
  while (val==0) {
    rhs_d = gkyl_array_fetch(rhs_h, i);
    val = rhs_d[0];
    if(val==0) i++;
  }
  TEST_CHECK(i == linl);

  // check data in first non-ghost cell
  rhs_d = gkyl_array_fetch(rhs_h, linl);

  //printf("first cell rhs\n");
  //for(int i=0; i<rhs->ncomp; i++) printf("%.16e\n", rhs_d[i]);
  TEST_CHECK( gkyl_compare_double(rhs_d[0],  1.1873679155168162e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[1],  2.2131609273296690e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[2],  -5.1421989058620037e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[3],  -2.1015559414749911e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[4],  -3.1274299616328875e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[5],  -5.7980678014735298e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[6],  -5.1794065426655207e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[7],  1.6555476342162777e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[8],  -3.7839306076374341e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[9],  6.0821819985513692e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[10],  -3.0765279980185458e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[11],  -2.4545628123190717e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[12],  -1.1912775540624752e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[13],  9.1266993445490581e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[14],  -4.5599700127748999e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[15],  1.1971791673432133e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[16],  3.6053906660129060e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[17],  -2.4385757910567278e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[18],  -1.1795732611711490e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[19],  1.8857755500762241e+01, 1e-12) );

  // get linear index of some other cell
  // 1-indexed for interfacing with G2 Lua layer
  int idx2[] = {6, 3, 5};
  int linl2 = gkyl_range_idx(&phaseRange, idx2);
  rhs_d = gkyl_array_fetch(rhs_h, linl2);

  //printf("second cell rhs\n");
  //for(int i=0; i<rhs->ncomp; i++) printf("%.16e\n", rhs_d[i]);
  TEST_CHECK( gkyl_compare_double(rhs_d[0],  7.9777292154304547e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[1],  -2.7182122357080729e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[2],  -3.1823319852967562e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[3],  -1.3560732323031073e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[4],  -9.3361523823747120e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[5],  -6.4546439524194774e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[6],  -3.0046857486230238e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[7],  5.7739758711920011e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[8],  -5.5255780047650804e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[9],  5.4569313596808011e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[10],  -9.3082985569032914e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[11],  -5.9264062443889756e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[12],  -3.4252162191461366e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[13],  5.2936631523369066e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[14],  -5.5255780047650571e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[15],  3.7711009940704315e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[16],  2.8771855264879638e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[17],  -5.8883885687489283e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[18],  -3.4263348605987986e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[19],  4.1758987770172105e+01, 1e-12) );

  // clean up
  gkyl_array_release(fin);
  gkyl_array_release(rhs);
  gkyl_array_release(rhs_h);
  gkyl_array_release(cflrate);
  gkyl_array_release(qmem);

  gkyl_hyper_dg_release(slvr);
  gkyl_dg_eqn_release(eqn);

  if (use_gpu) {
    gkyl_cu_free(cfl_ptr);
    
    gkyl_array_release(fin_h);
    gkyl_array_release(qmem_h);
  }
  else {
    gkyl_free(cfl_ptr);
  }
}

void
test_vlasov_2x3v_p1_(bool use_gpu)
{
  // initialize grid and ranges
  int cdim = 2, vdim = 3;
  int pdim = cdim+vdim;

  int cells[] = {8, 8, 8, 8, 8};
  int ghost[] = {1, 1, 0, 0, 0};
  double lower[] = {0., 0., -1., -1., -1.};
  double upper[] = {1., 1., 1., 1., 1.};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  // initialize basis
  int poly_order = 1;
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in velocity space). */
    gkyl_cart_modal_hybrid(&basis, cdim, vdim);
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // initialize eqn
  struct gkyl_dg_eqn *eqn;
  enum gkyl_field_id field_id = GKYL_FIELD_E_B;
  enum gkyl_model_id model_id = GKYL_MODEL_DEFAULT;
  eqn = gkyl_dg_vlasov_new(&confBasis, &basis, &confRange, &phaseRange, model_id, field_id, use_gpu);

  // initialize hyper_dg slvr
  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2, 3, 4};
  int zero_flux_flags[2*GKYL_MAX_DIM] = {0, 0, 1, 1, 1, 0, 0, 1, 1, 1};

  gkyl_hyper_dg *slvr;
  slvr = gkyl_hyper_dg_new(&phaseGrid, &basis, eqn, pdim, up_dirs, zero_flux_flags, 1, use_gpu);

  // initialize arrays
  struct gkyl_array *fin, *rhs, *cflrate, *qmem;
  struct gkyl_array *fin_h, *qmem_h, *rhs_h;
  
  fin = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  rhs = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  cflrate = mkarr1(use_gpu, 1, phaseRange_ext.volume);
  qmem = mkarr1(use_gpu, 8*confBasis.num_basis, confRange_ext.volume);

  // set initial condition
  int nf = phaseRange_ext.volume*basis.num_basis;
  double *fin_d;
  if (use_gpu) {
    fin_h = mkarr1(false, basis.num_basis, phaseRange_ext.volume);
    fin_d = fin_h->data;
  } else {
    fin_d = fin->data;
  }
  for(int i=0; i< nf; i++) {
    fin_d[i] = (double)(2*i+11 % nf) / nf  * ((i%2 == 0) ? 1 : -1);
  }
  if (use_gpu) gkyl_array_copy(fin, fin_h);

  int nem = confRange_ext.volume*confBasis.num_basis;
  double *qmem_d;
  if (use_gpu) {
    qmem_h = mkarr1(false, 8*confBasis.num_basis, confRange_ext.volume);
    qmem_d = qmem_h->data;
  } else {
    qmem_d = qmem->data;
  }
  for(int i=0; i< nem; i++) {
    qmem_d[i] = (double)(-i+27 % nem) / nem  * ((i%2 == 0) ? 1 : -1);
  }
  if (use_gpu) gkyl_array_copy(qmem, qmem_h);

  // run hyper_dg_advance
  int nrep = 10;
  for(int n=0; n<nrep; n++) {
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_vlasov_set_auxfields(eqn,
      (struct gkyl_dg_vlasov_auxfields) {.field = qmem, .cot_vec = 0, 
      .alpha_surf = 0, .sgn_alpha_surf = 0, .const_sgn_alpha = 0 }); // must set EM fields to use

    gkyl_hyper_dg_advance(slvr, &phaseRange, fin, cflrate, rhs);
  }

  // get linear index of first non-ghost cell
  // 1-indexed for interfacing with G2 Lua layer
  int idx[] = {1, 1, 1, 1, 1};
  int linl = gkyl_range_idx(&phaseRange, idx);

  rhs_h = mkarr1(false, basis.num_basis, phaseRange_ext.volume);
  gkyl_array_copy(rhs_h, rhs);

  // check that ghost cells are empty
  double val = 0;
  double *rhs_d;
  int i = 0;
  while (val==0) {
    rhs_d = gkyl_array_fetch(rhs_h, i);
    val = rhs_d[0];
    if(val==0) i++;
  }
  TEST_CHECK(i == linl);

  // check data in first non-ghost cell
  rhs_d = gkyl_array_fetch(rhs_h, linl);

//  printf("first cell rhs\n");
//  for(int i=0; i<rhs->ncomp; i++) printf("%.16e\n", rhs_d[i]);
  TEST_CHECK( gkyl_compare_double(rhs_d[0 ], 1.5583125504656174e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[1 ], 2.6228222340877672e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[2 ], 2.3673600178359444e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[3 ], 2.0445381301295842e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[4 ], 5.8009947700332258e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[5 ],-4.3951673951369203e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[6 ],-3.2740808507907353e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[7 ], 6.1685466145806656e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[8 ],-9.9225900289967033e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[9 ], 1.5556139476411834e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[10],-2.0649266813526133e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[11],-2.0472420155729774e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[12],-1.6861158492781747e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[13], 5.7160667881100276e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[14],-3.7334543012748869e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[15],-4.0631107293252215e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[16],-2.2368678996807198e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[17], 2.4975896418604588e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[18],-1.7766928843531296e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[19], 6.9681762700869800e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[20],-2.5402211148036020e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[21], 2.5963072988336688e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[22],-1.2065136079090388e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[23], 1.4311020988985428e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[24],-1.7307228807072821e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[25], 1.0648291736093665e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[26],-2.8942763774151025e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[27], 3.4876234221247344e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[28],-1.6213472698033865e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[29],-1.3152058297652811e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[30], 5.5160829801207081e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[31], 1.3248009285063111e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[32], 2.0626248834060192e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[33],-1.8675890171358858e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[34], 1.0700556896733145e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[35], 1.0267345121490179e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[36],-2.9814328597585448e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[37], 8.2810622260402198e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[38],-2.8000677504722145e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[39], 2.5161866916865655e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[40], 1.6582163376892304e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[41],-2.2287681627701883e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[42],-6.5143195699025309e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[43], 7.6733697089938193e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[44], 1.6105940484870842e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[45], 2.5284590774855953e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[46],-1.6803923116017092e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[47], 3.4993021915577259e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[48],-6.4429006203071723e-02, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[49], 2.1373817491282676e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[50],-9.0822803896705562e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[51], 6.3081769272408259e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[52],-1.1470246109901016e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[53], 2.7774567065033185e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[54],-2.7118393077830891e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[55], 1.9747275150178819e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[56],-1.8774237677748375e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[57],-4.0091254392929780e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[58], 1.1181326827391735e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[59], 1.0012387003982573e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[60],-1.5378762139131592e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[61], 1.9579525074768895e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[62],-7.2504813611281804e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[63], 2.8594267179149174e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[64], 3.3960460641818287e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[65],-1.2671905220786885e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[66], 2.3465937131321578e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[67],-2.5018574073021767e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[68],-7.6358915438484063e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[69],-4.8114279712708772e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[70], 4.2101152078896291e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[71],-1.0127724148100189e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[72],-5.7706791067934295e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[73],-1.0415188306383120e-01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[74],-5.6167239912421829e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[75], 3.9890443806074821e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[76],-2.0747338814618409e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[77], 2.6252175357150804e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[78],-1.5452629233319728e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[79], 3.5165336694630724e+01, 1e-12) );

  // get linear index of some other cell
  // 1-indexed for interfacing with G2 Lua layer
  int idx2[] = {6, 3, 5, 8, 2};
  int linl2 = gkyl_range_idx(&phaseRange, idx2);
  rhs_d = gkyl_array_fetch(rhs_h, linl2);

//  printf("second cell rhs\n");
//  for(int i=0; i<rhs->ncomp; i++) printf("%.16e\n", rhs_d[i]);
  TEST_CHECK( gkyl_compare_double(rhs_d[0 ], 4.7846096908267584e-02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[1 ], 1.1169581543468636e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[2 ],-5.7985569827277743e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[3 ], 6.3560588230111681e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[4 ],-3.8037040250162341e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[5 ], 3.1215390309173313e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[6 ],-5.9451450645977602e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[7 ], 1.4624629486937456e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[8 ],-4.9690031412084231e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[9 ], 3.5213109571593453e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[10],-6.0955021313238532e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[11], 2.6476452020050101e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[12],-3.5100041616767794e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[13], 4.8316996719345781e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[14],-6.1122488642299089e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[15], 7.9850069590596184e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[16],-6.4345760693277171e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[17], 5.5391125581775690e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[18],-7.3134438613368351e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[19], 5.2492380773845227e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[20],-6.0283166541522384e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[21], 7.9945201071625741e-02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[22],-5.7177859264242791e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[23], 1.0385552375417682e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[24],-4.5437126504480567e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[25], 7.7523547979949659e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[26],-5.2901252593887023e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[27], 5.6432376028981565e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[28],-6.3928564161796082e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[29], 1.4801034042997721e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[30],-6.1595791904709372e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[31], 6.8347478052667171e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[32], 7.7635047936775073e-02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[33], 1.1223967073645461e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[34],-4.8555037059507058e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[35], 1.2892406386095562e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[36],-3.9986417153680837e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[37], 6.1091198999195505e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[38],-1.1056795014088960e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[39], 4.8890197502476475e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[40],-3.2270859286785982e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[41], 4.8833583888325073e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[42], 5.4541037169184979e-02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[43], 6.0485697373310906e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[44],-5.1988442558875640e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[45], 1.1184035513346096e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[46],-4.8157508893945355e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[47], 5.9610415538418607e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[48],-3.4882339711782762e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[49], 3.5298202468643511e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[50],-5.6203287350628727e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[51], 3.4882339711782717e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[52],-3.6695369052217031e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[53], 6.7906001630206347e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[54],-3.0448691264619017e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[55], 5.7296223828675863e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[56],-2.8163526802784560e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[57], 4.9266605122525164e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[58],-5.5910759361390594e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[59], 6.0312249472335431e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[60],-5.9820310601514130e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[61], 3.9298520202026186e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[62],-5.6639340789701137e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[63], 5.9784613439394072e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[64],-3.1215390309173785e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[65], 3.4700765537768121e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[66],-4.7918000752070810e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[67], 2.8905989232414975e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[68],-3.9094010767585058e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[69], 5.9660773305385845e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[70],-2.9851254333743560e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[71], 4.9050915156881537e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[72],-2.7772648983706874e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[73], 4.8843130373068703e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[74],-5.2000000000000235e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[75], 6.8558300614068344e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[76],-5.1539246260931648e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[77], 3.9101083271150761e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[78],-4.8358250228732551e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[79], 6.0284545971110632e+01, 1e-12) );

  // clean up
  gkyl_array_release(fin);
  gkyl_array_release(rhs);
  gkyl_array_release(rhs_h);
  gkyl_array_release(cflrate);
  gkyl_array_release(qmem);

  gkyl_hyper_dg_release(slvr);
  gkyl_dg_eqn_release(eqn);

  if (use_gpu) {
    gkyl_array_release(fin_h);
    gkyl_array_release(qmem_h);
  }
}


void
test_vlasov_1x2v_p2()
{
  test_vlasov_1x2v_p2_(false);
}

void
test_vlasov_1x2v_p2_cu()
{
  test_vlasov_1x2v_p2_(true);
}

void
test_vlasov_2x3v_p1()
{
  test_vlasov_2x3v_p1_(false);
}

void
test_vlasov_2x3v_p1_cu()
{
  test_vlasov_2x3v_p1_(true);
}

#ifndef GKYL_HAVE_CUDA
int hyper_dg_kernel_test(const gkyl_hyper_dg *slvr) {
  return 0;
}
#endif

TEST_LIST = {
  { "test_vlasov_1x2v_p2", test_vlasov_1x2v_p2 },
  { "test_vlasov_2x3v_p1", test_vlasov_2x3v_p1 },
#ifdef GKYL_HAVE_CUDA
  { "test_vlasov_1x2v_p2_cu", test_vlasov_1x2v_p2_cu },
  { "test_vlasov_2x3v_p1_cu", test_vlasov_2x3v_p1_cu },
#endif
  { NULL, NULL },
};

