#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_basis.h>
#include <gkyl_dg_vlasov_calc_hamil.h>
#include <gkyl_dg_vlasov_vel_flux_surf.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_vlasov_velocity_map.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>

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

  int velcells[] = {12, 12};
  int velghost[] = {0, 0};
  double vellower[] = {-1., -1.};
  double velupper[] = {1., 1.};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid velGrid;
  struct gkyl_range velRange, velRange_ext;
  gkyl_rect_grid_init(&velGrid, vdim, vellower, velupper, velcells);
  gkyl_create_grid_ranges(&velGrid, velghost, &velRange_ext, &velRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  // initialize basis
  int poly_order = 2;
  struct gkyl_basis basis, confBasis, velBasis; // phase-space, conf-space basis

  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  enum gkyl_field_id field_id = GKYL_FIELD_E_B;
  enum gkyl_model_id model_id = GKYL_MODEL_DEFAULT;

  // initialize arrays
  struct gkyl_array *fin, *rhs, *cflrate, *qmem;
  struct gkyl_array *fin_h, *qmem_h, *rhs_h;
  fin = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  rhs = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  cflrate = mkarr1(use_gpu, 1, phaseRange_ext.volume);
  qmem = mkarr1(use_gpu, 8*confBasis.num_basis, confRange_ext.volume);

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

  // build hamil and gamma_inv
  struct gkyl_array *hamil = mkarr1(use_gpu, velBasis.num_basis, velRange.volume);
  struct gkyl_array *gamma_inv = mkarr1(use_gpu, velBasis.num_basis, velRange.volume);
  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = { 0 };
  struct gkyl_vlasov_velocity_map *vel_map = gkyl_vlasov_velocity_map_new(&velGrid,
    &velRange, &velBasis, inp_vmap, false, use_gpu);
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = { 0 };
  struct gkyl_vlasov_position_map *pos_map = gkyl_vlasov_position_map_new(&confGrid,
    &confRange, &confRange_ext, &confBasis, inp_pmap, use_gpu);

  gkyl_dg_vlasov_calc_hamil(&velGrid, &velBasis, &velRange, 
    GKYL_MODEL_DEFAULT, vel_map, hamil, gamma_inv, use_gpu);

  // Sturcture pointers for input objects (but not used)
  int num_pt_indices[3] = { 1 , 6, 18 }; 
  struct gkyl_array *poisson_tensor_conf = mkarr1(use_gpu, confBasis.num_basis*num_pt_indices[vdim-1], confRange.volume );
  struct gkyl_array *pot_tot = mkarr1(use_gpu, confBasis.num_basis*4, confRange_ext.volume );
  struct gkyl_array *vel_flux_surf = mkarr1(use_gpu, basis.num_basis*vdim, phaseRange_ext.volume );
  struct gkyl_array *f_no_J = mkarr1(use_gpu, fin->ncomp, fin->size); ;
  struct gkyl_array *rad = mkarr1(use_gpu, vdim*velBasis.num_basis, velRange.volume);

  struct gkyl_dg_vlasov_vel_flux_surf_inp inp_vel_flux = {
    .phase_grid = &phaseGrid, 
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .vel_map = vel_map,
    .pos_map = pos_map,
    .hamil_range = &velRange,
    .skip_cell_thresh = 0.0, 
    .model_id = model_id,
    .hamil_id = gkyl_hamil_id_from_model_id(model_id),
    .has_E = true, 
    .has_phi = false, 
    .has_B = true, 
    .has_rad = false, 
    .use_lo = false,
    .use_gpu = use_gpu,
  }; 
  struct gkyl_dg_vlasov_vel_flux_surf *calc_vel_flux = gkyl_dg_vlasov_vel_flux_surf_inew(&inp_vel_flux); 

  struct gkyl_dg_vlasov_inp inp_eqn = {
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confRange,
    .hamil_range = &velRange,
    .phase_range = &phaseRange,
    .vel_map = vel_map,
    .pos_map = pos_map,
    .skip_cell_thresh = 0.0, 
    .model_id = model_id,
    .hamil_id = gkyl_hamil_id_from_model_id(model_id),
    .has_E = true, 
    .has_phi = false, 
    .has_B = true, 
    .has_rad = false, 
    .poisson_tensor_conf = poisson_tensor_conf,
    .hamil = hamil,
    .qmem = qmem, 
    .pot_tot = pot_tot, 
    .vel_flux_surf = vel_flux_surf, 
    .f_no_J = f_no_J, 
    .rad = rad, 
    .use_lo = false,
    .use_gpu = use_gpu,
  };  
  // Construct Vlasov equation and Hyper DG object for updating equation. 
  struct gkyl_dg_eqn *eqn = gkyl_dg_vlasov_inew(&inp_eqn); 

  // initialize hyper_dg slvr
  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2};
  int zero_flux_flags[2*GKYL_MAX_DIM] = {0, 1, 1, 0, 1, 1};

  gkyl_hyper_dg *slvr;
  slvr = gkyl_hyper_dg_new(&phaseGrid, &basis, eqn, pdim, up_dirs, zero_flux_flags, 1, use_gpu);

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
  gkyl_array_set(f_no_J, 1.0, fin);

  // run hyper_dg_advance
  int nrep = 10;
  for(int n=0; n<nrep; n++) {
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);

    gkyl_dg_vlasov_vel_flux_surf_advance(calc_vel_flux, &confRange, &phaseRange,
     poisson_tensor_conf, hamil, qmem, pot_tot, rad, 
     f_no_J, cflrate, vel_flux_surf);    
    gkyl_hyper_dg_advance(slvr, &phaseRange, fin, cflrate, rhs);

    gkyl_array_reduce(cfl_ptr, cflrate, GKYL_MAX);
  }

  double cfl_ptr_h[1];
  if (use_gpu)
    gkyl_cu_memcpy(cfl_ptr_h, cfl_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    cfl_ptr_h[0] = cfl_ptr[0];
  TEST_CHECK( gkyl_compare_double(cfl_ptr_h[0], 2.0994949632236742e+02, 1e-12) );

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
  const double exp0[] = {
    1.2236372998891958e+00, 2.0105903561583763e-01, -5.1419282362326753e+00, -2.0352048757787333e+00,
    -3.1394800811609201e+01, -5.8085102299277480e+00, -5.1728224653103201e+00, 1.6650094262740243e+01,
    -5.0331198175173064e-01, 6.1804101443921522e-01, -3.0931309818069451e+01, -2.4166487365871351e+01,
    -1.1949509605874356e+01, 9.2106442019431967e+00, -6.6446814763038498e-01, 1.2077341612159545e+01,
    3.4825387939850700e+00, -2.3867523859361569e+01, -1.1804252614108304e+01, 1.8849579738214938e+01
  };
  for (int j = 0; j < rhs->ncomp; ++j) {
    TEST_CHECK(gkyl_compare_double(rhs_d[j], exp0[j], 1e-12));
    TEST_MSG("rhs_d[%d]=%.16e, expected=%.16e", j, rhs_d[j], exp0[j]);
  }

  // get linear index of some other cell
  // 1-indexed for interfacing with G2 Lua layer
  int idx2[] = {6, 3, 5};
  int linl2 = gkyl_range_idx(&phaseRange, idx2);
  rhs_d = gkyl_array_fetch(rhs_h, linl2);

  //printf("second cell rhs\n");
  //for(int i=0; i<rhs->ncomp; i++) printf("%.16e\n", rhs_d[i]);
  const double exp1[] = {
    7.9777292154304547e-01, -2.7182122357080729e+00, -3.1823319852967562e+00, -1.3560732323031073e+00,
    -9.3361523823747120e+01, -6.4546439524194774e+00, -3.0046857486230238e+00, 5.7739758711920011e+01,
    -5.5255780047650804e-01, 5.4569313596808011e-01, -9.3082985569032914e+01, -5.9264062443889756e+01,
    -3.4252162191461366e+01, 5.2936631523369066e+01, -5.5255780047650571e-01, 3.7711009940704315e+01,
    2.8771855264879638e+00, -5.8883885687489283e+01, -3.4263348605987986e+01, 4.1758987770172105e+01
  };
  for (int j = 0; j < rhs->ncomp; ++j) {
    TEST_CHECK(gkyl_compare_double(rhs_d[j], exp1[j], 1e-12));
    TEST_MSG("rhs_d[%d]=%.16e, expected=%.16e", j, rhs_d[j], exp1[j]);
  }

  // clean up
  gkyl_dg_vlasov_vel_flux_surf_release(calc_vel_flux);

  gkyl_array_release(hamil);
  gkyl_array_release(gamma_inv);
  gkyl_array_release(poisson_tensor_conf);
  gkyl_array_release(pot_tot);
  gkyl_array_release(vel_flux_surf);
  gkyl_array_release(f_no_J);
  gkyl_array_release(rad);

  gkyl_array_release(fin);
  gkyl_array_release(rhs);
  gkyl_array_release(rhs_h);
  gkyl_array_release(cflrate);
  gkyl_array_release(qmem);

  gkyl_hyper_dg_release(slvr);
  gkyl_dg_eqn_release(eqn);
  gkyl_vlasov_velocity_map_release(vel_map);
  gkyl_vlasov_position_map_release(pos_map);

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

  int velcells[] = { 8, 8, 8};
  int velghost[] = { 0, 0, 0};
  double vellower[] = { -1., -1., -1.};
  double velupper[] = { 1., 1., 1.};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid velGrid;
  struct gkyl_range velRange, velRange_ext;
  gkyl_rect_grid_init(&velGrid, vdim, vellower, velupper, velcells);
  gkyl_create_grid_ranges(&velGrid, velghost, &velRange_ext, &velRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  // initialize basis
  int poly_order = 1;
  struct gkyl_basis basis, confBasis, velBasis; // phase-space, conf-space basis
  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // initialize eqn
  enum gkyl_field_id field_id = GKYL_FIELD_E_B;
  enum gkyl_model_id model_id = GKYL_MODEL_DEFAULT;

  // initialize arrays
  struct gkyl_array *fin, *rhs, *cflrate, *qmem;
  struct gkyl_array *fin_h, *qmem_h, *rhs_h;
  
  fin = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  rhs = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  cflrate = mkarr1(use_gpu, 1, phaseRange_ext.volume);
  qmem = mkarr1(use_gpu, 8*confBasis.num_basis, confRange_ext.volume);

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

  // build hamil and gamma_inv
  struct gkyl_array *hamil = mkarr1(use_gpu, velBasis.num_basis, velRange.volume);
  struct gkyl_array *gamma_inv = mkarr1(use_gpu, velBasis.num_basis, velRange.volume);
  struct gkyl_vlasov_velocity_map_inp inp_vmap[GKYL_MAX_CDIM] = { 0 };
  struct gkyl_vlasov_velocity_map *vel_map = gkyl_vlasov_velocity_map_new(&velGrid,
    &velRange, &velBasis, inp_vmap, false, use_gpu);
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = { 0 };
  struct gkyl_vlasov_position_map *pos_map = gkyl_vlasov_position_map_new(&confGrid,
    &confRange, &confRange_ext, &confBasis, inp_pmap, use_gpu);

  gkyl_dg_vlasov_calc_hamil(&velGrid, &velBasis, &velRange, 
    GKYL_MODEL_DEFAULT, vel_map, hamil, gamma_inv, use_gpu);

  // Sturcture pointers for input objects (but not used)
  int num_pt_indices[3] = { 1 , 6, 18 }; 
  struct gkyl_array *poisson_tensor_conf = mkarr1(use_gpu, confBasis.num_basis*num_pt_indices[vdim-1], confRange.volume );
  struct gkyl_array *pot_tot = mkarr1(use_gpu, confBasis.num_basis*4, confRange_ext.volume );
  struct gkyl_array *vel_flux_surf = mkarr1(use_gpu, basis.num_basis*vdim, phaseRange_ext.volume );
  struct gkyl_array *f_no_J = mkarr1(use_gpu, fin->ncomp, fin->size); ;
  struct gkyl_array *rad = mkarr1(use_gpu, vdim*velBasis.num_basis, velRange.volume);

  struct gkyl_dg_vlasov_vel_flux_surf_inp inp_vel_flux = {
    .phase_grid = &phaseGrid, 
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .vel_map = vel_map,
    .pos_map = pos_map,
    .hamil_range = &velRange,
    .skip_cell_thresh = 0.0, 
    .model_id = model_id,
    .hamil_id = gkyl_hamil_id_from_model_id(model_id),
    .has_E = true, 
    .has_phi = false, 
    .has_B = true, 
    .has_rad = false, 
    .use_lo = false,
    .use_gpu = use_gpu,
  }; 
  struct gkyl_dg_vlasov_vel_flux_surf *calc_vel_flux = gkyl_dg_vlasov_vel_flux_surf_inew(&inp_vel_flux); 

  struct gkyl_dg_vlasov_inp inp_eqn = {
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confRange,
    .hamil_range = &velRange,
    .phase_range = &phaseRange,
    .vel_map = vel_map,
    .pos_map = pos_map,
    .skip_cell_thresh = 0.0, 
    .model_id = model_id,
    .hamil_id = gkyl_hamil_id_from_model_id(model_id),
    .has_E = true, 
    .has_phi = false, 
    .has_B = true, 
    .has_rad = false, 
    .poisson_tensor_conf = poisson_tensor_conf,
    .hamil = hamil,
    .qmem = qmem, 
    .pot_tot = pot_tot, 
    .vel_flux_surf = vel_flux_surf, 
    .f_no_J = f_no_J, 
    .rad = rad, 
    .use_lo = false,
    .use_gpu = use_gpu,
  };  
  // Construct Vlasov equation and Hyper DG object for updating equation. 
  struct gkyl_dg_eqn *eqn = gkyl_dg_vlasov_inew(&inp_eqn); 

  //eqn = gkyl_dg_vlasov_new(&confBasis, &basis, &confRange, &phaseRange, model_id, field_id, use_gpu);

  // initialize hyper_dg slvr
  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2, 3, 4};
  int zero_flux_flags[2*GKYL_MAX_DIM] = {0, 0, 1, 1, 1, 0, 0, 1, 1, 1};

  gkyl_hyper_dg *slvr;
  slvr = gkyl_hyper_dg_new(&phaseGrid, &basis, eqn, pdim, up_dirs, zero_flux_flags, 1, use_gpu);

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
   gkyl_array_set(f_no_J, 1.0, fin);

  // run hyper_dg_advance
  int nrep = 10;
  for(int n=0; n<nrep; n++) {
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_dg_vlasov_vel_flux_surf_advance(calc_vel_flux, &confRange, &phaseRange,
     poisson_tensor_conf, hamil, qmem, pot_tot, rad, 
     f_no_J, cflrate, vel_flux_surf);   
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
  //for(int i=0; i<rhs->ncomp; i++) printf("%.16e\n", rhs_d[i]);
  const double exp0[] = {
    5.8503928145486075e+00, 1.3348107973742653e+01, -8.3015588600670664e+00, -5.1601563129102523e-01,
    3.3662725617516118e+00, -6.0220290457402514e+00, -2.3066605359371479e+01, 1.1795492527988509e+01,
    -1.2094070618291957e+01, 1.6328023879915964e+01, -5.9571695141905820e+00, -9.4566262621926439e-02,
    -1.3664800178573332e+01, 3.2863247869002210e+00, -2.3221330832899003e+00, -5.5061118408615184e-01,
    -1.9749990082666812e+01, 2.5652234998826916e+01, -1.9592570200054443e+01, 9.2140945629575306e+00,
    -2.2973923054182546e+01, 2.2017046005263165e+01, -1.1646193018869390e+01, 1.1830908772108357e+01,
    -1.2129487374637437e+01, -5.9473538681355222e+00, -2.8902000068602643e+01, 3.1341669905686420e+01,
    -1.9716623728746605e+01, 1.7728926461597514e+01, -1.2708781703202801e+01, 3.1935276140332110e+01};

  for (int i = 0; i < rhs->ncomp; ++i) {
    TEST_CHECK(gkyl_compare_double(rhs_d[i], exp0[i], 1e-12));
  }

  // get linear index of some other cell
  // 1-indexed for interfacing with G2 Lua layer
  int idx2[] = {6, 3, 5, 8, 2};
  int linl2 = gkyl_range_idx(&phaseRange, idx2);
  rhs_d = gkyl_array_fetch(rhs_h, linl2);

  //  printf("second cell rhs\n");
  //for(int i=0; i<rhs->ncomp; i++) printf("%.16e\n", rhs_d[i]);
  const double exp1[] = {
    -2.3607695154586722e-01, 7.2773300308759890e+00, -5.3317823194938491e+01, 4.4392304845413322e-01,
    -2.3607695154587002e-01, 2.3607695154586483e-01, -5.9943625779670839e+01, 6.5845536531609365e+00,
    -5.3803105038245278e+01, 7.7623628507452853e+00, -5.2625413028160921e+01, 2.3607695154587174e-01,
    -7.7623848234015105e+00, 5.2625566836754672e+01, -2.3607695154586925e-01, 4.4392304845412967e-01,
    -5.9944211717170838e+01, 6.0152116407829062e+01, -7.7624287687140043e+00, 5.2625874453942174e+01,
    -6.0152292189079098e+01, 7.7624507413702837e+00, -5.2626028262535925e+01, 6.5846708406609364e+00,
    -5.3803925350745281e+01, 4.4392304845413200e-01, -6.0152643751579092e+01, 6.0152702345329089e+01,
    -5.9944914842170824e+01, 7.0696890120927236e+00, -5.3111412644905229e+01, 6.0152936720329066e+01};

  for (int i = 0; i < rhs->ncomp; ++i) {
    TEST_CHECK(gkyl_compare_double(rhs_d[i], exp1[i], 1e-12));
  }

  // clean up
  gkyl_dg_vlasov_vel_flux_surf_release(calc_vel_flux);
  
  gkyl_array_release(hamil);
  gkyl_array_release(gamma_inv);
  gkyl_array_release(poisson_tensor_conf);
  gkyl_array_release(pot_tot);
  gkyl_array_release(vel_flux_surf);
  gkyl_array_release(f_no_J);
  gkyl_array_release(rad);

  gkyl_array_release(fin);
  gkyl_array_release(rhs);
  gkyl_array_release(rhs_h);
  gkyl_array_release(cflrate);
  gkyl_array_release(qmem);

  gkyl_hyper_dg_release(slvr);
  gkyl_dg_eqn_release(eqn);
  gkyl_vlasov_velocity_map_release(vel_map);
  gkyl_vlasov_position_map_release(pos_map);

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

