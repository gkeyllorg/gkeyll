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
test_vlasov_3x3v_p1_(bool use_gpu)
{
  // initialize grid and ranges
  int cdim = 3, vdim = 3;
  int pdim = cdim+vdim;

  int cells[] = {8, 8, 8, 8, 8, 8};
  int ghost[] = {1, 1, 1, 0, 0, 0};
  double lower[] = {0., 0., 0., -1., -1., -1.};
  double upper[] = {1., 1., 1., 1., 1., 1.};

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
  enum gkyl_field_id field_id = GKYL_FIELD_NULL;
  enum gkyl_model_id model_id = GKYL_MODEL_DEFAULT;

  // initialize arrays
  struct gkyl_array *fin, *rhs, *cflrate;
  struct gkyl_array *fin_h, *rhs_h;
  
  fin = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  rhs = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  cflrate = mkarr1(use_gpu, 1, phaseRange_ext.volume);
  struct gkyl_array *qmem = mkarr1(use_gpu, 8*confBasis.num_basis, confRange_ext.volume);

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

  // Select the number of nodes, with case for hybrid-tensor.
  bool use_lo = false;
  int highorder = use_lo ? 0 : 1;
  int num_surf_vel_nodes = vdim*pow(poly_order+1+highorder,pdim - 1);
  if ((basis.b_type == GKYL_BASIS_MODAL_TENSOR) && (poly_order == 1)) {
    num_surf_vel_nodes = (int) vdim*(pow(poly_order+1+highorder,vdim - 1) + pow(poly_order,cdim));
  }

  // Sturcture pointers for input objects (but not used)
  int num_pt_indices[3] = { 1 , 6, 18 }; 
  struct gkyl_array *poisson_tensor_conf = mkarr1(use_gpu, confBasis.num_basis*num_pt_indices[vdim-1], confRange.volume );
  struct gkyl_array *pot_tot = mkarr1(use_gpu, confBasis.num_basis*4, confRange_ext.volume );
  struct gkyl_array *vel_flux_surf = mkarr1(use_gpu, num_surf_vel_nodes, phaseRange_ext.volume );
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
    .has_E = false, 
    .has_phi = false, 
    .has_B = false, 
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
    .has_E = false, 
    .has_phi = false, 
    .has_B = false, 
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
  // FIELD_NULL so only configuration space update, no velocity space update
  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2};
  int zero_flux_flags[2*GKYL_MAX_DIM] = {0, 0, 0, 0, 0, 0};
  int num_up_dirs = cdim; 

  gkyl_hyper_dg *slvr;
  slvr = gkyl_hyper_dg_new(&phaseGrid, &basis, eqn, num_up_dirs, up_dirs, zero_flux_flags, 1, use_gpu);

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
  int idx[] = {1, 1, 1, 1, 1, 1};
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

  static const double expected_cell1[] = {
  3.7606327288427490e+00,  1.6013126194947759e+01, -5.6384122711121769e+00,  5.6678864833199523e+00,
 -6.5263272884274992e-01,  6.5263272884274881e-01,  3.7606327288427495e+00,  2.5412199393144888e+01,
 -2.5441676168829229e+01,  1.5066964808470219e+01, -1.0726930201536613e+01,  1.1021649252325320e+01,
 -9.9841801670706705e+00,  1.0726937891966298e+01, -1.1021656942755007e+01,  9.9841878575003484e+00,
  3.7606327288427490e+00,  1.5576690407065403e+01, -1.0536692970541964e+01,  1.0032695534018536e+01,
  3.7606327288427477e+00, -3.7606327288427477e+00, -3.4840844215659160e+01,  2.1047513389321317e+01,
 -2.0543518516274446e+01,  1.5503523643227570e+01, -2.1047528770180691e+01,  2.0543533897133813e+01,
 -1.5503539024086948e+01,  1.5576721168784152e+01, -1.0536723732260723e+01,  1.0032726295737287e+01,
 -2.1047559531899442e+01,  2.0543564658852567e+01, -1.5503569785805698e+01,  1.5576736549643529e+01,
 -1.0536739113120099e+01,  1.0032741676596663e+01, -1.5576744240073216e+01,  1.0536746803549788e+01,
 -1.0032749367026351e+01, -3.7606327288427486e+00, -3.4840998024252904e+01,  3.4841005714682588e+01,
 -2.1047621055336929e+01,  2.0543626182290069e+01, -1.5503631309243186e+01,  3.4841036476401342e+01,
 -2.1096138985761371e+01,  2.0058672463983271e+01, -2.0353394078248542e+01,  2.1096154366620748e+01,
 -2.0058687844842659e+01,  2.0353409459107930e+01, -1.6013262059205573e+01,  5.6385481353700015e+00,
 -5.6680223475777760e+00,  3.4841113380698211e+01, -3.0427855613442404e+01,  3.0427863303872080e+01,
 -2.5412471121660499e+01,  2.5441947897344839e+01, -1.5067236536985821e+01,  3.4841159523276339e+01};

  const double expected_cell2[] = {
  9.8991476114550292e-01,  5.1674656662644999e+01, -3.7217126310250841e+01,  5.2529416951017183e+01,
 -3.5619147611455007e+00,  3.5619147611455038e+00,  9.8991476114550059e-01,  8.9881722270174166e+01,
 -1.0519401913652642e+02,  9.0736486220655692e+01, -4.6527024733839347e+01,  4.2461792635061585e+01,
 -5.6984274699373394e+01,  4.6527032424269045e+01, -4.2461798128225645e+01,  5.6984282389803070e+01,
  9.8991476114550359e-01,  5.1376784939368143e+01, -4.2115395290930621e+01,  5.7032790066321240e+01,
  9.8991476114550414e-01, -9.8991476114550359e-01, -1.4340116860636027e+02,  8.5378460482995067e+01,
 -1.0029586148397163e+02,  9.1034469272057564e+01, -8.5378473666588832e+01,  1.0029587686483099e+02,
 -9.1034482455651300e+01,  5.1376815701086905e+01, -4.2115417263586870e+01,  5.7032820828040009e+01,
 -8.5378500033776319e+01,  1.0029590762654976e+02, -9.1034508822838802e+01,  5.1376831081946278e+01,
 -4.2115428249914999e+01,  5.7032836208899369e+01, -5.1376838772375969e+01,  4.2115433743079059e+01,
 -5.7032843899329066e+01, -9.8991476114550037e-01, -1.4340130776651654e+02,  1.4340131472452435e+02,
 -8.5378552768151323e+01,  1.0029596914998726e+02, -9.1034561557213806e+01,  1.4340134255655559e+02,
 -8.5427067768888236e+01,  9.9949579496285992e+01, -9.5884321396531647e+01,  8.5427080952481987e+01,
 -9.9949594877145358e+01,  9.5884334580125412e+01, -5.1674792526902813e+01,  3.7217223356149283e+01,
 -5.2529552815274982e+01,  1.4340141213663372e+02, -1.3884958957235051e+02,  1.3884959653035835e+02,
 -8.9881955180330422e+01,  1.0519429086504206e+02, -9.0736719130811949e+01,  1.4340145388468056e+02};


  //printf("\nfirst cell rhs\n");
  //for(int i=0; i<rhs->ncomp; i++) printf("%.16e\n", rhs_d[i]);
  //for(int i=0; i<rhs->ncomp; i++) printf("  TEST_CHECK( gkyl_compare_double(rhs_d[%d], %.16e, 1e-12) );\n", i, rhs_d[i]);
  for (int i = 0; i < rhs->ncomp; ++i) {
    TEST_CHECK( gkyl_compare_double(rhs_d[i], expected_cell1[i], 1e-12) );
  }

/*   // get linear index of some other cell */
/*   // 1-indexed for interfacing with G2 Lua layer */
  int idx2[] = {6, 3, 5, 8, 2, 1};
  int linl2 = gkyl_range_idx(&phaseRange, idx2);
  rhs_d = gkyl_array_fetch(rhs_h, linl2);

  //printf("second cell rhs\n");
  //for(int i=0; i<rhs->ncomp; i++) printf("%.16e\n", rhs_d[i]);
  //for(int i=0; i<rhs->ncomp; i++) printf("  TEST_CHECK( gkyl_compare_double(rhs_d[%d], %.16e, 1e-12) );\n", i, rhs_d[i]);
  for (int i = 0; i < rhs->ncomp; ++i) {
    TEST_CHECK( gkyl_compare_double(rhs_d[i], expected_cell2[i], 1e-12) );
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
  gkyl_array_release(qmem);

  gkyl_array_release(fin);
  gkyl_array_release(rhs);
  gkyl_array_release(rhs_h);
  gkyl_array_release(cflrate);

  gkyl_hyper_dg_release(slvr);
  gkyl_dg_eqn_release(eqn);
  gkyl_vlasov_velocity_map_release(vel_map);
  gkyl_vlasov_position_map_release(pos_map);

  if (use_gpu) {
    gkyl_array_release(fin_h);
  }
}

void
test_vlasov_3x3v_p1()
{
  test_vlasov_3x3v_p1_(false);
}

TEST_LIST = {
  { "test_vlasov_3x3v_p1", test_vlasov_3x3v_p1 },
  { NULL, NULL },
};

