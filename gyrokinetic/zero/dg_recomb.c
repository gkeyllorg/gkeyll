#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_recomb.h>
#include <gkyl_dg_recomb_priv.h>
#include <gkyl_util.h>
#include <gkyl_const.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_rect_decomp.h>

struct gkyl_dg_recomb*
gkyl_dg_recomb_new(struct gkyl_dg_recomb_inp *inp, bool use_gpu)
{
  gkyl_dg_recomb *up = gkyl_malloc(sizeof(struct gkyl_dg_recomb));

  up->cbasis = inp->cbasis;
  up->pbasis = inp->pbasis;
  up->conf_rng = inp->conf_rng;
  up->conf_rng_ext = inp->conf_rng_ext;
  up->phase_rng = inp->phase_rng;
  up->grid = inp->grid;
  
  int cdim = up->cbasis->ndim;
  int pdim = up->pbasis->ndim;
  int vdim = pdim - cdim;
  int poly_order = up->cbasis->poly_order;
  up->cdim = cdim;
  up->vdim = vdim;
  up->use_gpu = use_gpu;
  
  up->elem_charge = GKYL_ELEMENTARY_CHARGE;
  up->mass_elc = GKYL_ELECTRON_MASS;

  int charge_state = inp->charge_state;
  enum gkyl_ion_type type_ion = inp->type_ion;
  
  // Project ADAS data
  struct adas_field data;

  read_adas_field_recomb(type_ion, &data);
  
  long sz = data.NT*data.NN;
  double minmax[2];

  if (data.logT == NULL) fprintf(stderr, "Unable to load ADAS 'logT_<elem>.npy' file.");
  if (data.logN == NULL) fprintf(stderr, "Unable to load ADAS 'logN_<elem>.npy' file.");
  if (data.logData == NULL) fprintf(stderr, "Unable to load ADAS 'recomb_<elem>.npy' file.");

  minmax_from_numpy(data.logT, data.NT, minmax);
  fclose(data.logT);
  double logTmin = minmax[0], logTmax = minmax[1];
  minmax_from_numpy(data.logN, data.NN, minmax);
  fclose(data.logN);
  double logNmin = minmax[0]+6., logNmax = minmax[1]+6.; //adjust for 1/cm^3 to 1/m^3 conversion

  struct gkyl_array *adas_nodal = gkyl_array_new(GKYL_DOUBLE, 1, sz);
  array_from_numpy(data.logData, sz, data.Zmax, charge_state, adas_nodal);
  fclose(data.logData);

  if (!adas_nodal) {
    fprintf(stderr, "Unable to read data from adas nodal numpy file!\n");
    return 0;
  }

  struct gkyl_range range_nodal;
  gkyl_range_init_from_shape(&range_nodal, 2, (int[]) { data.NT, data.NN } );
  
  // allocate grid and DG array
  struct gkyl_rect_grid tn_grid;
  gkyl_rect_grid_init(&tn_grid, 2,
    (double[]) { logTmin, logNmin},
    (double []) { logTmax, logNmax},
    (int[]) { data.NT-1, data.NN-1 }
  );

  if (use_gpu) {
    // allocate device basis if we are using GPUs
    up->basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  }
  else {
    up->basis_on_dev = &up->adas_basis;
  }
  gkyl_cart_modal_serendip(&up->adas_basis, 2, 1);
  if (use_gpu)
    gkyl_cart_modal_serendip_cu_dev(up->basis_on_dev, 2, 1);

  int ghost[GKYL_MAX_DIM] = { 1, 1};
  struct gkyl_range modal_range;
  struct gkyl_range modal_range_ext;
  gkyl_create_grid_ranges(&tn_grid, ghost, &modal_range_ext, &modal_range);

  struct gkyl_array *adas_dg = gkyl_array_new(GKYL_DOUBLE, up->adas_basis.num_basis, modal_range_ext.volume);

  //create_dg_from_nodal(&tn_grid, &range_node, adas_nodal, adas_dg, charge_state);

  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&up->adas_basis, &tn_grid, false);
  gkyl_nodal_ops_n2m(n2m, &up->adas_basis, &tn_grid, &range_nodal, &modal_range, 1, adas_nodal, adas_dg, false);
  gkyl_nodal_ops_release(n2m);

  // ADAS data pointers
  up->recomb_data = adas_dg;
  up->minLogM0 = logNmin;
  up->minLogTe = logTmin;
  up->maxLogM0 = logNmax;
  up->maxLogTe = logTmax;
  up->dlogTe = tn_grid.dx[0];
  up->dlogM0 = tn_grid.dx[1];
  up->resTe = tn_grid.cells[0];
  up->resM0 = tn_grid.cells[1];
  up->adas_rng = modal_range;

  if (use_gpu) {
    up->recomb_data = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->adas_basis.num_basis, modal_range_ext.volume);
    gkyl_array_copy(up->recomb_data, adas_dg);
  }
  else {
    up->recomb_data = gkyl_array_new(GKYL_DOUBLE, up->adas_basis.num_basis, modal_range_ext.volume);
    gkyl_array_copy(up->recomb_data, adas_dg);
  }
  
  up->on_dev = up; // CPU eqn obj points to itself

  gkyl_array_release(adas_nodal);
  gkyl_array_release(adas_dg);

  return up;
}

void gkyl_dg_recomb_coll(const struct gkyl_dg_recomb *up,
  const struct gkyl_array *prim_vars_elc, 
  struct gkyl_array *coef_recomb, struct gkyl_array *cflrate)
{
#ifdef GKYL_HAVE_CUDA
  if(gkyl_array_is_cu_dev(coef_recomb)) {
    return gkyl_dg_recomb_coll_cu(up, prim_vars_elc, coef_recomb, cflrate);
  }
#endif

  struct gkyl_range_iter conf_iter, vel_iter;
  int rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<up->conf_rng->ndim; ++d) rem_dir[d] = 1;
  gkyl_range_iter_init(&conf_iter, up->conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc = gkyl_range_idx(up->conf_rng, conf_iter.idx);
    long nc = coef_recomb->ncomp;

    const double *prim_vars_elc_d = gkyl_array_cfetch(prim_vars_elc, loc);
    double *coef_recomb_d = gkyl_array_fetch(coef_recomb, loc);

    //Find cell containing value of n,T
    double cell_av_fac = pow(1/sqrt(2),up->cdim);
    double m0_elc_av = prim_vars_elc_d[0]*cell_av_fac;
    double temp_elc_av = prim_vars_elc_d[2*nc]*cell_av_fac*up->mass_elc/up->elem_charge;
    double log_Te_av = log10(temp_elc_av);
    double log_m0_av = log10(m0_elc_av);
    double cell_val_t;
    double cell_val_m0;
    int m0_idx, t_idx;
    double cell_vals_2d[2];
    double cell_center; 

    if (log_Te_av < up->minLogTe) {
      t_idx=1;
      log_Te_av = up->minLogTe;
    }
    else if (log_Te_av > up->maxLogTe) {
      t_idx=up->resTe;
      log_Te_av = up->maxLogTe;
    }
    else t_idx = (log_Te_av - up->minLogTe)/(up->dlogTe)+1;
    cell_center = (t_idx - 0.5)*up->dlogTe + up->minLogTe;
    cell_vals_2d[0] = 2.0*(log_Te_av - cell_center)/up->dlogTe; // Te value on cell interval
      
    if (log_m0_av < up->minLogM0) {
      m0_idx=1;
      log_m0_av = up->minLogM0;
    }
    else if (log_m0_av > up->maxLogM0) {
      m0_idx=up->resM0;
      log_m0_av = up->maxLogM0;
    }
    else m0_idx = (log_m0_av - up->minLogM0)/(up->dlogM0)+1;
    cell_center = (m0_idx - 0.5)*up->dlogM0 + up->minLogM0;
    cell_vals_2d[1] = 2.0*(log_m0_av - cell_center)/up->dlogM0; // M0 value on cell interval

    if ((m0_elc_av <= 0.) || (temp_elc_av <= 0.)) {
      coef_recomb_d[0] = 0.0; 
    }
    else {
      double *recomb_dat_d = gkyl_array_fetch(up->recomb_data, gkyl_range_idx(&up->adas_rng, (int[2]) {t_idx,m0_idx}));
      double adas_eval = up->adas_basis.eval_expand(cell_vals_2d, recomb_dat_d);
      coef_recomb_d[0] = pow(10.0,adas_eval)/cell_av_fac;
    }
  }
  //gkyl_grid_sub_array_write(&s->grid, &s->local, react->coef_react[i], "coef_recomb.gkyl");

  // cfl calculation
  //struct gkyl_range vel_rng;
  /* gkyl_range_deflate(&vel_rng, up->phase_rng, rem_dir, conf_iter.idx); */
  /* gkyl_range_iter_no_split_init(&vel_iter, &vel_rng); */
  /* // cfl associated with reaction is a *phase space* cfl */
  /* // Need to loop over velocity space for each configuration space cell */
  /* // to get total cfl rate in each phase space cell */
  /* while (gkyl_range_iter_next(&vel_iter)) { */
  /*   long cfl_idx = gkyl_range_idx(&vel_rng, vel_iter.idx); */
  /*   double *cflrate_d = gkyl_array_fetch(cflrate, cfl_idx); */
  /*   cflrate_d[0] += cflr; // frequencies are additive */
  /* } */

}

void
gkyl_dg_recomb_release(struct gkyl_dg_recomb* up)
{
  gkyl_array_release(up->recomb_data);
  free(up);
}
