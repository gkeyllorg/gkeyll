/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_iz.h>
#include <gkyl_dg_iz_priv.h>
#include <gkyl_util.h>
#include <gkyl_const.h>
#include <gkyl_array_ops_priv.h>
}

__global__ static void
gkyl_iz_react_rate_cu_ker(const struct gkyl_dg_iz *up, 
  const struct gkyl_range conf_rng, const struct gkyl_range adas_rng,
  const struct gkyl_basis *adas_basis, const struct gkyl_array* maxwellian_moms_elc, 
  struct gkyl_array *vtSq_iz1, struct gkyl_array *vtSq_iz2, struct gkyl_array* coef_iz,
  struct gkyl_array* ioniz_data, enum gkyl_react_self_type type_self, 
  double mass_elc, double elem_charge, double E, 
  double maxLogTe, double minLogTe, double dlogTe, int resTe, 
  double maxLogM0, double minLogM0, double dlogM0, int resM0)
{
  int cidx[GKYL_MAX_CDIM];
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_rng.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_rng, tid, cidx);
    long loc = gkyl_range_idx(&conf_rng, cidx);
    long nc = coef_iz->ncomp;

    int cdim = conf_rng.ndim;
    const double *maxwellian_moms_elc_d = (const double*) gkyl_array_cfetch(maxwellian_moms_elc, loc);  
    double *vtSq_iz1_d = (double*) gkyl_array_fetch(vtSq_iz1, loc);
    double *vtSq_iz2_d = (double*) gkyl_array_fetch(vtSq_iz2, loc);
    double *coef_iz_d = (double*) gkyl_array_fetch(coef_iz, loc);
    
    //Find nearest neighbor for n, Te in ADAS interpolated data
    double cell_av_fac = pow(1.0/sqrt(2.0),cdim);
    double m0_elc_av = maxwellian_moms_elc_d[0]*cell_av_fac;
    double temp_elc_av = maxwellian_moms_elc_d[2*nc]*cell_av_fac*mass_elc/elem_charge;
    double log_Te_av = log10(temp_elc_av);
    double log_m0_av = log10(m0_elc_av);
    int m0_idx, t_idx;
    double cell_vals_2d[2];
    double cell_center;
    double temp_elc_2;
    double temp_flr = 3.0; 
    
    if (log_Te_av < minLogTe) {
      t_idx=1;
      log_Te_av = minLogTe;
    }
    else if (log_Te_av > maxLogTe) {
      t_idx=resTe;
      log_Te_av = maxLogTe;
    }
    else t_idx = (log_Te_av - minLogTe)/(dlogTe)+1;
    cell_center = (t_idx - 0.5)*dlogTe + minLogTe;
    cell_vals_2d[0] = 2.0*(log_Te_av - cell_center)/dlogTe; // Te value on cell interval
      
    if (log_m0_av < minLogM0) {
      m0_idx=1;
      log_m0_av = minLogM0;
    }
    else if (log_m0_av > maxLogM0) {
      m0_idx=resM0;
      log_m0_av = maxLogM0;
    }
    else m0_idx = (log_m0_av - minLogM0)/(dlogM0)+1;
    cell_center = (m0_idx - 0.5)*dlogM0 + minLogM0;
    cell_vals_2d[1] = 2.0*(log_m0_av - cell_center)/dlogM0; // M0 value on cell interval
 
    if ((temp_elc_av <= 0.) || (m0_elc_av <= 0.)) {
      coef_iz_d[0] = 0.0;
    }
    else {
      int ad_idx[2] = {t_idx, m0_idx};
      double *iz_dat_d = (double*) gkyl_array_fetch(ioniz_data, gkyl_range_idx(&adas_rng, ad_idx));
      double adas_eval = adas_basis->eval_expand(cell_vals_2d, iz_dat_d);
      coef_iz_d[0] = pow(10.0,adas_eval)/cell_av_fac;

      if (type_self == GKYL_SELF_ELC) {
        //calculate vtSq_iz at each cell for primary and secondary elc
        if ( 3./2.*temp_elc_av >= 2.*E) {
          // T_e2 = 1/3*sqrt(3./2.*T_e*E_iz - E_iz)
          temp_elc_2 = 1.0/3.0 * pow(E*(3./2.*temp_elc_av - E), 0.5);
          vtSq_iz2_d[0] = temp_elc_2*elem_charge/(mass_elc*cell_av_fac);
          // T_e1 = T_e - 2/3*E_iz - T_e2
          array_set2(nc, 2*nc, vtSq_iz1_d, 1.0, maxwellian_moms_elc_d);
          vtSq_iz1_d[0] = vtSq_iz1_d[0] - elem_charge*(2./3.*E + temp_elc_2)/(mass_elc*cell_av_fac);
        }
        else if (3./2.*temp_elc_av >= E) {
          // T_e2 = 1/3*(3/2*T_e - E_iz)
          temp_elc_2 = temp_elc_av/2.0 - E/3.0;
          vtSq_iz2_d[0] = temp_elc_2*elem_charge/(mass_elc*cell_av_fac);
          // T_e1 = T_e - 2/3*E_iz - T_e2
          array_set2(nc, 2*nc, vtSq_iz1_d, 1.0, maxwellian_moms_elc_d);
          vtSq_iz1_d[0] = vtSq_iz1_d[0] - elem_charge*(2./3.*E + temp_elc_2)/(mass_elc*cell_av_fac);
        }
        else {
          vtSq_iz2_d[0] = temp_flr*elem_charge/(mass_elc*cell_av_fac);
          vtSq_iz1_d[0] = temp_flr*elem_charge/(mass_elc*cell_av_fac);	  
        }
      }
    }
  }
}

void gkyl_dg_iz_coll_cu(const struct gkyl_dg_iz *up, 
  const struct gkyl_array *maxwellian_moms_elc, 
  struct gkyl_array *vtSq_iz1, struct gkyl_array *vtSq_iz2,
  struct gkyl_array *coef_iz, struct gkyl_array *cflrate) 
{
  gkyl_iz_react_rate_cu_ker<<<up->conf_rng->nblocks, up->conf_rng->nthreads>>>(up->on_dev, 
    *up->conf_rng, up->adas_rng, up->basis_on_dev, maxwellian_moms_elc->on_dev, 
    vtSq_iz1->on_dev, vtSq_iz2->on_dev, coef_iz->on_dev, 
    up->ioniz_data->on_dev, up->type_self, up->mass_elc, up->elem_charge, up->E, 
    up->maxLogTe, up->minLogTe, up->dlogTe, up->resTe, 
    up->maxLogM0, up->minLogM0, up->dlogM0, up->resM0);
}
