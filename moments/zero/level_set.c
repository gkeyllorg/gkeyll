#include <gkyl_level_set.h>

#include <gkyl_wv_euler_rgfm_priv.h>
#include <gkyl_wv_gr_maxwell_priv.h>
#include <gkyl_wv_gr_maxwell_tetrad_priv.h>
#include <gkyl_wv_gr_euler_priv.h>
#include <gkyl_wv_gr_euler_tetrad_priv.h>
#include <gkyl_wv_gr_ultra_rel_euler_priv.h>
#include <gkyl_wv_gr_ultra_rel_euler_tetrad_priv.h>
#include <gkyl_wv_gr_twofluid_priv.h>
#include <gkyl_wv_gr_twofluid_tetrad_priv.h>
#include <gkyl_wv_gr_mhd_priv.h>
#include <gkyl_wv_gr_mhd_tetrad_priv.h>
#include <gkyl_gr_blackhole.h>

struct gkyl_wave_prop {
  struct gkyl_rect_grid grid; // Grid object.
  int ndim; // Number of dimensions.
  int num_up_dirs; // Number of update directions.
  int update_dirs[GKYL_MAX_DIM]; // Directions to update.
  enum gkyl_wave_limiter limiter; // Limiter to use.
  double cfl; // CFL number.
  const struct gkyl_wv_eqn *equation; // Equation object.

  bool force_low_order_flux; // Only use Lax flux.
  bool check_inv_domain; // Flag to indicate if invariant domains are checked.

  enum gkyl_wave_split_type split_type; // Type of splitting to use.

  struct gkyl_wave_geom *geom; // Geometry object.
  struct gkyl_comm *comm; // Communicator.

  // Data for 1D slice update.
  struct gkyl_array *waves, *apdq, *amdq, *speeds, *flux2;
  // Flags to indicate if fluctuations should be recomputed.
  struct gkyl_array *redo_fluct;

  // Some stats.
  long n_calls; // Number of calls to updater.
  long n_bad_advance_calls; // Number of calls in which positivity had to be fixed.
  long n_bad_cells; // Number of cells fixed.
  long n_max_bad_cells; // Maximum number of cells fixed in a call.
};

void
euler_rgfm_reinit_level_set(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir)
{
  const struct gkyl_wv_eqn* eqn = wv->equation;
  const struct wv_euler_rgfm *euler_rgfm = container_of(eqn, struct wv_euler_rgfm, eqn);
  int num_species = euler_rgfm->num_species;
  int reinit_freq = euler_rgfm->reinit_freq;

  for (int i = loidx_c; i <= upidx_c; i++) {
    idxl[dir] = i;

    double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));

    double reinit_param = qnew[4 + (2 * num_species)];

    if (reinit_param > reinit_freq) {
      double rho_total = qnew[0];

      bool update_up = false;
      bool update_down = false;
      for (int j = 0; j < num_species - 1; j++) {
        if (qnew[5 + j] / rho_total >= 0.5) {
          idxl[dir] = i - 1;
          double *ql = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
          idxl[dir] = i - 2;
          double *qll = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
          idxl[dir] = i - 3;
          double *qlll = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
        
          double rho_total_l = ql[0];
          double rho_total_ll = qll[0];
          double rho_total_lll = qlll[0];

          idxl[dir] = i + 1;
          double *qr = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
          idxl[dir] = i + 2;
          double *qrr = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
          idxl[dir] = i + 3;
          double *qrrr = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
        
          double rho_total_r = qr[0];
          double rho_total_rr = qrr[0];
          double rho_total_rrr = qrrr[0];
        
          if (ql[5 + j] / rho_total_l < 0.5 || qll[5 + j] / rho_total_ll < 0.5 || qlll[5 + j] / rho_total_lll < 0.5 || qr[5 + j] / rho_total_r < 0.5 ||
            qrr[5 + j] / rho_total_rr < 0.5 || qrrr[5 + j] / rho_total_rrr < 0.5)  {
            qnew[5 + j] = 0.99999 * rho_total;
            qnew[4 + num_species + j] = 0.99999 * rho_total;
            update_up = true;
          }
        }
        
        if (qnew[5 + j] / rho_total < 0.5) {
          idxl[dir] = i - 1;
          double *ql = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
          idxl[dir] = i - 2;
          double *qll = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
          idxl[dir] = i - 3;
          double *qlll = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
        
          double rho_total_l = ql[0];
          double rho_total_ll = qll[0];
          double rho_total_lll = qlll[0];

          idxl[dir] = i + 1;
          double *qr = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
          idxl[dir] = i + 2;
          double *qrr = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
          idxl[dir] = i + 3;
          double *qrrr = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
        
          double rho_total_r = qr[0];
          double rho_total_rr = qrr[0];
          double rho_total_rrr = qrrr[0];

          if (qr[5 + j] / rho_total_r >= 0.5 || qrr[5 + j] / rho_total_rr >= 0.5 || qrrr[5 + j] / rho_total_rrr >= 0.5 || ql[5 + j] / rho_total_l >= 0.5 ||
            qll[5 + j] / rho_total_ll >= 0.5 || qlll[5 + j] / rho_total_lll >= 0.5) {
            qnew[5 + j] = 0.00001 * rho_total;
            qnew[4 + num_species + j] = 0.00001 * rho_total;
            update_down = true;
          }
        }
      }

      if (update_up) {
        qnew[3 + (2 * num_species)] = 0.00001 * rho_total;
      }
      if (update_down) {
        qnew[3 + (2 * num_species)] = 0.99999 * rho_total;
      }

      qnew[4 + (2 * num_species)] = 0.0;
    }
    else {
      qnew[4 + (2 * num_species)] += 1.0;
    }
  }
}

void
gr_maxwell_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir)
{
  const struct gkyl_wv_eqn* eqn = wv->equation;
  const struct wv_gr_maxwell *gr_maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);
  
  const enum gkyl_spacetime_gauge spacetime_gauge = gr_maxwell->spacetime_gauge;

  if (spacetime_gauge == GKYL_STATIC_GAUGE) {
    const struct gkyl_gr_spacetime* spacetime = gr_maxwell->spacetime;
    int reinit_freq = gr_maxwell->reinit_freq;
    
    for (int i = loidx_c; i<= upidx_c; i++) {
      idxl[dir] = i;

      double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
      double evol_param = qnew[22];

      if (evol_param > reinit_freq) {
        double x = qnew[23];
        double y = qnew[24];
        double z = qnew[25];

        double lapse;
        double *shift = gkyl_malloc(sizeof(double[3]));
        bool in_excision_region;

        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
        }

        spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
        spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
        spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);

        spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);

        qnew[8] = lapse;
        qnew[9] = shift[0]; qnew[10] = shift[1]; qnew[11] = shift[2];

        qnew[12] = spatial_metric[0][0]; qnew[13] = spatial_metric[0][1]; qnew[14] = spatial_metric[0][2];
        qnew[15] = spatial_metric[1][0]; qnew[16] = spatial_metric[1][1]; qnew[17] = spatial_metric[1][2];
        qnew[18] = spatial_metric[2][0]; qnew[19] = spatial_metric[2][1]; qnew[20] = spatial_metric[2][2];

        if (in_excision_region) {
          for (int i = 0; i < 22; i++) {
            qnew[i] = 0.0;
          }

          qnew[21] = -1.0;
        }
        else {
          qnew[21] = 1.0;
        }

        for (int i = 0; i < 3; i++) {
          gkyl_free(spatial_metric[i]);
        }
        gkyl_free(spatial_metric);
        gkyl_free(shift);
        
        qnew[22] = 0.0;
      }
      else {
        qnew[22] += 1.0;
      }
    }
  }
}

void
gr_maxwell_tetrad_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir)
{
  const struct gkyl_wv_eqn* eqn = wv->equation;
  const struct wv_gr_maxwell_tetrad *gr_maxwell_tetrad = container_of(eqn, struct wv_gr_maxwell_tetrad, eqn);
  
  const enum gkyl_spacetime_gauge spacetime_gauge = gr_maxwell_tetrad->spacetime_gauge;

  if (spacetime_gauge == GKYL_STATIC_GAUGE) {
    const struct gkyl_gr_spacetime* spacetime = gr_maxwell_tetrad->spacetime;
    int reinit_freq = gr_maxwell_tetrad->reinit_freq;
    
    for (int i = loidx_c; i<= upidx_c; i++) {
      idxl[dir] = i;

      double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
      double evol_param = qnew[22];

      if (evol_param > reinit_freq) {
        double x = qnew[23];
        double y = qnew[24];
        double z = qnew[25];

        double lapse;
        double *shift = gkyl_malloc(sizeof(double[3]));
        bool in_excision_region;

        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
        }

        spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
        spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
        spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);

        spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);

        qnew[8] = lapse;
        qnew[9] = shift[0]; qnew[10] = shift[1]; qnew[11] = shift[2];

        qnew[12] = spatial_metric[0][0]; qnew[13] = spatial_metric[0][1]; qnew[14] = spatial_metric[0][2];
        qnew[15] = spatial_metric[1][0]; qnew[16] = spatial_metric[1][1]; qnew[17] = spatial_metric[1][2];
        qnew[18] = spatial_metric[2][0]; qnew[19] = spatial_metric[2][1]; qnew[20] = spatial_metric[2][2];

        if (in_excision_region) {
          for (int i = 0; i < 22; i++) {
            qnew[i] = 0.0;
          }

          qnew[21] = -1.0;
        }
        else {
          qnew[21] = 1.0;
        }

        for (int i = 0; i < 3; i++) {
          gkyl_free(spatial_metric[i]);
        }
        gkyl_free(spatial_metric);
        gkyl_free(shift);
        
        qnew[22] = 0.0;
      }
      else {
        qnew[22] += 1.0;
      }
    }
  }
}

void
gr_euler_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir)
{
  const struct gkyl_wv_eqn* eqn = wv->equation;
  const struct wv_gr_euler *gr_euler = container_of(eqn, struct wv_gr_euler, eqn);
  
  const enum gkyl_spacetime_gauge spacetime_gauge = gr_euler->spacetime_gauge;

  if (spacetime_gauge == GKYL_STATIC_GAUGE) {
    const struct gkyl_gr_spacetime* spacetime = gr_euler->spacetime;
    int reinit_freq = gr_euler->reinit_freq;
    
    for (int i = loidx_c; i<= upidx_c; i++) {
      idxl[dir] = i;

      double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
      double evol_param = qnew[67];

      if (evol_param > reinit_freq) {
        double x = qnew[68];
        double y = qnew[69];
        double z = qnew[70];

        double lapse;
        double *shift = gkyl_malloc(sizeof(double[3]));
        bool in_excision_region;

        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
        }

        double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
        }

        double *lapse_der = gkyl_malloc(sizeof(double[3]));
        double **shift_der = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          shift_der[i] = gkyl_malloc(sizeof(double[3]));
        }

        double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
        for (int i = 0; i < 3; i++) {
          spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
        spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
        spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);

        spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
        spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

        spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
        spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
        spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

        qnew[5] = lapse;
        qnew[6] = shift[0]; qnew[7] = shift[1]; qnew[8] = shift[2];

        qnew[9] = spatial_metric[0][0]; qnew[10] = spatial_metric[0][1]; qnew[11] = spatial_metric[0][2];
        qnew[12] = spatial_metric[1][0]; qnew[13] = spatial_metric[1][1]; qnew[14] = spatial_metric[1][2];
        qnew[15] = spatial_metric[2][0]; qnew[16] = spatial_metric[2][1]; qnew[17] = spatial_metric[2][2];

        qnew[18] = extrinsic_curvature[0][0]; qnew[19] = extrinsic_curvature[0][1]; qnew[20] = extrinsic_curvature[0][2];
        qnew[21] = extrinsic_curvature[1][0]; qnew[22] = extrinsic_curvature[1][1]; qnew[23] = extrinsic_curvature[1][2];
        qnew[24] = extrinsic_curvature[2][0]; qnew[25] = extrinsic_curvature[2][1]; qnew[26] = extrinsic_curvature[2][2];

        if (in_excision_region) {
          qnew[27] = -1.0;
        }
        else {
          qnew[27] = 1.0;
        }

        qnew[28] = lapse_der[0]; qnew[29] = lapse_der[1]; qnew[30] = lapse_der[2];

        qnew[31] = shift_der[0][0]; qnew[32] = shift_der[0][1]; qnew[33] = shift_der[0][2];
        qnew[34] = shift_der[1][0]; qnew[35] = shift_der[1][1]; qnew[36] = shift_der[1][2];
        qnew[37] = shift_der[2][0]; qnew[38] = shift_der[2][1]; qnew[39] = shift_der[2][2];

        qnew[40] = spatial_metric_der[0][0][0]; qnew[41] = spatial_metric_der[0][0][1]; qnew[42] = spatial_metric_der[0][0][2];
        qnew[43] = spatial_metric_der[0][1][0]; qnew[44] = spatial_metric_der[0][1][1]; qnew[45] = spatial_metric_der[0][1][2];
        qnew[46] = spatial_metric_der[0][2][0]; qnew[47] = spatial_metric_der[0][2][1]; qnew[48] = spatial_metric_der[0][2][2];

        qnew[49] = spatial_metric_der[1][0][0]; qnew[50] = spatial_metric_der[1][0][1]; qnew[51] = spatial_metric_der[1][0][2];
        qnew[52] = spatial_metric_der[1][1][0]; qnew[53] = spatial_metric_der[1][1][1]; qnew[54] = spatial_metric_der[1][1][2];
        qnew[55] = spatial_metric_der[1][2][0]; qnew[56] = spatial_metric_der[1][2][1]; qnew[57] = spatial_metric_der[1][2][2];

        qnew[58] = spatial_metric_der[2][0][0]; qnew[59] = spatial_metric_der[2][0][1]; qnew[60] = spatial_metric_der[2][0][2];
        qnew[61] = spatial_metric_der[2][1][0]; qnew[62] = spatial_metric_der[2][1][1]; qnew[63] = spatial_metric_der[2][1][2];
        qnew[64] = spatial_metric_der[2][2][0]; qnew[65] = spatial_metric_der[2][2][1]; qnew[66] = spatial_metric_der[2][2][2];

        if (in_excision_region) {
          for (int i = 0; i < 67; i++) {
            qnew[i] = 0.0;
          }
          qnew[27] = -1.0;
        }

        for (int i = 0; i < 3; i++) {
          gkyl_free(spatial_metric[i]);
          gkyl_free(extrinsic_curvature[i]);
          gkyl_free(shift_der[i]);
      
          for (int j = 0; j < 3; j++) {
            gkyl_free(spatial_metric_der[i][j]);
          }
          gkyl_free(spatial_metric_der[i]);
        }
        gkyl_free(spatial_metric);
        gkyl_free(extrinsic_curvature);
        gkyl_free(shift);
        gkyl_free(lapse_der);
        gkyl_free(shift_der);
        gkyl_free(spatial_metric_der);
        
        qnew[67] = 0.0;
      }
      else {
        qnew[67] += 1.0;
      }
    }
  }
}

void
gr_euler_tetrad_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir)
{
  const struct gkyl_wv_eqn* eqn = wv->equation;
  const struct wv_gr_euler_tetrad *gr_euler_tetrad = container_of(eqn, struct wv_gr_euler_tetrad, eqn);
  
  const enum gkyl_spacetime_gauge spacetime_gauge = gr_euler_tetrad->spacetime_gauge;

  if (spacetime_gauge == GKYL_STATIC_GAUGE) {
    const struct gkyl_gr_spacetime* spacetime = gr_euler_tetrad->spacetime;
    int reinit_freq = gr_euler_tetrad->reinit_freq;
    
    for (int i = loidx_c; i<= upidx_c; i++) {
      idxl[dir] = i;

      double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
      double evol_param = qnew[67];

      if (evol_param > reinit_freq) {
        double x = qnew[68];
        double y = qnew[69];
        double z = qnew[70];

        double lapse;
        double *shift = gkyl_malloc(sizeof(double[3]));
        bool in_excision_region;

        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
        }

        double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
        }

        double *lapse_der = gkyl_malloc(sizeof(double[3]));
        double **shift_der = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          shift_der[i] = gkyl_malloc(sizeof(double[3]));
        }

        double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
        for (int i = 0; i < 3; i++) {
          spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
        spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
        spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);

        spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
        spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

        spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
        spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
        spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

        qnew[5] = lapse;
        qnew[6] = shift[0]; qnew[7] = shift[1]; qnew[8] = shift[2];

        qnew[9] = spatial_metric[0][0]; qnew[10] = spatial_metric[0][1]; qnew[11] = spatial_metric[0][2];
        qnew[12] = spatial_metric[1][0]; qnew[13] = spatial_metric[1][1]; qnew[14] = spatial_metric[1][2];
        qnew[15] = spatial_metric[2][0]; qnew[16] = spatial_metric[2][1]; qnew[17] = spatial_metric[2][2];

        qnew[18] = extrinsic_curvature[0][0]; qnew[19] = extrinsic_curvature[0][1]; qnew[20] = extrinsic_curvature[0][2];
        qnew[21] = extrinsic_curvature[1][0]; qnew[22] = extrinsic_curvature[1][1]; qnew[23] = extrinsic_curvature[1][2];
        qnew[24] = extrinsic_curvature[2][0]; qnew[25] = extrinsic_curvature[2][1]; qnew[26] = extrinsic_curvature[2][2];

        if (in_excision_region) {
          qnew[27] = -1.0;
        }
        else {
          qnew[27] = 1.0;
        }

        qnew[28] = lapse_der[0]; qnew[29] = lapse_der[1]; qnew[30] = lapse_der[2];

        qnew[31] = shift_der[0][0]; qnew[32] = shift_der[0][1]; qnew[33] = shift_der[0][2];
        qnew[34] = shift_der[1][0]; qnew[35] = shift_der[1][1]; qnew[36] = shift_der[1][2];
        qnew[37] = shift_der[2][0]; qnew[38] = shift_der[2][1]; qnew[39] = shift_der[2][2];

        qnew[40] = spatial_metric_der[0][0][0]; qnew[41] = spatial_metric_der[0][0][1]; qnew[42] = spatial_metric_der[0][0][2];
        qnew[43] = spatial_metric_der[0][1][0]; qnew[44] = spatial_metric_der[0][1][1]; qnew[45] = spatial_metric_der[0][1][2];
        qnew[46] = spatial_metric_der[0][2][0]; qnew[47] = spatial_metric_der[0][2][1]; qnew[48] = spatial_metric_der[0][2][2];

        qnew[49] = spatial_metric_der[1][0][0]; qnew[50] = spatial_metric_der[1][0][1]; qnew[51] = spatial_metric_der[1][0][2];
        qnew[52] = spatial_metric_der[1][1][0]; qnew[53] = spatial_metric_der[1][1][1]; qnew[54] = spatial_metric_der[1][1][2];
        qnew[55] = spatial_metric_der[1][2][0]; qnew[56] = spatial_metric_der[1][2][1]; qnew[57] = spatial_metric_der[1][2][2];

        qnew[58] = spatial_metric_der[2][0][0]; qnew[59] = spatial_metric_der[2][0][1]; qnew[60] = spatial_metric_der[2][0][2];
        qnew[61] = spatial_metric_der[2][1][0]; qnew[62] = spatial_metric_der[2][1][1]; qnew[63] = spatial_metric_der[2][1][2];
        qnew[64] = spatial_metric_der[2][2][0]; qnew[65] = spatial_metric_der[2][2][1]; qnew[66] = spatial_metric_der[2][2][2];

        if (in_excision_region) {
          for (int i = 0; i < 67; i++) {
            qnew[i] = 0.0;
          }
          qnew[27] = -1.0;
        }

        for (int i = 0; i < 3; i++) {
          gkyl_free(spatial_metric[i]);
          gkyl_free(extrinsic_curvature[i]);
          gkyl_free(shift_der[i]);
      
          for (int j = 0; j < 3; j++) {
            gkyl_free(spatial_metric_der[i][j]);
          }
          gkyl_free(spatial_metric_der[i]);
        }
        gkyl_free(spatial_metric);
        gkyl_free(extrinsic_curvature);
        gkyl_free(shift);
        gkyl_free(lapse_der);
        gkyl_free(shift_der);
        gkyl_free(spatial_metric_der);
        
        qnew[67] = 0.0;
      }
      else {
        qnew[67] += 1.0;
      }
    }
  }
}

void
gr_ultra_rel_euler_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir)
{
  const struct gkyl_wv_eqn* eqn = wv->equation;
  const struct wv_gr_ultra_rel_euler *gr_ultra_rel_euler = container_of(eqn, struct wv_gr_ultra_rel_euler, eqn);

  const enum gkyl_spacetime_gauge spacetime_gauge = gr_ultra_rel_euler->spacetime_gauge;

  if (spacetime_gauge == GKYL_STATIC_GAUGE) {
    const struct gkyl_gr_spacetime* spacetime = gr_ultra_rel_euler->spacetime;
    int reinit_freq = gr_ultra_rel_euler->reinit_freq;
    
    for (int i = loidx_c; i<= upidx_c; i++) {
      idxl[dir] = i;

      double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
      double evol_param = qnew[66];

      if (evol_param > reinit_freq) {
        double x = qnew[67];
        double y = qnew[68];
        double z = qnew[69];

        double lapse;
        double *shift = gkyl_malloc(sizeof(double[3]));
        bool in_excision_region;

        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
        }

        double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
        }

        double *lapse_der = gkyl_malloc(sizeof(double[3]));
        double **shift_der = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          shift_der[i] = gkyl_malloc(sizeof(double[3]));
        }

        double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
        for (int i = 0; i < 3; i++) {
          spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
        spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
        spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);

        spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
        spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

        spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
        spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
        spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

        qnew[4] = lapse;
        qnew[5] = shift[0]; qnew[6] = shift[1]; qnew[7] = shift[2];

        qnew[8] = spatial_metric[0][0]; qnew[9] = spatial_metric[0][1]; qnew[10] = spatial_metric[0][2];
        qnew[11] = spatial_metric[1][0]; qnew[12] = spatial_metric[1][1]; qnew[13] = spatial_metric[1][2];
        qnew[14] = spatial_metric[2][0]; qnew[15] = spatial_metric[2][1]; qnew[16] = spatial_metric[2][2];

        qnew[17] = extrinsic_curvature[0][0]; qnew[18] = extrinsic_curvature[0][1]; qnew[19] = extrinsic_curvature[0][2];
        qnew[20] = extrinsic_curvature[1][0]; qnew[21] = extrinsic_curvature[1][1]; qnew[22] = extrinsic_curvature[1][2];
        qnew[23] = extrinsic_curvature[2][0]; qnew[24] = extrinsic_curvature[2][1]; qnew[25] = extrinsic_curvature[2][2];

        if (in_excision_region) {
          qnew[26] = -1.0;
        }
        else {
          qnew[26] = 1.0;
        }

        qnew[27] = lapse_der[0]; qnew[28] = lapse_der[1]; qnew[29] = lapse_der[2];

        qnew[30] = shift_der[0][0]; qnew[31] = shift_der[0][1]; qnew[32] = shift_der[0][2];
        qnew[33] = shift_der[1][0]; qnew[34] = shift_der[1][1]; qnew[35] = shift_der[1][2];
        qnew[36] = shift_der[2][0]; qnew[37] = shift_der[2][1]; qnew[38] = shift_der[2][2];

        qnew[39] = spatial_metric_der[0][0][0]; qnew[40] = spatial_metric_der[0][0][1]; qnew[41] = spatial_metric_der[0][0][2];
        qnew[42] = spatial_metric_der[0][1][0]; qnew[43] = spatial_metric_der[0][1][1]; qnew[44] = spatial_metric_der[0][1][2];
        qnew[45] = spatial_metric_der[0][2][0]; qnew[46] = spatial_metric_der[0][2][1]; qnew[47] = spatial_metric_der[0][2][2];

        qnew[48] = spatial_metric_der[1][0][0]; qnew[49] = spatial_metric_der[1][0][1]; qnew[50] = spatial_metric_der[1][0][2];
        qnew[51] = spatial_metric_der[1][1][0]; qnew[52] = spatial_metric_der[1][1][1]; qnew[53] = spatial_metric_der[1][1][2];
        qnew[54] = spatial_metric_der[1][2][0]; qnew[55] = spatial_metric_der[1][2][1]; qnew[56] = spatial_metric_der[1][2][2];

        qnew[57] = spatial_metric_der[2][0][0]; qnew[58] = spatial_metric_der[2][0][1]; qnew[59] = spatial_metric_der[2][0][2];
        qnew[60] = spatial_metric_der[2][1][0]; qnew[61] = spatial_metric_der[2][1][1]; qnew[62] = spatial_metric_der[2][1][2];
        qnew[63] = spatial_metric_der[2][2][0]; qnew[64] = spatial_metric_der[2][2][1]; qnew[65] = spatial_metric_der[2][2][2];

        if (in_excision_region) {
          for (int i = 0; i < 66; i++) {
            qnew[i] = 0.0;
          }
          qnew[26] = -1.0;
        }

        for (int i = 0; i < 3; i++) {
          gkyl_free(spatial_metric[i]);
          gkyl_free(extrinsic_curvature[i]);
          gkyl_free(shift_der[i]);
      
          for (int j = 0; j < 3; j++) {
            gkyl_free(spatial_metric_der[i][j]);
          }
          gkyl_free(spatial_metric_der[i]);
        }
        gkyl_free(spatial_metric);
        gkyl_free(extrinsic_curvature);
        gkyl_free(shift);
        gkyl_free(lapse_der);
        gkyl_free(shift_der);
        gkyl_free(spatial_metric_der);

        qnew[66] = 0.0;
      }
      else {
        qnew[66] += 1.0;
      }
    }
  }
  else if (spacetime_gauge == GKYL_BLACKHOLE_COLLAPSE_GAUGE) {
    const struct gkyl_gr_spacetime* spacetime = gr_ultra_rel_euler->spacetime;
    const struct gr_blackhole *blackhole = container_of(spacetime, struct gr_blackhole, spacetime);

    double mass = blackhole->mass;
    double spin = blackhole->spin;

    double pos_x = blackhole->pos_x;
    double pos_y = blackhole->pos_y;
    double pos_z = blackhole->pos_z;

    for (int i = loidx_c; i <= upidx_c; i++) {
      idxl[dir] = i;

      double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
      double evol_param = qnew[66];

      if (evol_param < mass) {
        evol_param += 0.001 + (0.001 * evol_param);
      }
      qnew[66] = evol_param;

      double x = qnew[67];
      double y = qnew[68];
      double z = qnew[69];

      struct gkyl_gr_spacetime *new_spacetime = gkyl_gr_blackhole_new(false, fmin(evol_param, mass), spin, pos_x, pos_y, pos_z);

      double lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der = gkyl_malloc(sizeof(double[3]));
      double **shift_der = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      new_spacetime->lapse_function_func(new_spacetime, 0.0, x, y, z, &lapse);
      new_spacetime->shift_vector_func(new_spacetime, 0.0, x, y, z, &shift);
      new_spacetime->excision_region_func(new_spacetime, 0.0, x, y, z, &in_excision_region);

      new_spacetime->spatial_metric_tensor_func(new_spacetime, 0.0, x, y, z, &spatial_metric);
      new_spacetime->extrinsic_curvature_tensor_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

      new_spacetime->lapse_function_der_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
      new_spacetime->shift_vector_der_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
      new_spacetime->spatial_metric_tensor_der_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

      qnew[4] = lapse;
      qnew[5] = shift[0]; qnew[6] = shift[1]; qnew[7] = shift[2];

      qnew[8] = spatial_metric[0][0]; qnew[9] = spatial_metric[0][1]; qnew[10] = spatial_metric[0][2];
      qnew[11] = spatial_metric[1][0]; qnew[12] = spatial_metric[1][1]; qnew[13] = spatial_metric[1][2];
      qnew[14] = spatial_metric[2][0]; qnew[15] = spatial_metric[2][1]; qnew[16] = spatial_metric[2][2];

      qnew[17] = extrinsic_curvature[0][0]; qnew[18] = extrinsic_curvature[0][1]; qnew[19] = extrinsic_curvature[0][2];
      qnew[20] = extrinsic_curvature[1][0]; qnew[21] = extrinsic_curvature[1][1]; qnew[22] = extrinsic_curvature[1][2];
      qnew[23] = extrinsic_curvature[2][0]; qnew[24] = extrinsic_curvature[2][1]; qnew[25] = extrinsic_curvature[2][2];

      if (in_excision_region) {
        qnew[26] = -1.0;
      }
      else {
        qnew[26] = 1.0;
      }

      qnew[27] = lapse_der[0]; qnew[28] = lapse_der[1]; qnew[29] = lapse_der[2];

      qnew[30] = shift_der[0][0]; qnew[31] = shift_der[0][1]; qnew[32] = shift_der[0][2];
      qnew[33] = shift_der[1][0]; qnew[34] = shift_der[1][1]; qnew[35] = shift_der[1][2];
      qnew[36] = shift_der[2][0]; qnew[37] = shift_der[2][1]; qnew[38] = shift_der[2][2];

      qnew[39] = spatial_metric_der[0][0][0]; qnew[40] = spatial_metric_der[0][0][1]; qnew[41] = spatial_metric_der[0][0][2];
      qnew[42] = spatial_metric_der[0][1][0]; qnew[43] = spatial_metric_der[0][1][1]; qnew[44] = spatial_metric_der[0][1][2];
      qnew[45] = spatial_metric_der[0][2][0]; qnew[46] = spatial_metric_der[0][2][1]; qnew[47] = spatial_metric_der[0][2][2];

      qnew[48] = spatial_metric_der[1][0][0]; qnew[49] = spatial_metric_der[1][0][1]; qnew[50] = spatial_metric_der[1][0][2];
      qnew[51] = spatial_metric_der[1][1][0]; qnew[52] = spatial_metric_der[1][1][1]; qnew[53] = spatial_metric_der[1][1][2];
      qnew[54] = spatial_metric_der[1][2][0]; qnew[55] = spatial_metric_der[1][2][1]; qnew[56] = spatial_metric_der[1][2][2];

      qnew[57] = spatial_metric_der[2][0][0]; qnew[58] = spatial_metric_der[2][0][1]; qnew[59] = spatial_metric_der[2][0][2];
      qnew[60] = spatial_metric_der[2][1][0]; qnew[61] = spatial_metric_der[2][1][1]; qnew[62] = spatial_metric_der[2][1][2];
      qnew[63] = spatial_metric_der[2][2][0]; qnew[64] = spatial_metric_der[2][2][1]; qnew[65] = spatial_metric_der[2][2][2];

      if (in_excision_region) {
        for (int i = 0; i < 66; i++) {
          qnew[i] = 0.0;
        }
        qnew[26] = -1.0;
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(extrinsic_curvature[i]);
        gkyl_free(shift_der[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der[i][j]);
        }
        gkyl_free(spatial_metric_der[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(extrinsic_curvature);
      gkyl_free(shift);
      gkyl_free(lapse_der);
      gkyl_free(shift_der);
      gkyl_free(spatial_metric_der);
      gkyl_free(new_spacetime);
    }
  }
}

void
gr_ultra_rel_euler_tetrad_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir)
{
  const struct gkyl_wv_eqn* eqn = wv->equation;
  const struct wv_gr_ultra_rel_euler_tetrad *gr_ultra_rel_euler_tetrad = container_of(eqn, struct wv_gr_ultra_rel_euler_tetrad, eqn);

  const enum gkyl_spacetime_gauge spacetime_gauge = gr_ultra_rel_euler_tetrad->spacetime_gauge;

  if (spacetime_gauge == GKYL_STATIC_GAUGE) {
    const struct gkyl_gr_spacetime* spacetime = gr_ultra_rel_euler_tetrad->spacetime;
    int reinit_freq = gr_ultra_rel_euler_tetrad->reinit_freq;
    
    for (int i = loidx_c; i<= upidx_c; i++) {
      idxl[dir] = i;

      double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
      double evol_param = qnew[66];

      if (evol_param > reinit_freq) {
        double x = qnew[67];
        double y = qnew[68];
        double z = qnew[69];

        double lapse;
        double *shift = gkyl_malloc(sizeof(double[3]));
        bool in_excision_region;

        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
        }

        double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
        }

        double *lapse_der = gkyl_malloc(sizeof(double[3]));
        double **shift_der = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          shift_der[i] = gkyl_malloc(sizeof(double[3]));
        }

        double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
        for (int i = 0; i < 3; i++) {
          spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
        spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
        spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);

        spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
        spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

        spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
        spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
        spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

        qnew[4] = lapse;
        qnew[5] = shift[0]; qnew[6] = shift[1]; qnew[7] = shift[2];

        qnew[8] = spatial_metric[0][0]; qnew[9] = spatial_metric[0][1]; qnew[10] = spatial_metric[0][2];
        qnew[11] = spatial_metric[1][0]; qnew[12] = spatial_metric[1][1]; qnew[13] = spatial_metric[1][2];
        qnew[14] = spatial_metric[2][0]; qnew[15] = spatial_metric[2][1]; qnew[16] = spatial_metric[2][2];

        qnew[17] = extrinsic_curvature[0][0]; qnew[18] = extrinsic_curvature[0][1]; qnew[19] = extrinsic_curvature[0][2];
        qnew[20] = extrinsic_curvature[1][0]; qnew[21] = extrinsic_curvature[1][1]; qnew[22] = extrinsic_curvature[1][2];
        qnew[23] = extrinsic_curvature[2][0]; qnew[24] = extrinsic_curvature[2][1]; qnew[25] = extrinsic_curvature[2][2];

        if (in_excision_region) {
          qnew[26] = -1.0;
        }
        else {
          qnew[26] = 1.0;
        }

        qnew[27] = lapse_der[0]; qnew[28] = lapse_der[1]; qnew[29] = lapse_der[2];

        qnew[30] = shift_der[0][0]; qnew[31] = shift_der[0][1]; qnew[32] = shift_der[0][2];
        qnew[33] = shift_der[1][0]; qnew[34] = shift_der[1][1]; qnew[35] = shift_der[1][2];
        qnew[36] = shift_der[2][0]; qnew[37] = shift_der[2][1]; qnew[38] = shift_der[2][2];

        qnew[39] = spatial_metric_der[0][0][0]; qnew[40] = spatial_metric_der[0][0][1]; qnew[41] = spatial_metric_der[0][0][2];
        qnew[42] = spatial_metric_der[0][1][0]; qnew[43] = spatial_metric_der[0][1][1]; qnew[44] = spatial_metric_der[0][1][2];
        qnew[45] = spatial_metric_der[0][2][0]; qnew[46] = spatial_metric_der[0][2][1]; qnew[47] = spatial_metric_der[0][2][2];

        qnew[48] = spatial_metric_der[1][0][0]; qnew[49] = spatial_metric_der[1][0][1]; qnew[50] = spatial_metric_der[1][0][2];
        qnew[51] = spatial_metric_der[1][1][0]; qnew[52] = spatial_metric_der[1][1][1]; qnew[53] = spatial_metric_der[1][1][2];
        qnew[54] = spatial_metric_der[1][2][0]; qnew[55] = spatial_metric_der[1][2][1]; qnew[56] = spatial_metric_der[1][2][2];

        qnew[57] = spatial_metric_der[2][0][0]; qnew[58] = spatial_metric_der[2][0][1]; qnew[59] = spatial_metric_der[2][0][2];
        qnew[60] = spatial_metric_der[2][1][0]; qnew[61] = spatial_metric_der[2][1][1]; qnew[62] = spatial_metric_der[2][1][2];
        qnew[63] = spatial_metric_der[2][2][0]; qnew[64] = spatial_metric_der[2][2][1]; qnew[65] = spatial_metric_der[2][2][2];

        if (in_excision_region) {
          for (int i = 0; i < 66; i++) {
            qnew[i] = 0.0;
          }
          qnew[26] = -1.0;
        }

        for (int i = 0; i < 3; i++) {
          gkyl_free(spatial_metric[i]);
          gkyl_free(extrinsic_curvature[i]);
          gkyl_free(shift_der[i]);
      
          for (int j = 0; j < 3; j++) {
            gkyl_free(spatial_metric_der[i][j]);
          }
          gkyl_free(spatial_metric_der[i]);
        }
        gkyl_free(spatial_metric);
        gkyl_free(extrinsic_curvature);
        gkyl_free(shift);
        gkyl_free(lapse_der);
        gkyl_free(shift_der);
        gkyl_free(spatial_metric_der);

        qnew[66] = 0.0;
      }
      else {
        qnew[66] += 1.0;
      }
    }
  }
  else if (spacetime_gauge == GKYL_BLACKHOLE_COLLAPSE_GAUGE) {
    const struct gkyl_gr_spacetime* spacetime = gr_ultra_rel_euler_tetrad->spacetime;
    const struct gr_blackhole *blackhole = container_of(spacetime, struct gr_blackhole, spacetime);

    double mass = blackhole->mass;
    double spin = blackhole->spin;

    double pos_x = blackhole->pos_x;
    double pos_y = blackhole->pos_y;
    double pos_z = blackhole->pos_z;

    for (int i = loidx_c; i <= upidx_c; i++) {
      idxl[dir] = i;

      double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
      double evol_param = qnew[66];

      if (evol_param < mass) {
        evol_param += 0.001 + (0.001 * evol_param);
      }
      qnew[66] = evol_param;

      double x = qnew[67];
      double y = qnew[68];
      double z = qnew[69];

      struct gkyl_gr_spacetime *new_spacetime = gkyl_gr_blackhole_new(false, fmin(evol_param, mass), spin, pos_x, pos_y, pos_z);

      double lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der = gkyl_malloc(sizeof(double[3]));
      double **shift_der = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      new_spacetime->lapse_function_func(new_spacetime, 0.0, x, y, z, &lapse);
      new_spacetime->shift_vector_func(new_spacetime, 0.0, x, y, z, &shift);
      new_spacetime->excision_region_func(new_spacetime, 0.0, x, y, z, &in_excision_region);

      new_spacetime->spatial_metric_tensor_func(new_spacetime, 0.0, x, y, z, &spatial_metric);
      new_spacetime->extrinsic_curvature_tensor_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

      new_spacetime->lapse_function_der_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
      new_spacetime->shift_vector_der_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
      new_spacetime->spatial_metric_tensor_der_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

      qnew[4] = lapse;
      qnew[5] = shift[0]; qnew[6] = shift[1]; qnew[7] = shift[2];

      qnew[8] = spatial_metric[0][0]; qnew[9] = spatial_metric[0][1]; qnew[10] = spatial_metric[0][2];
      qnew[11] = spatial_metric[1][0]; qnew[12] = spatial_metric[1][1]; qnew[13] = spatial_metric[1][2];
      qnew[14] = spatial_metric[2][0]; qnew[15] = spatial_metric[2][1]; qnew[16] = spatial_metric[2][2];

      qnew[17] = extrinsic_curvature[0][0]; qnew[18] = extrinsic_curvature[0][1]; qnew[19] = extrinsic_curvature[0][2];
      qnew[20] = extrinsic_curvature[1][0]; qnew[21] = extrinsic_curvature[1][1]; qnew[22] = extrinsic_curvature[1][2];
      qnew[23] = extrinsic_curvature[2][0]; qnew[24] = extrinsic_curvature[2][1]; qnew[25] = extrinsic_curvature[2][2];

      if (in_excision_region) {
        qnew[26] = -1.0;
      }
      else {
        qnew[26] = 1.0;
      }

      qnew[27] = lapse_der[0]; qnew[28] = lapse_der[1]; qnew[29] = lapse_der[2];

      qnew[30] = shift_der[0][0]; qnew[31] = shift_der[0][1]; qnew[32] = shift_der[0][2];
      qnew[33] = shift_der[1][0]; qnew[34] = shift_der[1][1]; qnew[35] = shift_der[1][2];
      qnew[36] = shift_der[2][0]; qnew[37] = shift_der[2][1]; qnew[38] = shift_der[2][2];

      qnew[39] = spatial_metric_der[0][0][0]; qnew[40] = spatial_metric_der[0][0][1]; qnew[41] = spatial_metric_der[0][0][2];
      qnew[42] = spatial_metric_der[0][1][0]; qnew[43] = spatial_metric_der[0][1][1]; qnew[44] = spatial_metric_der[0][1][2];
      qnew[45] = spatial_metric_der[0][2][0]; qnew[46] = spatial_metric_der[0][2][1]; qnew[47] = spatial_metric_der[0][2][2];

      qnew[48] = spatial_metric_der[1][0][0]; qnew[49] = spatial_metric_der[1][0][1]; qnew[50] = spatial_metric_der[1][0][2];
      qnew[51] = spatial_metric_der[1][1][0]; qnew[52] = spatial_metric_der[1][1][1]; qnew[53] = spatial_metric_der[1][1][2];
      qnew[54] = spatial_metric_der[1][2][0]; qnew[55] = spatial_metric_der[1][2][1]; qnew[56] = spatial_metric_der[1][2][2];

      qnew[57] = spatial_metric_der[2][0][0]; qnew[58] = spatial_metric_der[2][0][1]; qnew[59] = spatial_metric_der[2][0][2];
      qnew[60] = spatial_metric_der[2][1][0]; qnew[61] = spatial_metric_der[2][1][1]; qnew[62] = spatial_metric_der[2][1][2];
      qnew[63] = spatial_metric_der[2][2][0]; qnew[64] = spatial_metric_der[2][2][1]; qnew[65] = spatial_metric_der[2][2][2];

      if (in_excision_region) {
        for (int i = 0; i < 66; i++) {
          qnew[i] = 0.0;
        }
        qnew[26] = -1.0;
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(extrinsic_curvature[i]);
        gkyl_free(shift_der[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der[i][j]);
        }
        gkyl_free(spatial_metric_der[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(extrinsic_curvature);
      gkyl_free(shift);
      gkyl_free(lapse_der);
      gkyl_free(shift_der);
      gkyl_free(spatial_metric_der);
      gkyl_free(new_spacetime);
    }
  }
}

void
gr_twofluid_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir)
{
  const struct gkyl_wv_eqn* eqn = wv->equation;
  const struct wv_gr_twofluid *gr_twofluid = container_of(eqn, struct wv_gr_twofluid, eqn);
  
  const enum gkyl_spacetime_gauge spacetime_gauge = gr_twofluid->spacetime_gauge;

  if (spacetime_gauge == GKYL_STATIC_GAUGE) {
    const struct gkyl_gr_spacetime* spacetime = gr_twofluid->spacetime;
    int reinit_freq = gr_twofluid->reinit_freq;
    
    for (int i = loidx_c; i<= upidx_c; i++) {
      idxl[dir] = i;

      double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
      double evol_param = qnew[80];

      if (evol_param > reinit_freq) {
        double x = qnew[81];
        double y = qnew[82];
        double z = qnew[83];

        double lapse;
        double *shift = gkyl_malloc(sizeof(double[3]));
        bool in_excision_region;

        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
        }

        double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
        }

        double *lapse_der = gkyl_malloc(sizeof(double[3]));
        double **shift_der = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          shift_der[i] = gkyl_malloc(sizeof(double[3]));
        }

        double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
        for (int i = 0; i < 3; i++) {
          spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
        spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
        spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);

        spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
        spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

        spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
        spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
        spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

        qnew[18] = lapse;
        qnew[19] = shift[0]; qnew[20] = shift[1]; qnew[21] = shift[2];

        qnew[22] = spatial_metric[0][0]; qnew[23] = spatial_metric[0][1]; qnew[24] = spatial_metric[0][2];
        qnew[25] = spatial_metric[1][0]; qnew[26] = spatial_metric[1][1]; qnew[27] = spatial_metric[1][2];
        qnew[28] = spatial_metric[2][0]; qnew[29] = spatial_metric[2][1]; qnew[30] = spatial_metric[2][2];

        qnew[31] = extrinsic_curvature[0][0]; qnew[32] = extrinsic_curvature[0][1]; qnew[33] = extrinsic_curvature[0][2];
        qnew[34] = extrinsic_curvature[1][0]; qnew[35] = extrinsic_curvature[1][1]; qnew[36] = extrinsic_curvature[1][2];
        qnew[37] = extrinsic_curvature[2][0]; qnew[38] = extrinsic_curvature[2][1]; qnew[39] = extrinsic_curvature[2][2];

        if (in_excision_region) {
          qnew[40] = -1.0;
        }
        else {
          qnew[40] = 1.0;
        }

        qnew[41] = lapse_der[0]; qnew[42] = lapse_der[1]; qnew[43] = lapse_der[2];

        qnew[44] = shift_der[0][0]; qnew[45] = shift_der[0][1]; qnew[46] = shift_der[0][2];
        qnew[47] = shift_der[1][0]; qnew[48] = shift_der[1][1]; qnew[49] = shift_der[1][2];
        qnew[50] = shift_der[2][0]; qnew[51] = shift_der[2][1]; qnew[52] = shift_der[2][2];

        qnew[53] = spatial_metric_der[0][0][0]; qnew[54] = spatial_metric_der[0][0][1]; qnew[55] = spatial_metric_der[0][0][2];
        qnew[56] = spatial_metric_der[0][1][0]; qnew[57] = spatial_metric_der[0][1][1]; qnew[58] = spatial_metric_der[0][1][2];
        qnew[59] = spatial_metric_der[0][2][0]; qnew[60] = spatial_metric_der[0][2][1]; qnew[61] = spatial_metric_der[0][2][2];

        qnew[62] = spatial_metric_der[1][0][0]; qnew[63] = spatial_metric_der[1][0][1]; qnew[64] = spatial_metric_der[1][0][2];
        qnew[65] = spatial_metric_der[1][1][0]; qnew[66] = spatial_metric_der[1][1][1]; qnew[67] = spatial_metric_der[1][1][2];
        qnew[68] = spatial_metric_der[1][2][0]; qnew[69] = spatial_metric_der[1][2][1]; qnew[70] = spatial_metric_der[1][2][2];

        qnew[71] = spatial_metric_der[2][0][0]; qnew[72] = spatial_metric_der[2][0][1]; qnew[73] = spatial_metric_der[2][0][2];
        qnew[74] = spatial_metric_der[2][1][0]; qnew[75] = spatial_metric_der[2][1][1]; qnew[76] = spatial_metric_der[2][1][2];
        qnew[77] = spatial_metric_der[2][2][0]; qnew[78] = spatial_metric_der[2][2][1]; qnew[79] = spatial_metric_der[2][2][2];

        if (in_excision_region) {
          for (int i = 0; i < 80; i++) {
            qnew[i] = 0.0;
          }
          qnew[40] = -1.0;
        }

        for (int i = 0; i < 3; i++) {
          gkyl_free(spatial_metric[i]);
          gkyl_free(extrinsic_curvature[i]);
          gkyl_free(shift_der[i]);
      
          for (int j = 0; j < 3; j++) {
            gkyl_free(spatial_metric_der[i][j]);
          }
          gkyl_free(spatial_metric_der[i]);
        }
        gkyl_free(spatial_metric);
        gkyl_free(extrinsic_curvature);
        gkyl_free(shift);
        gkyl_free(lapse_der);
        gkyl_free(shift_der);
        gkyl_free(spatial_metric_der);
        
        qnew[80] = 0.0;
      }
      else {
        qnew[80] += 1.0;
      }
    }
  }
}

void
gr_twofluid_tetrad_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir)
{
  const struct gkyl_wv_eqn* eqn = wv->equation;
  const struct wv_gr_twofluid_tetrad *gr_twofluid_tetrad = container_of(eqn, struct wv_gr_twofluid_tetrad, eqn);
  
  const enum gkyl_spacetime_gauge spacetime_gauge = gr_twofluid_tetrad->spacetime_gauge;

  if (spacetime_gauge == GKYL_STATIC_GAUGE) {
    const struct gkyl_gr_spacetime* spacetime = gr_twofluid_tetrad->spacetime;
    int reinit_freq = gr_twofluid_tetrad->reinit_freq;
    
    for (int i = loidx_c; i<= upidx_c; i++) {
      idxl[dir] = i;

      double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
      double evol_param = qnew[80];

      if (evol_param > reinit_freq) {
        double x = qnew[81];
        double y = qnew[82];
        double z = qnew[83];

        double lapse;
        double *shift = gkyl_malloc(sizeof(double[3]));
        bool in_excision_region;

        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
        }

        double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
        }

        double *lapse_der = gkyl_malloc(sizeof(double[3]));
        double **shift_der = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          shift_der[i] = gkyl_malloc(sizeof(double[3]));
        }

        double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
        for (int i = 0; i < 3; i++) {
          spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
        spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
        spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);

        spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
        spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

        spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
        spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
        spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

        qnew[18] = lapse;
        qnew[19] = shift[0]; qnew[20] = shift[1]; qnew[21] = shift[2];

        qnew[22] = spatial_metric[0][0]; qnew[23] = spatial_metric[0][1]; qnew[24] = spatial_metric[0][2];
        qnew[25] = spatial_metric[1][0]; qnew[26] = spatial_metric[1][1]; qnew[27] = spatial_metric[1][2];
        qnew[28] = spatial_metric[2][0]; qnew[29] = spatial_metric[2][1]; qnew[30] = spatial_metric[2][2];

        qnew[31] = extrinsic_curvature[0][0]; qnew[32] = extrinsic_curvature[0][1]; qnew[33] = extrinsic_curvature[0][2];
        qnew[34] = extrinsic_curvature[1][0]; qnew[35] = extrinsic_curvature[1][1]; qnew[36] = extrinsic_curvature[1][2];
        qnew[37] = extrinsic_curvature[2][0]; qnew[38] = extrinsic_curvature[2][1]; qnew[39] = extrinsic_curvature[2][2];

        if (in_excision_region) {
          qnew[40] = -1.0;
        }
        else {
          qnew[40] = 1.0;
        }

        qnew[41] = lapse_der[0]; qnew[42] = lapse_der[1]; qnew[43] = lapse_der[2];

        qnew[44] = shift_der[0][0]; qnew[45] = shift_der[0][1]; qnew[46] = shift_der[0][2];
        qnew[47] = shift_der[1][0]; qnew[48] = shift_der[1][1]; qnew[49] = shift_der[1][2];
        qnew[50] = shift_der[2][0]; qnew[51] = shift_der[2][1]; qnew[52] = shift_der[2][2];

        qnew[53] = spatial_metric_der[0][0][0]; qnew[54] = spatial_metric_der[0][0][1]; qnew[55] = spatial_metric_der[0][0][2];
        qnew[56] = spatial_metric_der[0][1][0]; qnew[57] = spatial_metric_der[0][1][1]; qnew[58] = spatial_metric_der[0][1][2];
        qnew[59] = spatial_metric_der[0][2][0]; qnew[60] = spatial_metric_der[0][2][1]; qnew[61] = spatial_metric_der[0][2][2];

        qnew[62] = spatial_metric_der[1][0][0]; qnew[63] = spatial_metric_der[1][0][1]; qnew[64] = spatial_metric_der[1][0][2];
        qnew[65] = spatial_metric_der[1][1][0]; qnew[66] = spatial_metric_der[1][1][1]; qnew[67] = spatial_metric_der[1][1][2];
        qnew[68] = spatial_metric_der[1][2][0]; qnew[69] = spatial_metric_der[1][2][1]; qnew[70] = spatial_metric_der[1][2][2];

        qnew[71] = spatial_metric_der[2][0][0]; qnew[72] = spatial_metric_der[2][0][1]; qnew[73] = spatial_metric_der[2][0][2];
        qnew[74] = spatial_metric_der[2][1][0]; qnew[75] = spatial_metric_der[2][1][1]; qnew[76] = spatial_metric_der[2][1][2];
        qnew[77] = spatial_metric_der[2][2][0]; qnew[78] = spatial_metric_der[2][2][1]; qnew[79] = spatial_metric_der[2][2][2];

        if (in_excision_region) {
          for (int i = 0; i < 80; i++) {
            qnew[i] = 0.0;
          }
          qnew[40] = -1.0;
        }

        for (int i = 0; i < 3; i++) {
          gkyl_free(spatial_metric[i]);
          gkyl_free(extrinsic_curvature[i]);
          gkyl_free(shift_der[i]);
      
          for (int j = 0; j < 3; j++) {
            gkyl_free(spatial_metric_der[i][j]);
          }
          gkyl_free(spatial_metric_der[i]);
        }
        gkyl_free(spatial_metric);
        gkyl_free(extrinsic_curvature);
        gkyl_free(shift);
        gkyl_free(lapse_der);
        gkyl_free(shift_der);
        gkyl_free(spatial_metric_der);
        
        qnew[80] = 0.0;
      }
      else {
        qnew[80] += 1.0;
      }
    }
  }
}

void
gr_mhd_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir)
{
  const struct gkyl_wv_eqn* eqn = wv->equation;
  const struct wv_gr_mhd *gr_mhd = container_of(eqn, struct wv_gr_mhd, eqn);
  
  const enum gkyl_spacetime_gauge spacetime_gauge = gr_mhd->spacetime_gauge;

  if (spacetime_gauge == GKYL_STATIC_GAUGE) {
    const struct gkyl_gr_spacetime* spacetime = gr_mhd->spacetime;
    int reinit_freq = gr_mhd->reinit_freq;
    
    for (int i = loidx_c; i<= upidx_c; i++) {
      idxl[dir] = i;

      double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
      double evol_param = qnew[71];

      if (evol_param > reinit_freq) {
        double x = qnew[72];
        double y = qnew[73];
        double z = qnew[74];

        double lapse;
        double *shift = gkyl_malloc(sizeof(double[3]));
        bool in_excision_region;

        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
        }

        double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
        }

        double *lapse_der = gkyl_malloc(sizeof(double[3]));
        double **shift_der = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          shift_der[i] = gkyl_malloc(sizeof(double[3]));
        }

        double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
        for (int i = 0; i < 3; i++) {
          spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
        spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
        spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);

        spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
        spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

        spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
        spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
        spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

        qnew[9] = lapse;
        qnew[10] = shift[0]; qnew[11] = shift[1]; qnew[12] = shift[2];

        qnew[13] = spatial_metric[0][0]; qnew[14] = spatial_metric[0][1]; qnew[15] = spatial_metric[0][2];
        qnew[16] = spatial_metric[1][0]; qnew[17] = spatial_metric[1][1]; qnew[18] = spatial_metric[1][2];
        qnew[19] = spatial_metric[2][0]; qnew[20] = spatial_metric[2][1]; qnew[21] = spatial_metric[2][2];

        qnew[22] = extrinsic_curvature[0][0]; qnew[23] = extrinsic_curvature[0][1]; qnew[24] = extrinsic_curvature[0][2];
        qnew[25] = extrinsic_curvature[1][0]; qnew[26] = extrinsic_curvature[1][1]; qnew[27] = extrinsic_curvature[1][2];
        qnew[28] = extrinsic_curvature[2][0]; qnew[29] = extrinsic_curvature[2][1]; qnew[30] = extrinsic_curvature[2][2];

        if (in_excision_region) {
          qnew[31] = -1.0;
        }
        else {
          qnew[31] = 1.0;
        }

        qnew[32] = lapse_der[0]; qnew[33] = lapse_der[1]; qnew[34] = lapse_der[2];

        qnew[35] = shift_der[0][0]; qnew[36] = shift_der[0][1]; qnew[37] = shift_der[0][2];
        qnew[38] = shift_der[1][0]; qnew[39] = shift_der[1][1]; qnew[40] = shift_der[1][2];
        qnew[41] = shift_der[2][0]; qnew[42] = shift_der[2][1]; qnew[43] = shift_der[2][2];

        qnew[44] = spatial_metric_der[0][0][0]; qnew[45] = spatial_metric_der[0][0][1]; qnew[46] = spatial_metric_der[0][0][2];
        qnew[47] = spatial_metric_der[0][1][0]; qnew[48] = spatial_metric_der[0][1][1]; qnew[49] = spatial_metric_der[0][1][2];
        qnew[50] = spatial_metric_der[0][2][0]; qnew[51] = spatial_metric_der[0][2][1]; qnew[52] = spatial_metric_der[0][2][2];

        qnew[53] = spatial_metric_der[1][0][0]; qnew[54] = spatial_metric_der[1][0][1]; qnew[55] = spatial_metric_der[1][0][2];
        qnew[56] = spatial_metric_der[1][1][0]; qnew[57] = spatial_metric_der[1][1][1]; qnew[58] = spatial_metric_der[1][1][2];
        qnew[59] = spatial_metric_der[1][2][0]; qnew[60] = spatial_metric_der[1][2][1]; qnew[61] = spatial_metric_der[1][2][2];

        qnew[62] = spatial_metric_der[2][0][0]; qnew[63] = spatial_metric_der[2][0][1]; qnew[64] = spatial_metric_der[2][0][2];
        qnew[65] = spatial_metric_der[2][1][0]; qnew[66] = spatial_metric_der[2][1][1]; qnew[67] = spatial_metric_der[2][1][2];
        qnew[68] = spatial_metric_der[2][2][0]; qnew[69] = spatial_metric_der[2][2][1]; qnew[70] = spatial_metric_der[2][2][2];

        if (in_excision_region) {
          for (int i = 0; i < 71; i++) {
            qnew[i] = 0.0;
          }
          qnew[31] = -1.0;
        }

        for (int i = 0; i < 3; i++) {
          gkyl_free(spatial_metric[i]);
          gkyl_free(extrinsic_curvature[i]);
          gkyl_free(shift_der[i]);
      
          for (int j = 0; j < 3; j++) {
            gkyl_free(spatial_metric_der[i][j]);
          }
          gkyl_free(spatial_metric_der[i]);
        }
        gkyl_free(spatial_metric);
        gkyl_free(extrinsic_curvature);
        gkyl_free(shift);
        gkyl_free(lapse_der);
        gkyl_free(shift_der);
        gkyl_free(spatial_metric_der);
        
        qnew[71] = 0.0;
      }
      else {
        qnew[71] += 1.0;
      }
    }
  }
  else if (spacetime_gauge == GKYL_BLACKHOLE_COLLAPSE_GAUGE) {
    const struct gkyl_gr_spacetime* spacetime = gr_mhd->spacetime;
    const struct gr_blackhole *blackhole = container_of(spacetime, struct gr_blackhole, spacetime);

    double mass = blackhole->mass;
    double spin = blackhole->spin;

    double pos_x = blackhole->pos_x;
    double pos_y = blackhole->pos_y;
    double pos_z = blackhole->pos_z;

    for (int i = loidx_c; i <= upidx_c; i++) {
      idxl[dir] = i;

      double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
      double evol_param = qnew[71];

      if (evol_param < mass) {
        evol_param += 0.001 + (0.001 * evol_param);
      }
      qnew[71] = evol_param;

      double x = qnew[72];
      double y = qnew[73];
      double z = qnew[74];

      struct gkyl_gr_spacetime *new_spacetime = gkyl_gr_blackhole_new(false, fmin(evol_param, mass), spin, pos_x, pos_y, pos_z);

      double lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der = gkyl_malloc(sizeof(double[3]));
      double **shift_der = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      new_spacetime->lapse_function_func(new_spacetime, 0.0, x, y, z, &lapse);
      new_spacetime->shift_vector_func(new_spacetime, 0.0, x, y, z, &shift);
      new_spacetime->excision_region_func(new_spacetime, 0.0, x, y, z, &in_excision_region);

      new_spacetime->spatial_metric_tensor_func(new_spacetime, 0.0, x, y, z, &spatial_metric);
      new_spacetime->extrinsic_curvature_tensor_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

      new_spacetime->lapse_function_der_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
      new_spacetime->shift_vector_der_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
      new_spacetime->spatial_metric_tensor_der_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

      qnew[9] = lapse;
      qnew[10] = shift[0]; qnew[11] = shift[1]; qnew[12] = shift[2];

      qnew[13] = spatial_metric[0][0]; qnew[14] = spatial_metric[0][1]; qnew[15] = spatial_metric[0][2];
      qnew[16] = spatial_metric[1][0]; qnew[17] = spatial_metric[1][1]; qnew[18] = spatial_metric[1][2];
      qnew[19] = spatial_metric[2][0]; qnew[20] = spatial_metric[2][1]; qnew[21] = spatial_metric[2][2];

      qnew[22] = extrinsic_curvature[0][0]; qnew[23] = extrinsic_curvature[0][1]; qnew[24] = extrinsic_curvature[0][2];
      qnew[25] = extrinsic_curvature[1][0]; qnew[26] = extrinsic_curvature[1][1]; qnew[27] = extrinsic_curvature[1][2];
      qnew[28] = extrinsic_curvature[2][0]; qnew[29] = extrinsic_curvature[2][1]; qnew[30] = extrinsic_curvature[2][2];

      if (in_excision_region) {
        qnew[31] = -1.0;
      }
      else {
        qnew[31] = 1.0;
      }

      qnew[32] = lapse_der[0]; qnew[33] = lapse_der[1]; qnew[34] = lapse_der[2];

      qnew[35] = shift_der[0][0]; qnew[36] = shift_der[0][1]; qnew[37] = shift_der[0][2];
      qnew[38] = shift_der[1][0]; qnew[39] = shift_der[1][1]; qnew[40] = shift_der[1][2];
      qnew[41] = shift_der[2][0]; qnew[42] = shift_der[2][1]; qnew[43] = shift_der[2][2];

      qnew[44] = spatial_metric_der[0][0][0]; qnew[45] = spatial_metric_der[0][0][1]; qnew[46] = spatial_metric_der[0][0][2];
      qnew[47] = spatial_metric_der[0][1][0]; qnew[48] = spatial_metric_der[0][1][1]; qnew[49] = spatial_metric_der[0][1][2];
      qnew[50] = spatial_metric_der[0][2][0]; qnew[51] = spatial_metric_der[0][2][1]; qnew[52] = spatial_metric_der[0][2][2];

      qnew[53] = spatial_metric_der[1][0][0]; qnew[54] = spatial_metric_der[1][0][1]; qnew[55] = spatial_metric_der[1][0][2];
      qnew[56] = spatial_metric_der[1][1][0]; qnew[57] = spatial_metric_der[1][1][1]; qnew[58] = spatial_metric_der[1][1][2];
      qnew[59] = spatial_metric_der[1][2][0]; qnew[60] = spatial_metric_der[1][2][1]; qnew[61] = spatial_metric_der[1][2][2];

      qnew[62] = spatial_metric_der[2][0][0]; qnew[63] = spatial_metric_der[2][0][1]; qnew[64] = spatial_metric_der[2][0][2];
      qnew[65] = spatial_metric_der[2][1][0]; qnew[66] = spatial_metric_der[2][1][1]; qnew[67] = spatial_metric_der[2][1][2];
      qnew[68] = spatial_metric_der[2][2][0]; qnew[69] = spatial_metric_der[2][2][1]; qnew[70] = spatial_metric_der[2][2][2];

      if (in_excision_region) {
        for (int i = 0; i < 71; i++) {
          qnew[i] = 0.0;
        }
        qnew[31] = -1.0;
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(extrinsic_curvature[i]);
        gkyl_free(shift_der[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der[i][j]);
        }
        gkyl_free(spatial_metric_der[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(extrinsic_curvature);
      gkyl_free(shift);
      gkyl_free(lapse_der);
      gkyl_free(shift_der);
      gkyl_free(spatial_metric_der);
      gkyl_free(new_spacetime);
    }
  }
}

void
gr_mhd_tetrad_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir)
{
  const struct gkyl_wv_eqn* eqn = wv->equation;
  const struct wv_gr_mhd_tetrad *gr_mhd_tetrad = container_of(eqn, struct wv_gr_mhd_tetrad, eqn);
  
  const enum gkyl_spacetime_gauge spacetime_gauge = gr_mhd_tetrad->spacetime_gauge;

  if (spacetime_gauge == GKYL_STATIC_GAUGE) {
    const struct gkyl_gr_spacetime* spacetime = gr_mhd_tetrad->spacetime;
    int reinit_freq = gr_mhd_tetrad->reinit_freq;
    
    for (int i = loidx_c; i<= upidx_c; i++) {
      idxl[dir] = i;

      double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
      double evol_param = qnew[71];

      if (evol_param > reinit_freq) {
        double x = qnew[72];
        double y = qnew[73];
        double z = qnew[74];

        double lapse;
        double *shift = gkyl_malloc(sizeof(double[3]));
        bool in_excision_region;

        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
        }

        double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
        }

        double *lapse_der = gkyl_malloc(sizeof(double[3]));
        double **shift_der = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
          shift_der[i] = gkyl_malloc(sizeof(double[3]));
        }

        double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
        for (int i = 0; i < 3; i++) {
          spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
        spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
        spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);

        spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
        spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

        spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
        spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
        spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

        qnew[9] = lapse;
        qnew[10] = shift[0]; qnew[11] = shift[1]; qnew[12] = shift[2];

        qnew[13] = spatial_metric[0][0]; qnew[14] = spatial_metric[0][1]; qnew[15] = spatial_metric[0][2];
        qnew[16] = spatial_metric[1][0]; qnew[17] = spatial_metric[1][1]; qnew[18] = spatial_metric[1][2];
        qnew[19] = spatial_metric[2][0]; qnew[20] = spatial_metric[2][1]; qnew[21] = spatial_metric[2][2];

        qnew[22] = extrinsic_curvature[0][0]; qnew[23] = extrinsic_curvature[0][1]; qnew[24] = extrinsic_curvature[0][2];
        qnew[25] = extrinsic_curvature[1][0]; qnew[26] = extrinsic_curvature[1][1]; qnew[27] = extrinsic_curvature[1][2];
        qnew[28] = extrinsic_curvature[2][0]; qnew[29] = extrinsic_curvature[2][1]; qnew[30] = extrinsic_curvature[2][2];

        if (in_excision_region) {
          qnew[31] = -1.0;
        }
        else {
          qnew[31] = 1.0;
        }

        qnew[32] = lapse_der[0]; qnew[33] = lapse_der[1]; qnew[34] = lapse_der[2];

        qnew[35] = shift_der[0][0]; qnew[36] = shift_der[0][1]; qnew[37] = shift_der[0][2];
        qnew[38] = shift_der[1][0]; qnew[39] = shift_der[1][1]; qnew[40] = shift_der[1][2];
        qnew[41] = shift_der[2][0]; qnew[42] = shift_der[2][1]; qnew[43] = shift_der[2][2];

        qnew[44] = spatial_metric_der[0][0][0]; qnew[45] = spatial_metric_der[0][0][1]; qnew[46] = spatial_metric_der[0][0][2];
        qnew[47] = spatial_metric_der[0][1][0]; qnew[48] = spatial_metric_der[0][1][1]; qnew[49] = spatial_metric_der[0][1][2];
        qnew[50] = spatial_metric_der[0][2][0]; qnew[51] = spatial_metric_der[0][2][1]; qnew[52] = spatial_metric_der[0][2][2];

        qnew[53] = spatial_metric_der[1][0][0]; qnew[54] = spatial_metric_der[1][0][1]; qnew[55] = spatial_metric_der[1][0][2];
        qnew[56] = spatial_metric_der[1][1][0]; qnew[57] = spatial_metric_der[1][1][1]; qnew[58] = spatial_metric_der[1][1][2];
        qnew[59] = spatial_metric_der[1][2][0]; qnew[60] = spatial_metric_der[1][2][1]; qnew[61] = spatial_metric_der[1][2][2];

        qnew[62] = spatial_metric_der[2][0][0]; qnew[63] = spatial_metric_der[2][0][1]; qnew[64] = spatial_metric_der[2][0][2];
        qnew[65] = spatial_metric_der[2][1][0]; qnew[66] = spatial_metric_der[2][1][1]; qnew[67] = spatial_metric_der[2][1][2];
        qnew[68] = spatial_metric_der[2][2][0]; qnew[69] = spatial_metric_der[2][2][1]; qnew[70] = spatial_metric_der[2][2][2];

        if (in_excision_region) {
          for (int i = 0; i < 71; i++) {
            qnew[i] = 0.0;
          }
          qnew[31] = -1.0;
        }

        for (int i = 0; i < 3; i++) {
          gkyl_free(spatial_metric[i]);
          gkyl_free(extrinsic_curvature[i]);
          gkyl_free(shift_der[i]);
      
          for (int j = 0; j < 3; j++) {
            gkyl_free(spatial_metric_der[i][j]);
          }
          gkyl_free(spatial_metric_der[i]);
        }
        gkyl_free(spatial_metric);
        gkyl_free(extrinsic_curvature);
        gkyl_free(shift);
        gkyl_free(lapse_der);
        gkyl_free(shift_der);
        gkyl_free(spatial_metric_der);
        
        qnew[71] = 0.0;
      }
      else {
        qnew[71] += 1.0;
      }
    }
  }
  else if (spacetime_gauge == GKYL_BLACKHOLE_COLLAPSE_GAUGE) {
    const struct gkyl_gr_spacetime* spacetime = gr_mhd_tetrad->spacetime;
    const struct gr_blackhole *blackhole = container_of(spacetime, struct gr_blackhole, spacetime);

    double mass = blackhole->mass;
    double spin = blackhole->spin;

    double pos_x = blackhole->pos_x;
    double pos_y = blackhole->pos_y;
    double pos_z = blackhole->pos_z;

    for (int i = loidx_c; i <= upidx_c; i++) {
      idxl[dir] = i;

      double *qnew = gkyl_array_fetch(qout, gkyl_range_idx(update_range, idxl));
      double evol_param = qnew[71];

      if (evol_param < mass) {
        evol_param += 0.001 + (0.001 * evol_param);
      }
      qnew[71] = evol_param;

      double x = qnew[72];
      double y = qnew[73];
      double z = qnew[74];

      struct gkyl_gr_spacetime *new_spacetime = gkyl_gr_blackhole_new(false, fmin(evol_param, mass), spin, pos_x, pos_y, pos_z);

      double lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der = gkyl_malloc(sizeof(double[3]));
      double **shift_der = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      new_spacetime->lapse_function_func(new_spacetime, 0.0, x, y, z, &lapse);
      new_spacetime->shift_vector_func(new_spacetime, 0.0, x, y, z, &shift);
      new_spacetime->excision_region_func(new_spacetime, 0.0, x, y, z, &in_excision_region);

      new_spacetime->spatial_metric_tensor_func(new_spacetime, 0.0, x, y, z, &spatial_metric);
      new_spacetime->extrinsic_curvature_tensor_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

      new_spacetime->lapse_function_der_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
      new_spacetime->shift_vector_der_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
      new_spacetime->spatial_metric_tensor_der_func(new_spacetime, 0.0, x, y, z, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

      qnew[9] = lapse;
      qnew[10] = shift[0]; qnew[11] = shift[1]; qnew[12] = shift[2];

      qnew[13] = spatial_metric[0][0]; qnew[14] = spatial_metric[0][1]; qnew[15] = spatial_metric[0][2];
      qnew[16] = spatial_metric[1][0]; qnew[17] = spatial_metric[1][1]; qnew[18] = spatial_metric[1][2];
      qnew[19] = spatial_metric[2][0]; qnew[20] = spatial_metric[2][1]; qnew[21] = spatial_metric[2][2];

      qnew[22] = extrinsic_curvature[0][0]; qnew[23] = extrinsic_curvature[0][1]; qnew[24] = extrinsic_curvature[0][2];
      qnew[25] = extrinsic_curvature[1][0]; qnew[26] = extrinsic_curvature[1][1]; qnew[27] = extrinsic_curvature[1][2];
      qnew[28] = extrinsic_curvature[2][0]; qnew[29] = extrinsic_curvature[2][1]; qnew[30] = extrinsic_curvature[2][2];

      if (in_excision_region) {
        qnew[31] = -1.0;
      }
      else {
        qnew[31] = 1.0;
      }

      qnew[32] = lapse_der[0]; qnew[33] = lapse_der[1]; qnew[34] = lapse_der[2];

      qnew[35] = shift_der[0][0]; qnew[36] = shift_der[0][1]; qnew[37] = shift_der[0][2];
      qnew[38] = shift_der[1][0]; qnew[39] = shift_der[1][1]; qnew[40] = shift_der[1][2];
      qnew[41] = shift_der[2][0]; qnew[42] = shift_der[2][1]; qnew[43] = shift_der[2][2];

      qnew[44] = spatial_metric_der[0][0][0]; qnew[45] = spatial_metric_der[0][0][1]; qnew[46] = spatial_metric_der[0][0][2];
      qnew[47] = spatial_metric_der[0][1][0]; qnew[48] = spatial_metric_der[0][1][1]; qnew[49] = spatial_metric_der[0][1][2];
      qnew[50] = spatial_metric_der[0][2][0]; qnew[51] = spatial_metric_der[0][2][1]; qnew[52] = spatial_metric_der[0][2][2];

      qnew[53] = spatial_metric_der[1][0][0]; qnew[54] = spatial_metric_der[1][0][1]; qnew[55] = spatial_metric_der[1][0][2];
      qnew[56] = spatial_metric_der[1][1][0]; qnew[57] = spatial_metric_der[1][1][1]; qnew[58] = spatial_metric_der[1][1][2];
      qnew[59] = spatial_metric_der[1][2][0]; qnew[60] = spatial_metric_der[1][2][1]; qnew[61] = spatial_metric_der[1][2][2];

      qnew[62] = spatial_metric_der[2][0][0]; qnew[63] = spatial_metric_der[2][0][1]; qnew[64] = spatial_metric_der[2][0][2];
      qnew[65] = spatial_metric_der[2][1][0]; qnew[66] = spatial_metric_der[2][1][1]; qnew[67] = spatial_metric_der[2][1][2];
      qnew[68] = spatial_metric_der[2][2][0]; qnew[69] = spatial_metric_der[2][2][1]; qnew[70] = spatial_metric_der[2][2][2];

      if (in_excision_region) {
        for (int i = 0; i < 71; i++) {
          qnew[i] = 0.0;
        }
        qnew[31] = -1.0;
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(extrinsic_curvature[i]);
        gkyl_free(shift_der[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der[i][j]);
        }
        gkyl_free(spatial_metric_der[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(extrinsic_curvature);
      gkyl_free(shift);
      gkyl_free(lapse_der);
      gkyl_free(shift_der);
      gkyl_free(spatial_metric_der);
      gkyl_free(new_spacetime);
    }
  }
}
