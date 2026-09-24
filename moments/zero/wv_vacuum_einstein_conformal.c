#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_vacuum_einstein_conformal.h>
#include <gkyl_wv_vacuum_einstein_conformal_priv.h>

void
gkyl_vacuum_einstein_conformal_flux(double excision_threshold, enum gkyl_spacetime_slicing spacetime_slicing, enum gkyl_spacetime_evolution spacetime_evolution,
  const double q[77], double flux[77])
{
  double bssn_conformal_fact = q[64];
  double conformal_fact = 1.0 / sqrt(bssn_conformal_fact);

  double bssn_conformal_fact_der[3];
  bssn_conformal_fact_der[0] = q[65]; bssn_conformal_fact_der[1] = q[66]; bssn_conformal_fact_der[2] = q[67];
  double conformal_fact_der[3];
  for (int i = 0; i < 3; i++) {
    conformal_fact_der[i] = -(0.5 * bssn_conformal_fact_der[i]) / (pow(bssn_conformal_fact, 1.5) * conformal_fact);
  }

  double conformal_spatial_metric[3][3];
  conformal_spatial_metric[0][0] = q[0]; conformal_spatial_metric[0][1] = q[1]; conformal_spatial_metric[0][2] = q[2];
  conformal_spatial_metric[1][0] = q[3]; conformal_spatial_metric[1][1] = q[4]; conformal_spatial_metric[1][2] = q[5];
  conformal_spatial_metric[2][0] = q[6]; conformal_spatial_metric[2][1] = q[7]; conformal_spatial_metric[2][2] = q[8];

  double conformal_lapse = q[9];

  double conformal_extrinsic_curvature[3][3];
  conformal_extrinsic_curvature[0][0] = q[10]; conformal_extrinsic_curvature[0][1] = q[11]; conformal_extrinsic_curvature[0][2] = q[12];
  conformal_extrinsic_curvature[1][0] = q[13]; conformal_extrinsic_curvature[1][1] = q[14]; conformal_extrinsic_curvature[1][2] = q[15];
  conformal_extrinsic_curvature[2][0] = q[16]; conformal_extrinsic_curvature[2][1] = q[17]; conformal_extrinsic_curvature[2][2] = q[18];

  double conformal_spatial_metric_der[3][3][3];
  conformal_spatial_metric_der[0][0][0] = q[19]; conformal_spatial_metric_der[0][0][1] = q[20]; conformal_spatial_metric_der[0][0][2] = q[21];
  conformal_spatial_metric_der[0][1][0] = q[22]; conformal_spatial_metric_der[0][1][1] = q[23]; conformal_spatial_metric_der[0][1][2] = q[24];
  conformal_spatial_metric_der[0][2][0] = q[25]; conformal_spatial_metric_der[0][2][1] = q[26]; conformal_spatial_metric_der[0][2][2] = q[27];

  conformal_spatial_metric_der[1][0][0] = q[28]; conformal_spatial_metric_der[1][0][1] = q[29]; conformal_spatial_metric_der[1][0][2] = q[30];
  conformal_spatial_metric_der[1][1][0] = q[31]; conformal_spatial_metric_der[1][1][1] = q[32]; conformal_spatial_metric_der[1][1][2] = q[33];
  conformal_spatial_metric_der[1][2][0] = q[34]; conformal_spatial_metric_der[1][2][1] = q[35]; conformal_spatial_metric_der[1][2][2] = q[36];

  conformal_spatial_metric_der[2][0][0] = q[37]; conformal_spatial_metric_der[2][0][1] = q[38]; conformal_spatial_metric_der[2][0][2] = q[39];
  conformal_spatial_metric_der[2][1][0] = q[40]; conformal_spatial_metric_der[2][1][1] = q[41]; conformal_spatial_metric_der[2][1][2] = q[42];
  conformal_spatial_metric_der[2][2][0] = q[43]; conformal_spatial_metric_der[2][2][1] = q[44]; conformal_spatial_metric_der[2][2][2] = q[45];

  double conformal_lapse_der[3];
  conformal_lapse_der[0] = q[46];
  conformal_lapse_der[1] = q[47];
  conformal_lapse_der[2] = q[48];

  double conformal_aux_vect[3];
  conformal_aux_vect[0] = q[49];
  conformal_aux_vect[1] = q[50];
  conformal_aux_vect[2] = q[51];

  double conformal_shift_vect[3];
  conformal_shift_vect[0] = q[52];
  conformal_shift_vect[1] = q[53];
  conformal_shift_vect[2] = q[54];

  double conformal_shift_vect_der[3][3];
  conformal_shift_vect_der[0][0] = q[55]; conformal_shift_vect_der[0][1] = q[56]; conformal_shift_vect_der[0][2] = q[57];
  conformal_shift_vect_der[1][0] = q[58]; conformal_shift_vect_der[1][1] = q[59]; conformal_shift_vect_der[1][2] = q[60];
  conformal_shift_vect_der[2][0] = q[61]; conformal_shift_vect_der[2][1] = q[62]; conformal_shift_vect_der[2][2] = q[63];

  bool in_excision_region = false;
  if (conformal_lapse < excision_threshold) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double **inv_conformal_spatial_metric = gkyl_malloc(sizeof(double*[3]));
    for (int i = 0; i < 3; i++) {
      inv_conformal_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
    }

    gkyl_vacuum_einstein_conformal_inv_spatial_metric(q, &inv_conformal_spatial_metric);

    double evolution_func = 0.0;
    if (spacetime_evolution == GKYL_RICCI_EVOLUTION) {
      evolution_func = 0.0;
    }
    else if (spacetime_evolution == GKYL_EINSTEIN_EVOLUTION) {
      evolution_func = 1.0;
    }

    double conformal_extrinsic_curvature_trace = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_extrinsic_curvature_trace += inv_conformal_spatial_metric[i][j] * conformal_extrinsic_curvature[i][j];
      }
    }

    double slicing_func = 0.0;
    if (spacetime_slicing == GKYL_GEODESIC_SLICING) {
      slicing_func = 0.0;
    }
    else if (spacetime_slicing == GKYL_HARMONIC_SLICING) {
      slicing_func = conformal_extrinsic_curvature_trace / (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
    }
    else if (spacetime_slicing == GKYL_1PLUSLOG_SLICING) {
      slicing_func = 2.0 * conformal_extrinsic_curvature_trace / (conformal_lapse * (conformal_fact * conformal_fact * conformal_fact * conformal_fact));
    }

    double conformal_spatial_metric_der_raised1[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          conformal_spatial_metric_der_raised1[k][i][j] = 0.0;
          
          for (int l = 0; l < 3; l++) {
            conformal_spatial_metric_der_raised1[k][i][j] += inv_conformal_spatial_metric[k][l] * conformal_spatial_metric_der[l][i][j];
          }
        }
      }
    }

    double conformal_spatial_metric_der_raised3[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          conformal_spatial_metric_der_raised3[i][j][k] = 0.0;

          for (int l = 0; l < 3; l++) {
            conformal_spatial_metric_der_raised3[i][j][k] += inv_conformal_spatial_metric[l][k] * conformal_spatial_metric_der[i][j][l];
          }
        }
      }
    }

    double conformal_aux_vect_raised[3];
    for (int k = 0; k < 3; k++) {
      conformal_aux_vect_raised[k] = 0.0;
        
      for (int l = 0; l < 3; l++) {
        conformal_aux_vect_raised[k] += inv_conformal_spatial_metric[k][l] * conformal_aux_vect[l];
      }
    }

    double conformal_shift_vect_der_lowered[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_shift_vect_der_lowered[i][j] = 0.0;

        for (int k = 0; k < 3; k++) {
          conformal_shift_vect_der_lowered[i][j] += conformal_spatial_metric[k][j] * conformal_shift_vect_der[i][k];
        }
      }
    }

    double conformal_shift_vect_der_switched[3][3];
    for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
        conformal_shift_vect_der_switched[i][k] = 0.0;

        for (int l = 0; l < 3; l++) {
          for (int m = 0; m < 3; m++) {
            conformal_shift_vect_der_switched[i][k] += inv_conformal_spatial_metric[i][l] * conformal_spatial_metric[m][k] * conformal_shift_vect_der[l][m];
          }
        }
      }
    }

    double conformal_symmetrized_shift[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_symmetrized_shift[i][j] = (1.0 / conformal_lapse) * (conformal_shift_vect_der_lowered[i][j] + conformal_shift_vect_der_lowered[j][i]);
      }
    }

    double conformal_extrinsic_curvature_flux[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_extrinsic_curvature_flux[i][j] = -conformal_shift_vect[0] * conformal_extrinsic_curvature[i][j];
        conformal_extrinsic_curvature_flux[i][j] += conformal_lapse * conformal_spatial_metric_der_raised1[0][i][j];
        conformal_extrinsic_curvature_flux[i][j] -= conformal_lapse * (0.5 * evolution_func) * conformal_aux_vect_raised[0] * conformal_spatial_metric[i][j];

        if (i == 0) {
          conformal_extrinsic_curvature_flux[i][j] += 0.5 * conformal_lapse * conformal_lapse_der[j];
          conformal_extrinsic_curvature_flux[i][j] += conformal_lapse * conformal_aux_vect[j];
          for (int r = 0; r < 3; r++) {
            conformal_extrinsic_curvature_flux[i][j] -= 0.5 * conformal_lapse * conformal_spatial_metric_der_raised3[j][r][r];
          }
        }

        if (j == 0) {
          conformal_extrinsic_curvature_flux[i][j] += 0.5 * conformal_lapse * conformal_lapse_der[i];
          conformal_extrinsic_curvature_flux[i][j] += conformal_lapse * conformal_aux_vect[i];
          for (int r = 0; r < 3; r++) {
            conformal_extrinsic_curvature_flux[i][j] -= 0.5 * conformal_lapse * conformal_spatial_metric_der_raised3[i][r][r];
          }
        }
      }
    }

    double conformal_spatial_metric_der_flux[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          conformal_spatial_metric_der_flux[k][i][j] = 0.0;
        }

        for (int r = 0; r < 3; r++) {
          conformal_spatial_metric_der_flux[0][i][j] -= conformal_shift_vect[r] * conformal_spatial_metric_der[r][i][j];
          conformal_spatial_metric_der_flux[0][i][j] -= 2.0 * conformal_shift_vect[r] * conformal_fact_der[r] * conformal_spatial_metric[i][j];
        }

        conformal_spatial_metric_der_flux[0][i][j] += conformal_lapse * (conformal_extrinsic_curvature[i][j] /
          (conformal_fact * conformal_fact * conformal_fact * conformal_fact));
        conformal_spatial_metric_der_flux[0][i][j] -= conformal_lapse * conformal_symmetrized_shift[i][j];
      }
    }

    double conformal_lapse_der_flux[3];
    for (int i = 0; i < 3; i++) {
      conformal_lapse_der_flux[i] = 0.0;
    }
    for (int r = 0; r < 3; r++) {
      conformal_lapse_der_flux[0] -= conformal_shift_vect[r] * conformal_lapse_der[r];
    }
    conformal_lapse_der_flux[0] += conformal_lapse * slicing_func;

    double conformal_aux_vect_flux[3];
    for (int i = 0; i < 3; i++) {
      conformal_aux_vect_flux[i] = -conformal_shift_vect[0] * conformal_aux_vect[i];
      conformal_aux_vect_flux[i] -= 4.0 * conformal_shift_vect[0] * conformal_fact_der[i];
      conformal_aux_vect_flux[i] += conformal_shift_vect_der_switched[0][i];
      conformal_aux_vect_flux[i] -= conformal_shift_vect_der[i][0];
    }
    
    for (int i = 0; i < 10; i++) {
      flux[i] = 0.0;
    }

    flux[10] = conformal_extrinsic_curvature_flux[0][0]; flux[11] = conformal_extrinsic_curvature_flux[0][1]; flux[12] = conformal_extrinsic_curvature_flux[0][2];
    flux[13] = conformal_extrinsic_curvature_flux[1][0]; flux[14] = conformal_extrinsic_curvature_flux[1][1]; flux[15] = conformal_extrinsic_curvature_flux[1][2];
    flux[16] = conformal_extrinsic_curvature_flux[2][0]; flux[17] = conformal_extrinsic_curvature_flux[2][1]; flux[18] = conformal_extrinsic_curvature_flux[2][2];

    flux[19] = conformal_spatial_metric_der_flux[0][0][0]; flux[20] = conformal_spatial_metric_der_flux[0][0][1]; flux[21] = conformal_spatial_metric_der_flux[0][0][2];
    flux[22] = conformal_spatial_metric_der_flux[0][1][0]; flux[23] = conformal_spatial_metric_der_flux[0][1][1]; flux[24] = conformal_spatial_metric_der_flux[0][1][2];
    flux[25] = conformal_spatial_metric_der_flux[0][2][0]; flux[26] = conformal_spatial_metric_der_flux[0][2][1]; flux[27] = conformal_spatial_metric_der_flux[0][2][2];

    flux[28] = conformal_spatial_metric_der_flux[1][0][0]; flux[29] = conformal_spatial_metric_der_flux[1][0][1]; flux[30] = conformal_spatial_metric_der_flux[1][0][2];
    flux[31] = conformal_spatial_metric_der_flux[1][1][0]; flux[32] = conformal_spatial_metric_der_flux[1][1][1]; flux[33] = conformal_spatial_metric_der_flux[1][1][2];
    flux[34] = conformal_spatial_metric_der_flux[1][2][0]; flux[35] = conformal_spatial_metric_der_flux[1][2][1]; flux[36] = conformal_spatial_metric_der_flux[1][2][2];

    flux[37] = conformal_spatial_metric_der_flux[2][0][0]; flux[38] = conformal_spatial_metric_der_flux[2][0][1]; flux[39] = conformal_spatial_metric_der_flux[2][0][2];
    flux[40] = conformal_spatial_metric_der_flux[2][1][0]; flux[41] = conformal_spatial_metric_der_flux[2][1][1]; flux[42] = conformal_spatial_metric_der_flux[2][1][2];
    flux[43] = conformal_spatial_metric_der_flux[2][2][0]; flux[44] = conformal_spatial_metric_der_flux[2][2][1]; flux[45] = conformal_spatial_metric_der_flux[2][2][2];

    flux[46] = conformal_lapse_der_flux[0];
    flux[47] = conformal_lapse_der_flux[1];
    flux[48] = conformal_lapse_der_flux[2];

    flux[49] = conformal_aux_vect_flux[0];
    flux[50] = conformal_aux_vect_flux[1];
    flux[51] = conformal_aux_vect_flux[2];

    for (int i = 52; i < 77; i++) {
      flux[i] = 0.0;
    }
    
    for (int i = 0; i < 3; i++) {
      gkyl_free(inv_conformal_spatial_metric[i]);
    }
    gkyl_free(inv_conformal_spatial_metric);
  }
  else {
    for (int i = 0; i < 77; i++) {
      flux[i] = 0.0;
    }
  }
}

void 
gkyl_vacuum_einstein_conformal_inv_spatial_metric(const double q[77], double ***inv_conformal_spatial_metric)
{
  double conformal_spatial_metric[3][3];
  conformal_spatial_metric[0][0] = q[0]; conformal_spatial_metric[0][1] = q[1]; conformal_spatial_metric[0][2] = q[2];
  conformal_spatial_metric[1][0] = q[3]; conformal_spatial_metric[1][1] = q[4]; conformal_spatial_metric[1][2] = q[5];
  conformal_spatial_metric[2][0] = q[6]; conformal_spatial_metric[2][1] = q[7]; conformal_spatial_metric[2][2] = q[8];

  double conformal_spatial_det =
    (conformal_spatial_metric[0][0] * ((conformal_spatial_metric[1][1] * conformal_spatial_metric[2][2]) - (conformal_spatial_metric[2][1] * conformal_spatial_metric[1][2]))) -
    (conformal_spatial_metric[0][1] * ((conformal_spatial_metric[1][0] * conformal_spatial_metric[2][2]) - (conformal_spatial_metric[1][2] * conformal_spatial_metric[2][0]))) +
    (conformal_spatial_metric[0][2] * ((conformal_spatial_metric[1][0] * conformal_spatial_metric[2][1]) - (conformal_spatial_metric[1][1] * conformal_spatial_metric[2][0])));
  
  double trace = 0.0;
  for (int i = 0; i < 3; i++) {
    trace += conformal_spatial_metric[i][i];
  }

  double conformal_spatial_metric_sq[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      conformal_spatial_metric_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        conformal_spatial_metric_sq[i][j] += conformal_spatial_metric[i][k] * conformal_spatial_metric[k][j];
      }
    }
  }

  double sq_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace += conformal_spatial_metric_sq[i][i];
  }

  double euclidean_metric[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        euclidean_metric[i][j] = 1.0;
      }
      else {
        euclidean_metric[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*inv_conformal_spatial_metric)[i][j] = (1.0 / conformal_spatial_det) *
        ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) - (trace * conformal_spatial_metric[i][j]) + conformal_spatial_metric_sq[i][j]);
    }
  }
}

static inline double
gkyl_vacuum_einstein_conformal_max_abs_speed(double excision_threshold, enum gkyl_spacetime_slicing spacetime_slicing, const double q[77])
{
  double bssn_conformal_fact = q[64];
  double conformal_fact = 1.0 / sqrt(bssn_conformal_fact);

  double conformal_lapse = q[9];

  double conformal_shift_vect[3];
  conformal_shift_vect[0] = q[52];
  conformal_shift_vect[1] = q[53];
  conformal_shift_vect[2] = q[54];

  double **inv_conformal_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_conformal_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_vacuum_einstein_conformal_inv_spatial_metric(q, &inv_conformal_spatial_metric);

  double slicing_func = 0.0;
  if (spacetime_slicing == GKYL_GEODESIC_SLICING) {
    slicing_func = 0.0;
  }
  else if (spacetime_slicing == GKYL_HARMONIC_SLICING) {
    slicing_func = 1.0 / (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
  }
  else if (spacetime_slicing == GKYL_1PLUSLOG_SLICING) {
    slicing_func = 2.0 / (conformal_lapse * (conformal_fact * conformal_fact * conformal_fact * conformal_fact));
  }

  bool in_excision_region = false;
  if (conformal_lapse < excision_threshold) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double max_eig = 0.0;

    for (int i = 0; i < 3; i++) {
      if (fabs(conformal_shift_vect[i]) > max_eig) {
        max_eig = fabs(conformal_shift_vect[i]);
      }
      if (fabs(-conformal_shift_vect[i] + (conformal_lapse * sqrt(slicing_func * inv_conformal_spatial_metric[i][i]))) > max_eig) {
        max_eig = fabs(-conformal_shift_vect[i] + (conformal_lapse * sqrt(slicing_func * inv_conformal_spatial_metric[i][i])));
      }
      if (fabs(-conformal_shift_vect[i] - (conformal_lapse * sqrt(slicing_func * inv_conformal_spatial_metric[i][i]))) > max_eig) {
        max_eig = fabs(-conformal_shift_vect[i] - (conformal_lapse * sqrt(slicing_func * inv_conformal_spatial_metric[i][i])));
      }
    }
    
    for (int i = 0; i < 3; i++) {
      gkyl_free(inv_conformal_spatial_metric[i]);
    }
    gkyl_free(inv_conformal_spatial_metric);

    return max_eig;
  }
  else {
    for (int i = 0; i < 3; i++) {
      gkyl_free(inv_conformal_spatial_metric[i]);
    }
    gkyl_free(inv_conformal_spatial_metric);

    return pow(10.0, -1.0);
  }
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 77; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 77; i++) {
    qout[i] = win[i];
  }
}

static void
vacuum_einstein_conformal_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  // Set spatial metric tensor.
  ghost[0] = 1.0; ghost[1] = 0.0; ghost[2] = 0.0;
  ghost[3] = 0.0; ghost[4] = 1.0; ghost[5] = 0.0;
  ghost[6] = 0.0; ghost[7] = 0.0; ghost[8] = 1.0;

  // Set lapse gauge variable.
  ghost[9] = 1.0;

  // Set extrinsic curvature tensor.
  ghost[10] = 0.0; ghost[11] = 0.0; ghost[12] = 0.0;
  ghost[13] = 0.0; ghost[14] = 0.0; ghost[15] = 0.0;
  ghost[16] = 0.0; ghost[17] = 0.0; ghost[18] = 0.0;

  // Set spatial metric tensor derivatives.
  ghost[19] = 0.0; ghost[20] = 0.0; ghost[21] = 0.0;
  ghost[22] = 0.0; ghost[23] = 0.0; ghost[24] = 0.0;
  ghost[25] = 0.0; ghost[26] = 0.0; ghost[27] = 0.0;

  ghost[28] = 0.0; ghost[29] = 0.0; ghost[30] = 0.0;
  ghost[31] = 0.0; ghost[32] = 0.0; ghost[33] = 0.0;
  ghost[34] = 0.0; ghost[35] = 0.0; ghost[36] = 0.0;

  ghost[37] = 0.0; ghost[38] = 0.0; ghost[39] = 0.0;
  ghost[40] = 0.0; ghost[41] = 0.0; ghost[42] = 0.0;
  ghost[43] = 0.0; ghost[44] = 0.0; ghost[45] = 0.0;

  // Set lapse function derivatives.
  ghost[46] = 0.0; ghost[47] = 0.0; ghost[48] = 0.0;

  // Set auxiliary vector.
  ghost[49] = 0.0; ghost[50] = 0.0; ghost[51] = 0.0;

  // Set shift gauge variables.
  ghost[52] = 0.0; ghost[53] = 0.0; ghost[54] = 0.0;

  // Set shift vector derivatives.
  ghost[55] = 0.0; ghost[56] = 0.0; ghost[57] = 0.0;
  ghost[58] = 0.0; ghost[59] = 0.0; ghost[60] = 0.0;
  ghost[61] = 0.0; ghost[62] = 0.0; ghost[63] = 0.0;

  // Set BSSN conformal factor and derivatives.
  ghost[64] = 1.0;
  ghost[65] = 0.0; ghost[66] = 0.0; ghost[67] = 0.0;
}

static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal,
  double* GKYL_RESTRICT qlocal)
{
  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qglobal[0] * norm[0]) + (qglobal[1] * norm[1]) + (qglobal[2] * norm[2]);
  r1[1] = (qglobal[0] * tau1[0]) + (qglobal[1] * tau1[1]) + (qglobal[2] * tau1[2]);
  r1[2] = (qglobal[0] * tau2[0]) + (qglobal[1] * tau2[1]) + (qglobal[2] * tau2[2]);

  r2[0] = (qglobal[3] * norm[0]) + (qglobal[4] * norm[1]) + (qglobal[5] * norm[2]);
  r2[1] = (qglobal[3] * tau1[0]) + (qglobal[4] * tau1[1]) + (qglobal[5] * tau1[2]);
  r2[2] = (qglobal[3] * tau2[0]) + (qglobal[4] * tau2[1]) + (qglobal[5] * tau2[2]);

  r3[0] = (qglobal[6] * norm[0]) + (qglobal[7] * norm[1]) + (qglobal[8] * norm[2]);
  r3[1] = (qglobal[6] * tau1[0]) + (qglobal[7] * tau1[1]) + (qglobal[8] * tau1[2]);
  r3[2] = (qglobal[6] * tau2[0]) + (qglobal[7] * tau2[1]) + (qglobal[8] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double v1[3], v2[3], v3[3];
  v1[0] = (r1[0] * norm[0]) + (r2[0] * norm[1]) + (r3[0] * norm[2]);
  v1[1] = (r1[0] * tau1[0]) + (r2[0] * tau1[1]) + (r3[0] * tau1[2]);
  v1[2] = (r1[0] * tau2[0]) + (r2[0] * tau2[1]) + (r3[0] * tau2[2]);

  v2[0] = (r1[1] * norm[0]) + (r2[1] * norm[1]) + (r3[1] * norm[2]);
  v2[1] = (r1[1] * tau1[0]) + (r2[1] * tau1[1]) + (r3[1] * tau1[2]);
  v2[2] = (r1[1] * tau2[0]) + (r2[1] * tau2[1]) + (r3[1] * tau2[2]);

  v3[0] = (r1[2] * norm[0]) + (r2[2] * norm[1]) + (r3[2] * norm[2]);
  v3[1] = (r1[2] * tau1[0]) + (r2[2] * tau1[1]) + (r3[2] * tau1[2]);
  v3[2] = (r1[2] * tau2[0]) + (r2[2] * tau2[1]) + (r3[2] * tau2[2]);

  // Rotate spatial metric tensor to local coordinate frame.
  qlocal[0] = v1[0]; qlocal[1] = v1[1]; qlocal[2] = v1[2];
  qlocal[3] = v2[0]; qlocal[4] = v2[1]; qlocal[5] = v2[2];
  qlocal[6] = v3[0]; qlocal[7] = v3[1]; qlocal[8] = v3[2];

  qlocal[9] = qglobal[9];

  // Temporary arrays to store rotated extrinsic column vectors.
  double extr_r1[3], extr_r2[3], extr_r3[3];
  extr_r1[0] = (qglobal[10] * norm[0]) + (qglobal[11] * norm[1]) + (qglobal[12] * norm[2]);
  extr_r1[1] = (qglobal[10] * tau1[0]) + (qglobal[11] * tau1[1]) + (qglobal[12] * tau1[2]);
  extr_r1[2] = (qglobal[10] * tau2[0]) + (qglobal[11] * tau2[1]) + (qglobal[12] * tau2[2]);

  extr_r2[0] = (qglobal[13] * norm[0]) + (qglobal[14] * norm[1]) + (qglobal[15] * norm[2]);
  extr_r2[1] = (qglobal[13] * tau1[0]) + (qglobal[14] * tau1[1]) + (qglobal[15] * tau1[2]);
  extr_r2[2] = (qglobal[13] * tau2[0]) + (qglobal[14] * tau2[1]) + (qglobal[15] * tau2[2]);

  extr_r3[0] = (qglobal[16] * norm[0]) + (qglobal[17] * norm[1]) + (qglobal[18] * norm[2]);
  extr_r3[1] = (qglobal[16] * tau1[0]) + (qglobal[17] * tau1[1]) + (qglobal[18] * tau1[2]);
  extr_r3[2] = (qglobal[16] * tau2[0]) + (qglobal[17] * tau2[1]) + (qglobal[18] * tau2[2]);

  // Temporary arrays to store rotated extrinsic row vectors.
  double inv_v1[3], inv_v2[3], inv_v3[3];
  inv_v1[0] = (extr_r1[0] * norm[0]) + (extr_r2[0] * norm[1]) + (extr_r3[0] * norm[2]);
  inv_v1[1] = (extr_r1[0] * tau1[0]) + (extr_r2[0] * tau1[1]) + (extr_r3[0] * tau1[2]);
  inv_v1[2] = (extr_r1[0] * tau2[0]) + (extr_r2[0] * tau2[1]) + (extr_r3[0] * tau2[2]);

  inv_v2[0] = (extr_r1[1] * norm[0]) + (extr_r2[1] * norm[1]) + (extr_r3[1] * norm[2]);
  inv_v2[1] = (extr_r1[1] * tau1[0]) + (extr_r2[1] * tau1[1]) + (extr_r3[1] * tau1[2]);
  inv_v2[2] = (extr_r1[1] * tau2[0]) + (extr_r2[1] * tau2[1]) + (extr_r3[1] * tau2[2]);

  inv_v3[0] = (extr_r1[2] * norm[0]) + (extr_r2[2] * norm[1]) + (extr_r3[2] * norm[2]);
  inv_v3[1] = (extr_r1[2] * tau1[0]) + (extr_r2[2] * tau1[1]) + (extr_r3[2] * tau1[2]);
  inv_v3[2] = (extr_r1[2] * tau2[0]) + (extr_r2[2] * tau2[1]) + (extr_r3[2] * tau2[2]);

  // Rotate extrinsic curvature tensor to local coordinate frame.
  qlocal[10] = inv_v1[0]; qlocal[11] = inv_v1[1]; qlocal[12] = inv_v1[2];
  qlocal[13] = inv_v2[0]; qlocal[14] = inv_v2[1]; qlocal[15] = inv_v2[2];
  qlocal[16] = inv_v3[0]; qlocal[17] = inv_v3[1]; qlocal[18] = inv_v3[2];

  // Temporary arrays to store rotated column vectors.
  double r11[3], r12[3], r13[3];
  double r21[3], r22[3], r23[3];
  double r31[3], r32[3], r33[3];

  r11[0] = (qglobal[19] * norm[0]) + (qglobal[20] * norm[1]) + (qglobal[21] * norm[2]);
  r11[1] = (qglobal[19] * tau1[0]) + (qglobal[20] * tau1[1]) + (qglobal[21] * tau1[2]);
  r11[2] = (qglobal[19] * tau2[0]) + (qglobal[20] * tau2[1]) + (qglobal[21] * tau2[2]);

  r12[0] = (qglobal[22] * norm[0]) + (qglobal[23] * norm[1]) + (qglobal[24] * norm[2]);
  r12[1] = (qglobal[22] * tau1[0]) + (qglobal[23] * tau1[1]) + (qglobal[24] * tau1[2]);
  r12[2] = (qglobal[22] * tau2[0]) + (qglobal[23] * tau2[1]) + (qglobal[24] * tau2[2]);

  r13[0] = (qglobal[25] * norm[0]) + (qglobal[26] * norm[1]) + (qglobal[27] * norm[2]);
  r13[1] = (qglobal[25] * tau1[0]) + (qglobal[26] * tau1[1]) + (qglobal[27] * tau1[2]);
  r13[2] = (qglobal[25] * tau2[0]) + (qglobal[26] * tau2[1]) + (qglobal[27] * tau2[2]);

  r21[0] = (qglobal[28] * norm[0]) + (qglobal[29] * norm[1]) + (qglobal[30] * norm[2]);
  r21[1] = (qglobal[28] * tau1[0]) + (qglobal[29] * tau1[1]) + (qglobal[30] * tau1[2]);
  r21[2] = (qglobal[28] * tau2[0]) + (qglobal[29] * tau2[1]) + (qglobal[30] * tau2[2]);

  r22[0] = (qglobal[31] * norm[0]) + (qglobal[32] * norm[1]) + (qglobal[33] * norm[2]);
  r22[1] = (qglobal[31] * tau1[0]) + (qglobal[32] * tau1[1]) + (qglobal[33] * tau1[2]);
  r22[2] = (qglobal[31] * tau2[0]) + (qglobal[32] * tau2[1]) + (qglobal[33] * tau2[2]);

  r23[0] = (qglobal[34] * norm[0]) + (qglobal[35] * norm[1]) + (qglobal[36] * norm[2]);
  r23[1] = (qglobal[34] * tau1[0]) + (qglobal[35] * tau1[1]) + (qglobal[36] * tau1[2]);
  r23[2] = (qglobal[34] * tau2[0]) + (qglobal[35] * tau2[1]) + (qglobal[36] * tau2[2]);

  r31[0] = (qglobal[37] * norm[0]) + (qglobal[38] * norm[1]) + (qglobal[39] * norm[2]);
  r31[1] = (qglobal[37] * tau1[0]) + (qglobal[38] * tau1[1]) + (qglobal[39] * tau1[2]);
  r31[2] = (qglobal[37] * tau2[0]) + (qglobal[38] * tau2[1]) + (qglobal[39] * tau2[2]);

  r32[0] = (qglobal[40] * norm[0]) + (qglobal[41] * norm[1]) + (qglobal[42] * norm[2]);
  r32[1] = (qglobal[40] * tau1[0]) + (qglobal[41] * tau1[1]) + (qglobal[42] * tau1[2]);
  r32[2] = (qglobal[40] * tau2[0]) + (qglobal[41] * tau2[1]) + (qglobal[42] * tau2[2]);

  r33[0] = (qglobal[43] * norm[0]) + (qglobal[44] * norm[1]) + (qglobal[45] * norm[2]);
  r33[1] = (qglobal[43] * tau1[0]) + (qglobal[44] * tau1[1]) + (qglobal[45] * tau1[2]);
  r33[2] = (qglobal[43] * tau2[0]) + (qglobal[44] * tau2[1]) + (qglobal[45] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double s11[3], s12[3], s13[3];
  double s21[3], s22[3], s23[3];
  double s31[3], s32[3], s33[3];

  s11[0] = (r11[0] * norm[0]) + (r12[0] * norm[1]) + (r13[0] * norm[2]);
  s11[1] = (r11[1] * norm[0]) + (r12[1] * norm[1]) + (r13[1] * norm[2]);
  s11[2] = (r11[2] * norm[0]) + (r12[2] * norm[1]) + (r13[2] * norm[2]);

  s12[0] = (r11[0] * tau1[0]) + (r12[0] * tau1[1]) + (r13[0] * tau1[2]);
  s12[1] = (r11[1] * tau1[0]) + (r12[1] * tau1[1]) + (r13[1] * tau1[2]);
  s12[2] = (r11[2] * tau1[0]) + (r12[2] * tau1[1]) + (r13[2] * tau1[2]);

  s13[0] = (r11[0] * tau2[0]) + (r12[0] * tau2[1]) + (r13[0] * tau2[2]);
  s13[1] = (r11[1] * tau2[0]) + (r12[1] * tau2[1]) + (r13[1] * tau2[2]);
  s13[2] = (r11[2] * tau2[0]) + (r12[2] * tau2[1]) + (r13[2] * tau2[2]);

  s21[0] = (r21[0] * norm[0]) + (r22[0] * norm[1]) + (r23[0] * norm[2]);
  s21[1] = (r21[1] * norm[0]) + (r22[1] * norm[1]) + (r23[1] * norm[2]);
  s21[2] = (r21[2] * norm[0]) + (r22[2] * norm[1]) + (r23[2] * norm[2]);

  s22[0] = (r21[0] * tau1[0]) + (r22[0] * tau1[1]) + (r23[0] * tau1[2]);
  s22[1] = (r21[1] * tau1[0]) + (r22[1] * tau1[1]) + (r23[1] * tau1[2]);
  s22[2] = (r21[2] * tau1[0]) + (r22[2] * tau1[1]) + (r23[2] * tau1[2]);

  s23[0] = (r21[0] * tau2[0]) + (r22[0] * tau2[1]) + (r23[0] * tau2[2]);
  s23[1] = (r21[1] * tau2[0]) + (r22[1] * tau2[1]) + (r23[1] * tau2[2]);
  s23[2] = (r21[2] * tau2[0]) + (r22[2] * tau2[1]) + (r23[2] * tau2[2]);

  s31[0] = (r31[0] * norm[0]) + (r32[0] * norm[1]) + (r33[0] * norm[2]);
  s31[1] = (r31[1] * norm[0]) + (r32[1] * norm[1]) + (r33[1] * norm[2]);
  s31[2] = (r31[2] * norm[0]) + (r32[2] * norm[1]) + (r33[2] * norm[2]);

  s32[0] = (r31[0] * tau1[0]) + (r32[0] * tau1[1]) + (r33[0] * tau1[2]);
  s32[1] = (r31[1] * tau1[0]) + (r32[1] * tau1[1]) + (r33[1] * tau1[2]);
  s32[2] = (r31[2] * tau1[0]) + (r32[2] * tau1[1]) + (r33[2] * tau1[2]);

  s33[0] = (r31[0] * tau2[0]) + (r32[0] * tau2[1]) + (r33[0] * tau2[2]);
  s33[1] = (r31[1] * tau2[0]) + (r32[1] * tau2[1]) + (r33[1] * tau2[2]);
  s33[2] = (r31[2] * tau2[0]) + (r32[2] * tau2[1]) + (r33[2] * tau2[2]);
  
  // Rotate spatial metric tensor derivative to local coordinate frame.
  qlocal[19] = (s11[0] * norm[0]) + (s21[0] * norm[1]) + (s31[0] * norm[2]);
  qlocal[20] = (s11[1] * norm[0]) + (s21[1] * norm[1]) + (s31[1] * norm[2]);
  qlocal[21] = (s11[2] * norm[0]) + (s21[2] * norm[1]) + (s31[2] * norm[2]);

  qlocal[22] = (s12[0] * norm[0]) + (s22[0] * norm[1]) + (s32[0] * norm[2]);
  qlocal[23] = (s12[1] * norm[0]) + (s22[1] * norm[1]) + (s32[1] * norm[2]);
  qlocal[24] = (s12[2] * norm[0]) + (s22[2] * norm[1]) + (s32[2] * norm[2]);

  qlocal[25] = (s13[0] * norm[0]) + (s23[0] * norm[1]) + (s33[0] * norm[2]);
  qlocal[26] = (s13[1] * norm[0]) + (s23[1] * norm[1]) + (s33[1] * norm[2]);
  qlocal[27] = (s13[2] * norm[0]) + (s23[2] * norm[1]) + (s33[2] * norm[2]);

  qlocal[28] = (s11[0] * tau1[0]) + (s21[0] * tau1[1]) + (s31[0] * tau1[2]);
  qlocal[29] = (s11[1] * tau1[0]) + (s21[1] * tau1[1]) + (s31[1] * tau1[2]);
  qlocal[30] = (s11[2] * tau1[0]) + (s21[2] * tau1[1]) + (s31[2] * tau1[2]);

  qlocal[31] = (s12[0] * tau1[0]) + (s22[0] * tau1[1]) + (s32[0] * tau1[2]);
  qlocal[32] = (s12[1] * tau1[0]) + (s22[1] * tau1[1]) + (s32[1] * tau1[2]);
  qlocal[33] = (s12[2] * tau1[0]) + (s22[2] * tau1[1]) + (s32[2] * tau1[2]);

  qlocal[34] = (s13[0] * tau1[0]) + (s23[0] * tau1[1]) + (s33[0] * tau1[2]);
  qlocal[35] = (s13[1] * tau1[0]) + (s23[1] * tau1[1]) + (s33[1] * tau1[2]);
  qlocal[36] = (s13[2] * tau1[0]) + (s23[2] * tau1[1]) + (s33[2] * tau1[2]);

  qlocal[37] = (s11[0] * tau2[0]) + (s21[0] * tau2[1]) + (s31[0] * tau2[2]);
  qlocal[38] = (s11[1] * tau2[0]) + (s21[1] * tau2[1]) + (s31[1] * tau2[2]);
  qlocal[39] = (s11[2] * tau2[0]) + (s21[2] * tau2[1]) + (s31[2] * tau2[2]);

  qlocal[40] = (s12[0] * tau2[0]) + (s22[0] * tau2[1]) + (s32[0] * tau2[2]);
  qlocal[41] = (s12[1] * tau2[0]) + (s22[1] * tau2[1]) + (s32[1] * tau2[2]);
  qlocal[42] = (s12[2] * tau2[0]) + (s22[2] * tau2[1]) + (s32[2] * tau2[2]);

  qlocal[43] = (s13[0] * tau2[0]) + (s23[0] * tau2[1]) + (s33[0] * tau2[2]);
  qlocal[44] = (s13[1] * tau2[0]) + (s23[1] * tau2[1]) + (s33[1] * tau2[2]);
  qlocal[45] = (s13[2] * tau2[0]) + (s23[2] * tau2[1]) + (s33[2] * tau2[2]);

  qlocal[46] = (qglobal[46] * norm[0]) + (qglobal[47] * norm[1]) + (qglobal[48] * norm[2]);
  qlocal[47] = (qglobal[46] * tau1[0]) + (qglobal[47] * tau1[1]) + (qglobal[48] * tau1[2]);
  qlocal[48] = (qglobal[46] * tau2[0]) + (qglobal[47] * tau2[1]) + (qglobal[48] * tau2[2]);

  qlocal[49] = (qglobal[49] * norm[0]) + (qglobal[50] * norm[1]) + (qglobal[51] * norm[2]);
  qlocal[50] = (qglobal[49] * tau1[0]) + (qglobal[50] * tau1[1]) + (qglobal[51] * tau1[2]);
  qlocal[51] = (qglobal[49] * tau2[0]) + (qglobal[50] * tau2[1]) + (qglobal[51] * tau2[2]);

  qlocal[52] = (qglobal[52] * norm[0]) + (qglobal[53] * norm[1]) + (qglobal[54] * norm[2]);
  qlocal[53] = (qglobal[52] * tau1[0]) + (qglobal[53] * tau1[1]) + (qglobal[54] * tau1[2]);
  qlocal[54] = (qglobal[52] * tau2[0]) + (qglobal[53] * tau2[1]) + (qglobal[54] * tau2[2]);

  // Temporary arrays to store rotated shift derivative column vectors.
  double shiftder_r1[3], shiftder_r2[3], shiftder_r3[3];
  shiftder_r1[0] = (qglobal[55] * norm[0]) + (qglobal[56] * norm[1]) + (qglobal[57] * norm[2]);
  shiftder_r1[1] = (qglobal[55] * tau1[0]) + (qglobal[56] * tau1[1]) + (qglobal[57] * tau1[2]);
  shiftder_r1[2] = (qglobal[55] * tau2[0]) + (qglobal[56] * tau2[1]) + (qglobal[57] * tau2[2]);

  shiftder_r2[0] = (qglobal[58] * norm[0]) + (qglobal[59] * norm[1]) + (qglobal[60] * norm[2]);
  shiftder_r2[1] = (qglobal[58] * tau1[0]) + (qglobal[59] * tau1[1]) + (qglobal[60] * tau1[2]);
  shiftder_r2[2] = (qglobal[58] * tau2[0]) + (qglobal[59] * tau2[1]) + (qglobal[60] * tau2[2]);

  shiftder_r3[0] = (qglobal[61] * norm[0]) + (qglobal[62] * norm[1]) + (qglobal[63] * norm[2]);
  shiftder_r3[1] = (qglobal[61] * tau1[0]) + (qglobal[62] * tau1[1]) + (qglobal[63] * tau1[2]);
  shiftder_r3[2] = (qglobal[61] * tau2[0]) + (qglobal[62] * tau2[1]) + (qglobal[63] * tau2[2]);

  // Temporary arrays to store rotated shift derivative row vectors.
  double shiftder_v1[3], shiftder_v2[3], shiftder_v3[3];
  shiftder_v1[0] = (shiftder_r1[0] * norm[0]) + (shiftder_r2[0] * norm[1]) + (shiftder_r3[0] * norm[2]);
  shiftder_v1[1] = (shiftder_r1[0] * tau1[0]) + (shiftder_r2[0] * tau1[1]) + (shiftder_r3[0] * tau1[2]);
  shiftder_v1[2] = (shiftder_r1[0] * tau2[0]) + (shiftder_r2[0] * tau2[1]) + (shiftder_r3[0] * tau2[2]);

  shiftder_v2[0] = (shiftder_r1[1] * norm[0]) + (shiftder_r2[1] * norm[1]) + (shiftder_r3[1] * norm[2]);
  shiftder_v2[1] = (shiftder_r1[1] * tau1[0]) + (shiftder_r2[1] * tau1[1]) + (shiftder_r3[1] * tau1[2]);
  shiftder_v2[2] = (shiftder_r1[1] * tau2[0]) + (shiftder_r2[1] * tau2[1]) + (shiftder_r3[1] * tau2[2]);

  shiftder_v3[0] = (shiftder_r1[2] * norm[0]) + (shiftder_r2[2] * norm[1]) + (shiftder_r3[2] * norm[2]);
  shiftder_v3[1] = (shiftder_r1[2] * tau1[0]) + (shiftder_r2[2] * tau1[1]) + (shiftder_r3[2] * tau1[2]);
  shiftder_v3[2] = (shiftder_r1[2] * tau2[0]) + (shiftder_r2[2] * tau2[1]) + (shiftder_r3[2] * tau2[2]);

  // Rotate shift vector derivative to local coordinate frame.
  qlocal[55] = shiftder_v1[0]; qlocal[56] = shiftder_v1[1]; qlocal[57] = shiftder_v1[2];
  qlocal[58] = shiftder_v2[0]; qlocal[59] = shiftder_v2[1]; qlocal[60] = shiftder_v2[2];
  qlocal[61] = shiftder_v3[0]; qlocal[62] = shiftder_v3[1]; qlocal[63] = shiftder_v3[2];

  qlocal[64] = qglobal[64];
  qlocal[65] = (qglobal[65] * norm[0]) + (qglobal[66] * norm[1]) + (qglobal[67] * norm[2]);
  qlocal[66] = (qglobal[65] * tau1[0]) + (qglobal[66] * tau1[1]) + (qglobal[67] * tau1[2]);
  qlocal[67] = (qglobal[65] * tau2[0]) + (qglobal[66] * tau2[1]) + (qglobal[67] * tau2[2]);
}

static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal,
  double* GKYL_RESTRICT qglobal)
{
  // Temporary arrays to store rotated column vectors.
  double r1[3], r2[3], r3[3];
  r1[0] = (qlocal[0] * norm[0]) + (qlocal[1] * tau1[0]) + (qlocal[2] * tau2[0]);
  r1[1] = (qlocal[0] * norm[1]) + (qlocal[1] * tau1[1]) + (qlocal[2] * tau2[1]);
  r1[2] = (qlocal[0] * norm[2]) + (qlocal[1] * tau1[2]) + (qlocal[2] * tau2[2]);

  r2[0] = (qlocal[3] * norm[0]) + (qlocal[4] * tau1[0]) + (qlocal[5] * tau2[0]);
  r2[1] = (qlocal[3] * norm[1]) + (qlocal[4] * tau1[1]) + (qlocal[5] * tau2[1]);
  r2[2] = (qlocal[3] * norm[2]) + (qlocal[4] * tau1[2]) + (qlocal[5] * tau2[2]);

  r3[0] = (qlocal[6] * norm[0]) + (qlocal[7] * tau1[0]) + (qlocal[8] * tau2[0]);
  r3[1] = (qlocal[6] * norm[1]) + (qlocal[7] * tau1[1]) + (qlocal[8] * tau2[1]);
  r3[2] = (qlocal[6] * norm[2]) + (qlocal[7] * tau1[2]) + (qlocal[8] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double v1[3], v2[3], v3[3];
  v1[0] = (r1[0] * norm[0]) + (r2[0] * tau1[0]) + (r3[0] * tau2[0]);
  v1[1] = (r1[0] * norm[1]) + (r2[0] * tau1[1]) + (r3[0] * tau2[1]);
  v1[2] = (r1[0] * norm[2]) + (r2[0] * tau1[2]) + (r3[0] * tau2[2]);

  v2[0] = (r1[1] * norm[0]) + (r2[1] * tau1[0]) + (r3[1] * tau2[0]);
  v2[1] = (r1[1] * norm[1]) + (r2[1] * tau1[1]) + (r3[1] * tau2[1]);
  v2[2] = (r1[1] * norm[2]) + (r2[1] * tau1[2]) + (r3[1] * tau2[2]);

  v3[0] = (r1[2] * norm[0]) + (r2[2] * tau1[0]) + (r3[2] * tau2[0]);
  v3[1] = (r1[2] * norm[1]) + (r2[2] * tau1[1]) + (r3[2] * tau2[1]);
  v3[2] = (r1[2] * norm[2]) + (r2[2] * tau1[2]) + (r3[2] * tau2[2]);

  // Rotate spatial metric tensor back to global coordinate frame.
  qglobal[0] = v1[0]; qglobal[1] = v1[1]; qglobal[2] = v1[2];
  qglobal[3] = v2[0]; qglobal[4] = v2[1]; qglobal[5] = v2[2];
  qglobal[6] = v3[0]; qglobal[7] = v3[1]; qglobal[8] = v3[2];

  qglobal[9] = qlocal[9];

  // Temporary arrays to store rotated extrinsic column vectors.
  double extr_r1[3], extr_r2[3], extr_r3[3];
  extr_r1[0] = (qlocal[10] * norm[0]) + (qlocal[11] * tau1[0]) + (qlocal[12] * tau2[0]);
  extr_r1[1] = (qlocal[10] * norm[1]) + (qlocal[11] * tau1[1]) + (qlocal[12] * tau2[1]);
  extr_r1[2] = (qlocal[10] * norm[2]) + (qlocal[11] * tau1[2]) + (qlocal[12] * tau2[2]);

  extr_r2[0] = (qlocal[13] * norm[0]) + (qlocal[14] * tau1[0]) + (qlocal[15] * tau2[0]);
  extr_r2[1] = (qlocal[13] * norm[1]) + (qlocal[14] * tau1[1]) + (qlocal[15] * tau2[1]);
  extr_r2[2] = (qlocal[13] * norm[2]) + (qlocal[14] * tau1[2]) + (qlocal[15] * tau2[2]);

  extr_r3[0] = (qlocal[16] * norm[0]) + (qlocal[17] * tau1[0]) + (qlocal[18] * tau2[0]);
  extr_r3[1] = (qlocal[16] * norm[1]) + (qlocal[17] * tau1[1]) + (qlocal[18] * tau2[1]);
  extr_r3[2] = (qlocal[16] * norm[2]) + (qlocal[17] * tau1[2]) + (qlocal[18] * tau2[2]);

  // Temporary arrays to store rotated extrinsic row vectors.
  double inv_v1[3], inv_v2[3], inv_v3[3];
  inv_v1[0] = (extr_r1[0] * norm[0]) + (extr_r2[0] * tau1[0]) + (extr_r3[0] * tau2[0]);
  inv_v1[1] = (extr_r1[0] * norm[1]) + (extr_r2[0] * tau1[1]) + (extr_r3[0] * tau2[1]);
  inv_v1[2] = (extr_r1[0] * norm[2]) + (extr_r2[0] * tau1[2]) + (extr_r3[0] * tau2[2]);

  inv_v2[0] = (extr_r1[1] * norm[0]) + (extr_r2[1] * tau1[0]) + (extr_r3[1] * tau2[0]);
  inv_v2[1] = (extr_r1[1] * norm[1]) + (extr_r2[1] * tau1[1]) + (extr_r3[1] * tau2[1]);
  inv_v2[2] = (extr_r1[1] * norm[2]) + (extr_r2[1] * tau1[2]) + (extr_r3[1] * tau2[2]);

  inv_v3[0] = (extr_r1[2] * norm[0]) + (extr_r2[2] * tau1[0]) + (extr_r3[2] * tau2[0]);
  inv_v3[1] = (extr_r1[2] * norm[1]) + (extr_r2[2] * tau1[1]) + (extr_r3[2] * tau2[1]);
  inv_v3[2] = (extr_r1[2] * norm[2]) + (extr_r2[2] * tau1[2]) + (extr_r3[2] * tau2[2]);

  // Rotate extrinsic curvature tensor back to global coordinate frame.
  qglobal[10] = inv_v1[0]; qglobal[11] = inv_v1[1]; qglobal[12] = inv_v1[2];
  qglobal[13] = inv_v2[0]; qglobal[14] = inv_v2[1]; qglobal[15] = inv_v2[2];
  qglobal[16] = inv_v3[0]; qglobal[17] = inv_v3[1]; qglobal[18] = inv_v3[2];

  // Temporary arrays to store rotated column vectors.
  double r11[3], r12[3], r13[3];
  double r21[3], r22[3], r23[3];
  double r31[3], r32[3], r33[3];

  r11[0] = (qlocal[19] * norm[0]) + (qlocal[20] * tau1[0]) + (qlocal[21] * tau2[0]);
  r11[1] = (qlocal[19] * norm[1]) + (qlocal[20] * tau1[1]) + (qlocal[21] * tau2[1]);
  r11[2] = (qlocal[19] * norm[2]) + (qlocal[20] * tau1[2]) + (qlocal[21] * tau2[2]);

  r12[0] = (qlocal[22] * norm[0]) + (qlocal[23] * tau1[0]) + (qlocal[24] * tau2[0]);
  r12[1] = (qlocal[22] * norm[1]) + (qlocal[23] * tau1[1]) + (qlocal[24] * tau2[1]);
  r12[2] = (qlocal[22] * norm[2]) + (qlocal[23] * tau1[2]) + (qlocal[24] * tau2[2]);

  r13[0] = (qlocal[25] * norm[0]) + (qlocal[26] * tau1[0]) + (qlocal[27] * tau2[0]);
  r13[1] = (qlocal[25] * norm[1]) + (qlocal[26] * tau1[1]) + (qlocal[27] * tau2[1]);
  r13[2] = (qlocal[25] * norm[2]) + (qlocal[26] * tau1[2]) + (qlocal[27] * tau2[2]);

  r21[0] = (qlocal[28] * norm[0]) + (qlocal[29] * tau1[0]) + (qlocal[30] * tau2[0]);
  r21[1] = (qlocal[28] * norm[1]) + (qlocal[29] * tau1[1]) + (qlocal[30] * tau2[1]);
  r21[2] = (qlocal[28] * norm[2]) + (qlocal[29] * tau1[2]) + (qlocal[30] * tau2[2]);

  r22[0] = (qlocal[31] * norm[0]) + (qlocal[32] * tau1[0]) + (qlocal[33] * tau2[0]);
  r22[1] = (qlocal[31] * norm[1]) + (qlocal[32] * tau1[1]) + (qlocal[33] * tau2[1]);
  r22[2] = (qlocal[31] * norm[2]) + (qlocal[32] * tau1[2]) + (qlocal[33] * tau2[2]);

  r23[0] = (qlocal[34] * norm[0]) + (qlocal[35] * tau1[0]) + (qlocal[36] * tau2[0]);
  r23[1] = (qlocal[34] * norm[1]) + (qlocal[35] * tau1[1]) + (qlocal[36] * tau2[1]);
  r23[2] = (qlocal[34] * norm[2]) + (qlocal[35] * tau1[2]) + (qlocal[36] * tau2[2]);

  r31[0] = (qlocal[37] * norm[0]) + (qlocal[38] * tau1[0]) + (qlocal[39] * tau2[0]);
  r31[1] = (qlocal[37] * norm[1]) + (qlocal[38] * tau1[1]) + (qlocal[39] * tau2[1]);
  r31[2] = (qlocal[37] * norm[2]) + (qlocal[38] * tau1[2]) + (qlocal[39] * tau2[2]);

  r32[0] = (qlocal[40] * norm[0]) + (qlocal[41] * tau1[0]) + (qlocal[42] * tau2[0]);
  r32[1] = (qlocal[40] * norm[1]) + (qlocal[41] * tau1[1]) + (qlocal[42] * tau2[1]);
  r32[2] = (qlocal[40] * norm[2]) + (qlocal[41] * tau1[2]) + (qlocal[42] * tau2[2]);

  r33[0] = (qlocal[43] * norm[0]) + (qlocal[44] * tau1[0]) + (qlocal[45] * tau2[0]);
  r33[1] = (qlocal[43] * norm[1]) + (qlocal[44] * tau1[1]) + (qlocal[45] * tau2[1]);
  r33[2] = (qlocal[43] * norm[2]) + (qlocal[44] * tau1[2]) + (qlocal[45] * tau2[2]);

  // Temporary arrays to store rotated row vectors.
  double s11[3], s12[3], s13[3];
  double s21[3], s22[3], s23[3];
  double s31[3], s32[3], s33[3];

  s11[0] = (r11[0] * norm[0]) + (r12[0] * tau1[0]) + (r13[0] * tau2[0]);
  s11[1] = (r11[1] * norm[0]) + (r12[1] * tau1[0]) + (r13[1] * tau2[0]);
  s11[2] = (r11[2] * norm[0]) + (r12[2] * tau1[0]) + (r13[2] * tau2[0]);

  s12[0] = (r11[0] * norm[1]) + (r12[0] * tau1[1]) + (r13[0] * tau2[1]);
  s12[1] = (r11[1] * norm[1]) + (r12[1] * tau1[1]) + (r13[1] * tau2[1]);
  s12[2] = (r11[2] * norm[1]) + (r12[2] * tau1[1]) + (r13[2] * tau2[1]);

  s13[0] = (r11[0] * norm[2]) + (r12[0] * tau1[2]) + (r13[0] * tau2[2]);
  s13[1] = (r11[1] * norm[2]) + (r12[1] * tau1[2]) + (r13[1] * tau2[2]);
  s13[2] = (r11[2] * norm[2]) + (r12[2] * tau1[2]) + (r13[2] * tau2[2]);

  s21[0] = (r21[0] * norm[0]) + (r22[0] * tau1[0]) + (r23[0] * tau2[0]);
  s21[1] = (r21[1] * norm[0]) + (r22[1] * tau1[0]) + (r23[1] * tau2[0]);
  s21[2] = (r21[2] * norm[0]) + (r22[2] * tau1[0]) + (r23[2] * tau2[0]);

  s22[0] = (r21[0] * norm[1]) + (r22[0] * tau1[1]) + (r23[0] * tau2[1]);
  s22[1] = (r21[1] * norm[1]) + (r22[1] * tau1[1]) + (r23[1] * tau2[1]);
  s22[2] = (r21[2] * norm[1]) + (r22[2] * tau1[1]) + (r23[2] * tau2[1]);

  s23[0] = (r21[0] * norm[2]) + (r22[0] * tau1[2]) + (r23[0] * tau2[2]);
  s23[1] = (r21[1] * norm[2]) + (r22[1] * tau1[2]) + (r23[1] * tau2[2]);
  s23[2] = (r21[2] * norm[2]) + (r22[2] * tau1[2]) + (r23[2] * tau2[2]);

  s31[0] = (r31[0] * norm[0]) + (r32[0] * tau1[0]) + (r33[0] * tau2[0]);
  s31[1] = (r31[1] * norm[0]) + (r32[1] * tau1[0]) + (r33[1] * tau2[0]);
  s31[2] = (r31[2] * norm[0]) + (r32[2] * tau1[0]) + (r33[2] * tau2[0]);

  s32[0] = (r31[0] * norm[1]) + (r32[0] * tau1[1]) + (r33[0] * tau2[1]);
  s32[1] = (r31[1] * norm[1]) + (r32[1] * tau1[1]) + (r33[1] * tau2[1]);
  s32[2] = (r31[2] * norm[1]) + (r32[2] * tau1[1]) + (r33[2] * tau2[1]);

  s33[0] = (r31[0] * norm[2]) + (r32[0] * tau1[2]) + (r33[0] * tau2[2]);
  s33[1] = (r31[1] * norm[2]) + (r32[1] * tau1[2]) + (r33[1] * tau2[2]);
  s33[2] = (r31[2] * norm[2]) + (r32[2] * tau1[2]) + (r33[2] * tau2[2]);

  // Rotate spatial metric tensor derivative back to global coordinate frame.
  qglobal[19] = (s11[0] * norm[0]) + (s21[0] * tau1[0]) + (s31[0] * tau2[0]);
  qglobal[20] = (s11[1] * norm[0]) + (s21[1] * tau1[0]) + (s31[1] * tau2[0]);
  qglobal[21] = (s11[2] * norm[0]) + (s21[2] * tau1[0]) + (s31[2] * tau2[0]);

  qglobal[22] = (s12[0] * norm[0]) + (s22[0] * tau1[0]) + (s32[0] * tau2[0]);
  qglobal[23] = (s12[1] * norm[0]) + (s22[1] * tau1[0]) + (s32[1] * tau2[0]);
  qglobal[24] = (s12[2] * norm[0]) + (s22[2] * tau1[0]) + (s32[2] * tau2[0]);

  qglobal[25] = (s13[0] * norm[0]) + (s23[0] * tau1[0]) + (s33[0] * tau2[0]);
  qglobal[26] = (s13[1] * norm[0]) + (s23[1] * tau1[0]) + (s33[1] * tau2[0]);
  qglobal[27] = (s13[2] * norm[0]) + (s23[2] * tau1[0]) + (s33[2] * tau2[0]);

  qglobal[28] = (s11[0] * norm[1]) + (s21[0] * tau1[1]) + (s31[0] * tau2[1]);
  qglobal[29] = (s11[1] * norm[1]) + (s21[1] * tau1[1]) + (s31[1] * tau2[1]);
  qglobal[30] = (s11[2] * norm[1]) + (s21[2] * tau1[1]) + (s31[2] * tau2[1]);

  qglobal[31] = (s12[0] * norm[1]) + (s22[0] * tau1[1]) + (s32[0] * tau2[1]);
  qglobal[32] = (s12[1] * norm[1]) + (s22[1] * tau1[1]) + (s32[1] * tau2[1]);
  qglobal[33] = (s12[2] * norm[1]) + (s22[2] * tau1[1]) + (s32[2] * tau2[1]);

  qglobal[34] = (s13[0] * norm[1]) + (s23[0] * tau1[1]) + (s33[0] * tau2[1]);
  qglobal[35] = (s13[1] * norm[1]) + (s23[1] * tau1[1]) + (s33[1] * tau2[1]);
  qglobal[36] = (s13[2] * norm[1]) + (s23[2] * tau1[1]) + (s33[2] * tau2[1]);

  qglobal[37] = (s11[0] * norm[2]) + (s21[0] * tau1[2]) + (s31[0] * tau2[2]);
  qglobal[38] = (s11[1] * norm[2]) + (s21[1] * tau1[2]) + (s31[1] * tau2[2]);
  qglobal[39] = (s11[2] * norm[2]) + (s21[2] * tau1[2]) + (s31[2] * tau2[2]);

  qglobal[40] = (s12[0] * norm[2]) + (s22[0] * tau1[2]) + (s32[0] * tau2[2]);
  qglobal[41] = (s12[1] * norm[2]) + (s22[1] * tau1[2]) + (s32[1] * tau2[2]);
  qglobal[42] = (s12[2] * norm[2]) + (s22[2] * tau1[2]) + (s32[2] * tau2[2]);

  qglobal[43] = (s13[0] * norm[2]) + (s23[0] * tau1[2]) + (s33[0] * tau2[2]);
  qglobal[44] = (s13[1] * norm[2]) + (s23[1] * tau1[2]) + (s33[1] * tau2[2]);
  qglobal[45] = (s13[2] * norm[2]) + (s23[2] * tau1[2]) + (s33[2] * tau2[2]);

  qglobal[46] = (qlocal[46] * norm[0]) + (qlocal[47] * tau1[0]) + (qlocal[48] * tau2[0]);
  qglobal[47] = (qlocal[46] * norm[1]) + (qlocal[47] * tau1[1]) + (qlocal[48] * tau2[1]);
  qglobal[48] = (qlocal[46] * norm[2]) + (qlocal[47] * tau1[2]) + (qlocal[48] * tau2[2]);

  qglobal[49] = (qlocal[49] * norm[0]) + (qlocal[50] * tau1[0]) + (qlocal[51] * tau2[0]);
  qglobal[50] = (qlocal[49] * norm[1]) + (qlocal[50] * tau1[1]) + (qlocal[51] * tau2[1]);
  qglobal[51] = (qlocal[49] * norm[2]) + (qlocal[50] * tau1[2]) + (qlocal[51] * tau2[2]);

  qglobal[52] = (qlocal[52] * norm[0]) + (qlocal[53] * tau1[0]) + (qlocal[54] * tau2[0]);
  qglobal[53] = (qlocal[52] * norm[1]) + (qlocal[53] * tau1[1]) + (qlocal[54] * tau2[1]);
  qglobal[54] = (qlocal[52] * norm[2]) + (qlocal[53] * tau1[2]) + (qlocal[54] * tau2[2]);

  // Temporary arrays to store rotated shift derivative column vectors.
  double shiftder_r1[3], shiftder_r2[3], shiftder_r3[3];
  shiftder_r1[0] = (qlocal[55] * norm[0]) + (qlocal[56] * tau1[0]) + (qlocal[57] * tau2[0]);
  shiftder_r1[1] = (qlocal[55] * norm[1]) + (qlocal[56] * tau1[1]) + (qlocal[57] * tau2[1]);
  shiftder_r1[2] = (qlocal[55] * norm[2]) + (qlocal[56] * tau1[2]) + (qlocal[57] * tau2[2]);

  shiftder_r2[0] = (qlocal[58] * norm[0]) + (qlocal[59] * tau1[0]) + (qlocal[60] * tau2[0]);
  shiftder_r2[1] = (qlocal[58] * norm[1]) + (qlocal[59] * tau1[1]) + (qlocal[60] * tau2[1]);
  shiftder_r2[2] = (qlocal[58] * norm[2]) + (qlocal[59] * tau1[2]) + (qlocal[60] * tau2[2]);

  shiftder_r3[0] = (qlocal[61] * norm[0]) + (qlocal[62] * tau1[0]) + (qlocal[63] * tau2[0]);
  shiftder_r3[1] = (qlocal[61] * norm[1]) + (qlocal[62] * tau1[1]) + (qlocal[63] * tau2[1]);
  shiftder_r3[2] = (qlocal[61] * norm[2]) + (qlocal[62] * tau1[2]) + (qlocal[63] * tau2[2]);

  // Temporary arrays to store rotated shift derivative row vectors.
  double shiftder_v1[3], shiftder_v2[3], shiftder_v3[3];
  shiftder_v1[0] = (shiftder_r1[0] * norm[0]) + (shiftder_r2[0] * tau1[0]) + (shiftder_r3[0] * tau2[0]);
  shiftder_v1[1] = (shiftder_r1[0] * norm[1]) + (shiftder_r2[0] * tau1[1]) + (shiftder_r3[0] * tau2[1]);
  shiftder_v1[2] = (shiftder_r1[0] * norm[2]) + (shiftder_r2[0] * tau1[2]) + (shiftder_r3[0] * tau2[2]);

  shiftder_v2[0] = (shiftder_r1[1] * norm[0]) + (shiftder_r2[1] * tau1[0]) + (shiftder_r3[1] * tau2[0]);
  shiftder_v2[1] = (shiftder_r1[1] * norm[1]) + (shiftder_r2[1] * tau1[1]) + (shiftder_r3[1] * tau2[1]);
  shiftder_v2[2] = (shiftder_r1[1] * norm[2]) + (shiftder_r2[1] * tau1[2]) + (shiftder_r3[1] * tau2[2]);

  shiftder_v3[0] = (shiftder_r1[2] * norm[0]) + (shiftder_r2[2] * tau1[0]) + (shiftder_r3[2] * tau2[0]);
  shiftder_v3[1] = (shiftder_r1[2] * norm[1]) + (shiftder_r2[2] * tau1[1]) + (shiftder_r3[2] * tau2[1]);
  shiftder_v3[2] = (shiftder_r1[2] * norm[2]) + (shiftder_r2[2] * tau1[2]) + (shiftder_r3[2] * tau2[2]);

  // Rotate shift vector derivative back to global coordinate frame.
  qglobal[55] = shiftder_v1[0]; qglobal[56] = shiftder_v1[1]; qglobal[57] = shiftder_v1[2];
  qglobal[58] = shiftder_v2[0]; qglobal[59] = shiftder_v2[1]; qglobal[60] = shiftder_v2[2];
  qglobal[61] = shiftder_v3[0]; qglobal[62] = shiftder_v3[1]; qglobal[63] = shiftder_v3[2];

  qglobal[64] = qlocal[64];
  qglobal[65] = (qlocal[65] * norm[0]) + (qlocal[66] * tau1[0]) + (qlocal[67] * tau2[0]);
  qglobal[66] = (qlocal[65] * norm[1]) + (qlocal[66] * tau1[1]) + (qlocal[67] * tau2[1]);
  qglobal[67] = (qlocal[65] * norm[2]) + (qlocal[66] * tau1[2]) + (qlocal[67] * tau2[2]);
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_vacuum_einstein_conformal *vacuum_einstein_conformal = container_of(eqn, struct wv_vacuum_einstein_conformal, eqn);
  double excision_threshold = vacuum_einstein_conformal->excision_threshold;
  enum gkyl_spacetime_slicing spacetime_slicing = vacuum_einstein_conformal->spacetime_slicing;
  enum gkyl_spacetime_evolution spacetime_evolution = vacuum_einstein_conformal->spacetime_evolution;

  double sl = gkyl_vacuum_einstein_conformal_max_abs_speed(excision_threshold, spacetime_slicing, ql);
  double sr = gkyl_vacuum_einstein_conformal_max_abs_speed(excision_threshold, spacetime_slicing, qr);
  double amax = fmax(sl, sr);

  double fl[77], fr[77];
  gkyl_vacuum_einstein_conformal_flux(excision_threshold, spacetime_slicing, spacetime_evolution, ql, fl);
  gkyl_vacuum_einstein_conformal_flux(excision_threshold, spacetime_slicing, spacetime_evolution, qr, fr);

  bool in_excision_region_l = false;
  if (ql[9] < excision_threshold) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[9] < excision_threshold) {
    in_excision_region_r = true;
  }

  double *w0 = &waves[0], *w1 = &waves[77];
  if (!in_excision_region_l && !in_excision_region_r) {
    for (int i = 0; i < 77; i++) {
      w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr[i] - fl[i]) / amax);
      w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr[i] - fl[i]) / amax);
    }
  }
  else {
    for (int i = 0; i < 77; i++) {
      w0[i] = 0.0;
      w1[i] = 0.0;
    }
  }

  s[0] = -amax;
  s[1] = amax;

  return s[1];
}

static void
qfluct_lax(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[77];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 77; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr,  const double phil, const double phir,
  double* waves, double* s)
{
  return wave_lax(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double phil, const double phir,
  const double* waves, const double* s, double* amdq, double* apdq)
{
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static double
wave_hll(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_vacuum_einstein_conformal *vacuum_einstein_conformal = container_of(eqn, struct wv_vacuum_einstein_conformal, eqn);
  double excision_threshold = vacuum_einstein_conformal->excision_threshold;
  enum gkyl_spacetime_slicing spacetime_slicing = vacuum_einstein_conformal->spacetime_slicing;
  enum gkyl_spacetime_evolution spacetime_evolution = vacuum_einstein_conformal->spacetime_evolution;

  double conformal_lapse_l = ql[9];
  double conformal_shiftx_l = ql[52];

  double bssn_conformal_fact_l = ql[64];
  double conformal_fact_l = 1.0 / sqrt(bssn_conformal_fact_l);

  double slicing_func_l = 0.0;
  if (spacetime_slicing == GKYL_GEODESIC_SLICING) {
    slicing_func_l = 0.0;
  }
  else if (spacetime_slicing == GKYL_HARMONIC_SLICING) {
    slicing_func_l = 1.0 / (conformal_fact_l * conformal_fact_l * conformal_fact_l * conformal_fact_l);
  }
  else if (spacetime_slicing == GKYL_1PLUSLOG_SLICING) {
    slicing_func_l = 2.0 / (conformal_lapse_l * (conformal_fact_l * conformal_fact_l * conformal_fact_l * conformal_fact_l));
  }

  double **inv_conformal_spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_conformal_spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_vacuum_einstein_conformal_inv_spatial_metric(ql, &inv_conformal_spatial_metric_l);

  bool in_excision_region_l = false;
  if (conformal_lapse_l < excision_threshold) {
    in_excision_region_l = true;
  }

  double conformal_lapse_r = qr[9];
  double conformal_shiftx_r = qr[52];

  double bssn_conformal_fact_r = qr[64];
  double conformal_fact_r = 1.0 / sqrt(bssn_conformal_fact_r);

  double slicing_func_r = 0.0;
  if (spacetime_slicing == GKYL_GEODESIC_SLICING) {
    slicing_func_r = 0.0;
  }
  else if (spacetime_slicing == GKYL_HARMONIC_SLICING) {
    slicing_func_r = 1.0 / (conformal_fact_r * conformal_fact_r * conformal_fact_r * conformal_fact_r);
  }
  else if (spacetime_slicing == GKYL_1PLUSLOG_SLICING) {
    slicing_func_r = 2.0 / (conformal_lapse_l / (conformal_fact_r * conformal_fact_r * conformal_fact_r * conformal_fact_r));
  }

  double **inv_conformal_spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_conformal_spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_vacuum_einstein_conformal_inv_spatial_metric(qr, &inv_conformal_spatial_metric_r);

  bool in_excision_region_r = false;
  if (conformal_lapse_r < excision_threshold) {
    in_excision_region_r = true;
  }

  double vx_l = -conformal_shiftx_l;
  double vx_r = -conformal_shiftx_r;
  double c_sl = conformal_lapse_l * sqrt(slicing_func_l * inv_conformal_spatial_metric_l[0][0]);
  double c_sr = conformal_lapse_r * sqrt(slicing_func_r * inv_conformal_spatial_metric_r[0][0]);

  double vx_avg = 0.5 * (vx_l + vx_r);
  double cs_avg = 0.5 * (c_sl + c_sr);

  double sl = (vx_avg - cs_avg) / (1.0 - (vx_avg * cs_avg));
  double sr = (vx_avg + cs_avg) / (1.0 + (vx_avg * cs_avg));

  double fl[77], fr[77];
  gkyl_vacuum_einstein_conformal_flux(excision_threshold, spacetime_slicing, spacetime_evolution, ql, fl);
  gkyl_vacuum_einstein_conformal_flux(excision_threshold, spacetime_slicing, spacetime_evolution, qr, fr);

  double qm[77];
  for (int i = 0; i < 77; i++) {
    qm[i] = ((sr * qr[i]) - (sl * ql[i]) + (fl[i] - fr[i])) / (sr - sl);
  }

  double *w0 = &waves[0], *w1 = &waves[77];
  if (!in_excision_region_l && !in_excision_region_r) {
    for (int i = 0; i < 77; i++) {
      w0[i] = qm[i] - ql[i];
      w1[i] = qr[i] - qm[i];
    }
  }
  else {
    for (int i = 0; i < 77; i++) {
      w0[i] = 0.0;
      w1[i] = 0.0;
    }
  }

  s[0] = sl;
  s[1] = sr;

  for (int i = 0; i < 3; i++) {
    gkyl_free(inv_conformal_spatial_metric_l[i]);
    gkyl_free(inv_conformal_spatial_metric_r[i]);
  }
  gkyl_free(inv_conformal_spatial_metric_l);
  gkyl_free(inv_conformal_spatial_metric_r);

  return fmax(fabs(sl), fabs(sr));
}

static void
qfluct_hll(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[77];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 77; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_hll_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, const double phil, const double phir,
  double* waves, double* s)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return wave_hll(eqn, delta, ql, qr, waves, s);
  }
  else {
    return wave_lax(eqn, delta, ql, qr, waves, s);
  }

  return 0.0; // Unreachable code.
}

static void
qfluct_hll_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double phil, const double phir,
  const double* waves, const double* s, double* amdq, double* apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return qfluct_hll(eqn, ql, qr, waves, s, amdq, apdq);
  }
  else {
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
  }
}

static double
flux_jump(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, double* flux_jump)
{
  const struct wv_vacuum_einstein_conformal *vacuum_einstein_conformal = container_of(eqn, struct wv_vacuum_einstein_conformal, eqn);
  double excision_threshold = vacuum_einstein_conformal->excision_threshold;
  enum gkyl_spacetime_slicing spacetime_slicing = vacuum_einstein_conformal->spacetime_slicing;
  enum gkyl_spacetime_evolution spacetime_evolution = vacuum_einstein_conformal->spacetime_evolution;

  double fr[77], fl[77];
  gkyl_vacuum_einstein_conformal_flux(excision_threshold, spacetime_slicing, spacetime_evolution, ql, fl);
  gkyl_vacuum_einstein_conformal_flux(excision_threshold, spacetime_slicing, spacetime_evolution, qr, fr);

  bool in_excision_region_l = false;
  if (ql[9] < excision_threshold) {
    in_excision_region_l = true;
  }

  bool in_excision_region_r = false;
  if (qr[9] < excision_threshold) {
    in_excision_region_r = true;
  }

  if (!in_excision_region_l && !in_excision_region_r) {
    for (int m = 0; m < 77; m++) {
      flux_jump[m] = fr[m] - fl[m];
    }
  }
  else {
    for (int m = 0; m < 77; m++) {
      flux_jump[m] = 0.0;
    }
  }

  double amaxl = gkyl_vacuum_einstein_conformal_max_abs_speed(excision_threshold, spacetime_slicing, ql);
  double amaxr = gkyl_vacuum_einstein_conformal_max_abs_speed(excision_threshold, spacetime_slicing, qr);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  if (q[9] < 0.0) {
    return false;
  }
  else {
    return true;
  }
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_vacuum_einstein_conformal *vacuum_einstein_conformal = container_of(eqn, struct wv_vacuum_einstein_conformal, eqn);
  double excision_threshold = vacuum_einstein_conformal->excision_threshold;
  enum gkyl_spacetime_slicing spacetime_slicing = vacuum_einstein_conformal->spacetime_slicing;

  return gkyl_vacuum_einstein_conformal_max_abs_speed(excision_threshold, spacetime_slicing, q);
}

static inline void
vacuum_einstein_conformal_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  diag[0] = qin[9];
}

static inline void
vacuum_einstein_conformal_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  const struct wv_vacuum_einstein_conformal *vacuum_einstein_conformal = container_of(eqn, struct wv_vacuum_einstein_conformal, eqn);
  double excision_threshold = vacuum_einstein_conformal->excision_threshold;
  enum gkyl_spacetime_slicing spacetime_slicing = vacuum_einstein_conformal->spacetime_slicing;
  enum gkyl_spacetime_evolution spacetime_evolution = vacuum_einstein_conformal->spacetime_evolution;

  double bssn_conformal_fact = qin[64];
  double conformal_fact = 1.0 / sqrt(bssn_conformal_fact);

  double bssn_conformal_fact_der[3];
  bssn_conformal_fact_der[0] = qin[65]; bssn_conformal_fact_der[1] = qin[66]; bssn_conformal_fact_der[2] = qin[67];
  double conformal_fact_der[3];
  for (int i = 0; i < 3; i++) {
    conformal_fact_der[i] = -(0.5 * bssn_conformal_fact_der[i]) / (pow(bssn_conformal_fact, 1.5) * conformal_fact);
  }

  double bssn_conformal_fact_der2[3][3];
  bssn_conformal_fact_der2[0][0] = qin[68]; bssn_conformal_fact_der2[0][1] = qin[69]; bssn_conformal_fact_der2[0][2] = qin[70];
  bssn_conformal_fact_der2[1][0] = qin[71]; bssn_conformal_fact_der2[1][1] = qin[72]; bssn_conformal_fact_der2[1][2] = qin[73];
  bssn_conformal_fact_der2[2][0] = qin[74]; bssn_conformal_fact_der2[2][1] = qin[75]; bssn_conformal_fact_der2[2][2] = qin[76];
  double conformal_fact_der2[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      conformal_fact_der2[i][j] = -(0.5 * bssn_conformal_fact_der2[i][j]) / (pow(bssn_conformal_fact, 1.5) * conformal_fact);
      conformal_fact_der2[i][j] += (0.75 * bssn_conformal_fact_der[i] * bssn_conformal_fact_der[j]) / (pow(bssn_conformal_fact, 2.5) * conformal_fact);
    }
  }

  double conformal_spatial_metric[3][3];
  conformal_spatial_metric[0][0] = qin[0]; conformal_spatial_metric[0][1] = qin[1]; conformal_spatial_metric[0][2] = qin[2];
  conformal_spatial_metric[1][0] = qin[3]; conformal_spatial_metric[1][1] = qin[4]; conformal_spatial_metric[1][2] = qin[5];
  conformal_spatial_metric[2][0] = qin[6]; conformal_spatial_metric[2][1] = qin[7]; conformal_spatial_metric[2][2] = qin[8];

  double conformal_lapse = qin[9];

  double conformal_extrinsic_curvature[3][3];
  conformal_extrinsic_curvature[0][0] = qin[10]; conformal_extrinsic_curvature[0][1] = qin[11]; conformal_extrinsic_curvature[0][2] = qin[12];
  conformal_extrinsic_curvature[1][0] = qin[13]; conformal_extrinsic_curvature[1][1] = qin[14]; conformal_extrinsic_curvature[1][2] = qin[15];
  conformal_extrinsic_curvature[2][0] = qin[16]; conformal_extrinsic_curvature[2][1] = qin[17]; conformal_extrinsic_curvature[2][2] = qin[18];

  double conformal_spatial_metric_der[3][3][3];
  conformal_spatial_metric_der[0][0][0] = qin[19]; conformal_spatial_metric_der[0][0][1] = qin[20]; conformal_spatial_metric_der[0][0][2] = qin[21];
  conformal_spatial_metric_der[0][1][0] = qin[22]; conformal_spatial_metric_der[0][1][1] = qin[23]; conformal_spatial_metric_der[0][1][2] = qin[24];
  conformal_spatial_metric_der[0][2][0] = qin[25]; conformal_spatial_metric_der[0][2][1] = qin[26]; conformal_spatial_metric_der[0][2][2] = qin[27];

  conformal_spatial_metric_der[1][0][0] = qin[28]; conformal_spatial_metric_der[1][0][1] = qin[29]; conformal_spatial_metric_der[1][0][2] = qin[30];
  conformal_spatial_metric_der[1][1][0] = qin[31]; conformal_spatial_metric_der[1][1][1] = qin[32]; conformal_spatial_metric_der[1][1][2] = qin[33];
  conformal_spatial_metric_der[1][2][0] = qin[34]; conformal_spatial_metric_der[1][2][1] = qin[35]; conformal_spatial_metric_der[1][2][2] = qin[36];

  conformal_spatial_metric_der[2][0][0] = qin[37]; conformal_spatial_metric_der[2][0][1] = qin[38]; conformal_spatial_metric_der[2][0][2] = qin[39];
  conformal_spatial_metric_der[2][1][0] = qin[40]; conformal_spatial_metric_der[2][1][1] = qin[41]; conformal_spatial_metric_der[2][1][2] = qin[42];
  conformal_spatial_metric_der[2][2][0] = qin[43]; conformal_spatial_metric_der[2][2][1] = qin[44]; conformal_spatial_metric_der[2][2][2] = qin[45];

  double conformal_lapse_der[3];
  conformal_lapse_der[0] = qin[46];
  conformal_lapse_der[1] = qin[47];
  conformal_lapse_der[2] = qin[48];

  double conformal_aux_vect[3];
  conformal_aux_vect[0] = qin[49];
  conformal_aux_vect[1] = qin[50];
  conformal_aux_vect[2] = qin[51];

  double conformal_shift_vect[3];
  conformal_shift_vect[0] = qin[52];
  conformal_shift_vect[1] = qin[53];
  conformal_shift_vect[2] = qin[54];

  double conformal_shift_vect_der[3][3];
  conformal_shift_vect_der[0][0] = qin[55]; conformal_shift_vect_der[0][1] = qin[56]; conformal_shift_vect_der[0][2] = qin[57];
  conformal_shift_vect_der[1][0] = qin[58]; conformal_shift_vect_der[1][1] = qin[59]; conformal_shift_vect_der[1][2] = qin[60];
  conformal_shift_vect_der[2][0] = qin[61]; conformal_shift_vect_der[2][1] = qin[62]; conformal_shift_vect_der[2][2] = qin[63];

  bool in_excision_region = false;
  if (conformal_lapse < excision_threshold) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double **inv_conformal_spatial_metric = gkyl_malloc(sizeof(double*[3]));
    for (int i = 0; i < 3; i++) {
      inv_conformal_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
    }

    gkyl_vacuum_einstein_conformal_inv_spatial_metric(qin, &inv_conformal_spatial_metric);

    double evolution_func = 0.0;
    if (spacetime_evolution == GKYL_RICCI_EVOLUTION) {
      evolution_func = 0.0;
    }
    else if (spacetime_evolution == GKYL_EINSTEIN_EVOLUTION) {
      evolution_func = 1.0;
    }

    double conformal_extrinsic_curvature_trace = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_extrinsic_curvature_trace += inv_conformal_spatial_metric[i][j] * conformal_extrinsic_curvature[i][j];
      }
    }

    double conformal_extrinsic_curvature_mixed[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_extrinsic_curvature_mixed[i][j] = 0.0;

        for (int l = 0; l < 3; l++) {
          conformal_extrinsic_curvature_mixed[i][j] += inv_conformal_spatial_metric[l][j] * conformal_extrinsic_curvature[i][l];
        }
      }
    }

    double conformal_extrinsic_curvature_raised[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_extrinsic_curvature_raised[i][j] = 0.0;

        for (int l = 0; l < 3; l++) {
          for (int m = 0; m < 3; m++) {
            conformal_extrinsic_curvature_raised[i][j] += inv_conformal_spatial_metric[i][l] * inv_conformal_spatial_metric[m][j] * conformal_extrinsic_curvature[l][m];
          }
        }
      }
    }

    double conformal_shift_vect_der_lowered[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_shift_vect_der_lowered[i][j] = 0.0;

        for (int k = 0; k < 3; k++) {
          conformal_shift_vect_der_lowered[i][j] += conformal_spatial_metric[k][j] * conformal_shift_vect_der[i][k];
        }
      }
    }

    double conformal_shift_vect_der_switched[3][3];
    for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
        conformal_shift_vect_der_switched[i][k] = 0.0;

        for (int l = 0; l < 3; l++) {
          for (int m = 0; m < 3; m++) {
            conformal_shift_vect_der_switched[i][k] += inv_conformal_spatial_metric[i][l] * conformal_spatial_metric[m][k] * conformal_shift_vect_der[l][m];
          }
        }
      }
    }

    double conformal_shift_vect_der_trace = 0.0;
    for (int i = 0; i < 3; i++) {
      conformal_shift_vect_der_trace += conformal_shift_vect_der[i][i];
    }

    double conformal_symmetrized_shift[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_symmetrized_shift[i][j] = (1.0 / conformal_lapse) * (conformal_shift_vect_der_lowered[i][j] + conformal_shift_vect_der_lowered[j][i]);
      }
    }

    double slicing_func = 0.0;
    if (spacetime_slicing == GKYL_GEODESIC_SLICING) {
      slicing_func = 0.0;
    }
    else if (spacetime_slicing == GKYL_HARMONIC_SLICING) {
      slicing_func = conformal_extrinsic_curvature_trace / (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
    }
    else if (spacetime_slicing == GKYL_1PLUSLOG_SLICING) {
      slicing_func = 2.0 * conformal_extrinsic_curvature_trace / (conformal_lapse * (conformal_fact * conformal_fact * conformal_fact * conformal_fact));
    }

    double conformal_spatial_metric_der_raised1[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          conformal_spatial_metric_der_raised1[k][i][j] = 0.0;
          
          for (int l = 0; l < 3; l++) {
            conformal_spatial_metric_der_raised1[k][i][j] += inv_conformal_spatial_metric[k][l] * conformal_spatial_metric_der[l][i][j];
          }
        }
      }
    }

    double conformal_spatial_metric_der_raised3[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          conformal_spatial_metric_der_raised3[i][j][k] = 0.0;

          for (int l = 0; l < 3; l++) {
            conformal_spatial_metric_der_raised3[i][j][k] += inv_conformal_spatial_metric[l][k] * conformal_spatial_metric_der[i][j][l];
          }
        }
      }
    }

    double conformal_spatial_metric_der_lowered1[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          conformal_spatial_metric_der_lowered1[i][j][k] = 0.0;

          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              conformal_spatial_metric_der_lowered1[i][j][k] += inv_conformal_spatial_metric[j][l] * inv_conformal_spatial_metric[m][k] * conformal_spatial_metric_der[i][l][m];
            }
          }
        }
      }
    }

    double conformal_spatial_metric_der_lowered3[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          conformal_spatial_metric_der_lowered3[i][j][k] = 0.0;

          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              conformal_spatial_metric_der_lowered3[i][j][k] += inv_conformal_spatial_metric[i][l] * inv_conformal_spatial_metric[m][j] * conformal_spatial_metric_der[l][m][k];
            }
          }
        }
      }
    }

    double conformal_aux_vect_raised[3];
    for (int k = 0; k < 3; k++) {
      conformal_aux_vect_raised[k] = 0.0;
        
      for (int l = 0; l < 3; l++) {
        conformal_aux_vect_raised[k] += inv_conformal_spatial_metric[k][l] * conformal_aux_vect[l];
      }
    }

    double conformal_spatial_christoffel[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          conformal_spatial_christoffel[i][k][l] = 0.0;

          for (int m = 0; m < 3; m++) {
            conformal_spatial_christoffel[i][k][l] += inv_conformal_spatial_metric[i][m] * conformal_spatial_metric_der[l][m][k];
            conformal_spatial_christoffel[i][k][l] += inv_conformal_spatial_metric[i][m] * conformal_spatial_metric_der[k][m][l];
            conformal_spatial_christoffel[i][k][l] -= inv_conformal_spatial_metric[i][m] * conformal_spatial_metric_der[m][k][l];
          }
        }
      }
    }

    double conformal_lapse_der_raised[3];
    for (int k = 0; k < 3; k++) {
      conformal_lapse_der_raised[k] = 0.0;

      for (int l = 0; l < 3; l++) {
        conformal_lapse_der_raised[k] += inv_conformal_spatial_metric[k][l] * conformal_lapse_der[l];
      }
    }

    double conformal_fact_der_raised[3];
    for (int k = 0; k < 3; k++) {
      conformal_fact_der_raised[k] = 0.0;

      for (int l = 0; l < 3; l++) {
        conformal_fact_der_raised[k] += inv_conformal_spatial_metric[k][l] * conformal_fact_der[l];
      }
    }

    double Y_tensor[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        Y_tensor[i][j] = 2.0 * conformal_fact_der2[i][j];

        for (int k = 0; k < 3; k++) {
          Y_tensor[i][j] += 2.0 * conformal_spatial_metric[i][j] * conformal_fact_der_raised[k] * conformal_fact_der[k];
        }

        for (int r = 0; r < 3; r++) {
          Y_tensor[i][j] -= 2.0 * conformal_fact_der[r] * conformal_spatial_christoffel[r][i][j];
        }
        Y_tensor[i][j] -= 6.0 * conformal_fact_der[i] * conformal_fact_der[j];
      }
    }

    double Y_tensor_mixed[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        Y_tensor_mixed[i][j] = 0.0;

        for (int l = 0; l < 3; l++) {
          Y_tensor_mixed[i][j] += inv_conformal_spatial_metric[l][j] * Y_tensor[i][l];
        }
      }
    }

    double conformal_spatial_metric_source[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_spatial_metric_source[i][j] = -2.0 * conformal_lapse * (conformal_extrinsic_curvature[i][j] / (conformal_fact * conformal_fact * conformal_fact * conformal_fact));
        conformal_spatial_metric_source[i][j] += 2.0 * conformal_lapse * conformal_symmetrized_shift[i][j];

        for (int r = 0; r < 3; r++) {
          conformal_spatial_metric_source[i][j] += 2.0 * conformal_shift_vect[r] * conformal_spatial_metric_der[r][i][j];
          conformal_spatial_metric_source[i][j] += 4.0 * conformal_shift_vect[r] * conformal_fact_der[r] * conformal_spatial_metric[i][j];
        }
      }
    }

    double conformal_lapse_source = -(conformal_lapse * conformal_lapse) * slicing_func;
    for (int r = 0; r < 3; r++) {
      conformal_lapse_source += conformal_lapse * conformal_shift_vect[r] * conformal_lapse_der[r];
    }

    double conformal_extrinsic_curvature_source[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_extrinsic_curvature_source[i][j] = 0.0;

        for (int r = 0; r < 3; r++) {
          conformal_extrinsic_curvature_source[i][j] += 2.0 * conformal_extrinsic_curvature[i][r] * conformal_shift_vect_der[j][r];
          conformal_extrinsic_curvature_source[i][j] += 2.0 * conformal_extrinsic_curvature[j][r] * conformal_shift_vect_der[i][r];
          conformal_extrinsic_curvature_source[i][j] -= 2.0 * conformal_extrinsic_curvature[i][j] * conformal_shift_vect_der[r][r];
        }

        for (int k = 0; k < 3; k++) {
          conformal_extrinsic_curvature_source[i][j] -= (2.0 * conformal_lapse * conformal_extrinsic_curvature_mixed[i][k] * conformal_extrinsic_curvature[k][j]) /
            (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
          conformal_extrinsic_curvature_source[i][j] += (conformal_lapse * conformal_extrinsic_curvature_trace * conformal_extrinsic_curvature[i][j]) /
            (conformal_fact * conformal_fact * conformal_fact * conformal_fact);

          for (int r = 0; r < 3; r++) {
            conformal_extrinsic_curvature_source[i][j] -= conformal_lapse * conformal_spatial_christoffel[k][r][i] * conformal_spatial_christoffel[r][k][j];

            conformal_extrinsic_curvature_source[i][j] += 2.0 * conformal_lapse * conformal_spatial_metric_der_raised3[i][k][r] * conformal_spatial_metric_der_raised3[r][j][k];
            conformal_extrinsic_curvature_source[i][j] += 2.0 * conformal_lapse * conformal_spatial_metric_der_raised3[j][k][r] * conformal_spatial_metric_der_raised3[r][i][k];
            conformal_extrinsic_curvature_source[i][j] += conformal_lapse * conformal_spatial_christoffel[k][k][r] * conformal_spatial_christoffel[r][i][j];

            conformal_extrinsic_curvature_source[i][j] -= conformal_lapse * (2.0 * conformal_spatial_metric_der_raised3[k][r][k] - conformal_lapse_der[r]) *
              (conformal_spatial_metric_der_raised3[i][j][r] + conformal_spatial_metric_der_raised3[j][i][r]);
          }

          conformal_extrinsic_curvature_source[i][j] += conformal_lapse * conformal_lapse_der[i] * (conformal_aux_vect[j] - (0.5 * conformal_spatial_metric_der_raised3[j][k][k]));
          conformal_extrinsic_curvature_source[i][j] += conformal_lapse * conformal_lapse_der[j] * (conformal_aux_vect[i] - (0.5 * conformal_spatial_metric_der_raised3[i][k][k]));

          conformal_extrinsic_curvature_source[i][j] -= conformal_lapse * evolution_func * conformal_aux_vect_raised[k] * conformal_spatial_metric_der[k][i][j];
        }

        for (int k = 0; k < 3; k++) {
          for (int r = 0; r < 3; r++) {
            for (int s = 0; s < 3; s++) {
              conformal_extrinsic_curvature_source[i][j] -= (0.25 * evolution_func * conformal_lapse * conformal_spatial_metric[i][j]) *
                conformal_spatial_metric_der_lowered1[k][r][s] * conformal_spatial_christoffel[k][r][s];
              conformal_extrinsic_curvature_source[i][j] += (0.25 * evolution_func * conformal_lapse * conformal_spatial_metric[i][j]) *
                conformal_spatial_metric_der_raised3[k][r][r] * conformal_spatial_metric_der_lowered3[k][s][s];
            }
          }

          conformal_extrinsic_curvature_source[i][j] -= (0.5 * evolution_func * conformal_lapse * conformal_spatial_metric[i][j]) * conformal_aux_vect_raised[k] * conformal_lapse_der[k];
        }

        for (int r = 0; r < 3; r++) {
          for (int s = 0; s < 3; s++) {
            conformal_extrinsic_curvature_source[i][j] += ((0.25 * evolution_func * conformal_lapse * conformal_spatial_metric[i][j]) *
              conformal_extrinsic_curvature_raised[r][s] * conformal_extrinsic_curvature[r][s]) / (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
          }
        }

        conformal_extrinsic_curvature_source[i][j] -= ((0.25 * evolution_func * conformal_lapse * conformal_spatial_metric[i][j]) *
          (conformal_extrinsic_curvature_trace * conformal_extrinsic_curvature_trace)) / (conformal_fact * conformal_fact * conformal_fact * conformal_fact);

        conformal_extrinsic_curvature_source[i][j] -= Y_tensor[i][j];
        conformal_extrinsic_curvature_source[i][j] += 2.0 * conformal_lapse_der[i] * conformal_fact_der[j];
        conformal_extrinsic_curvature_source[i][j] += 2.0 * conformal_lapse_der[j] * conformal_fact_der[i];
        
        for (int k = 0; k < 3; k++) {
          conformal_extrinsic_curvature_source[i][j] += conformal_spatial_metric[i][j] * (evolution_func - 1.0) * Y_tensor_mixed[k][k];
          conformal_extrinsic_curvature_source[i][j] -= 2.0 * conformal_spatial_metric[i][j] * conformal_lapse_der_raised[k] * conformal_fact_der[k];
        }
      }
    }

    double conformal_aux_vect_source[3];
    for (int i = 0; i < 3; i++) {
      conformal_aux_vect_source[i] = 0.0;

      for (int r = 0; r < 3; r++) {
        conformal_aux_vect_source[i] += (conformal_lapse * conformal_lapse_der[r] * conformal_extrinsic_curvature_mixed[i][r]) /
          (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
        
        if (i == r) {
          conformal_aux_vect_source[i] -= (conformal_lapse * conformal_lapse_der[r] * conformal_extrinsic_curvature_trace) /
            (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
        }

        for (int s = 0; s < 3; s++) {
          conformal_aux_vect_source[i] += (conformal_lapse * conformal_extrinsic_curvature_mixed[s][r] * conformal_spatial_metric_der_raised3[i][r][s]) /
            (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
          conformal_aux_vect_source[i] -= (2.0 * conformal_lapse * conformal_extrinsic_curvature_mixed[s][r] * conformal_spatial_metric_der_raised3[r][i][s]) /
            (conformal_fact * conformal_fact * conformal_fact * conformal_fact);

          conformal_aux_vect_source[i] -= (conformal_lapse * conformal_extrinsic_curvature_mixed[i][r] * conformal_spatial_metric_der_raised3[r][s][s]) /
            (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
          conformal_aux_vect_source[i] += (2.0 * conformal_lapse * conformal_extrinsic_curvature_mixed[i][r] * conformal_spatial_metric_der_raised3[s][r][s]) /
            (conformal_fact * conformal_fact * conformal_fact * conformal_fact);

          conformal_aux_vect_source[i] -= (6.0 * conformal_lapse * conformal_fact_der[r] * conformal_extrinsic_curvature_mixed[s][r]) /
            (conformal_fact * conformal_fact * conformal_fact * conformal_fact);

          if (s == r) {
            conformal_aux_vect_source[i] += (2.0 * conformal_lapse * conformal_fact_der[r] * conformal_extrinsic_curvature_trace) /
              (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
          }
        }

        conformal_aux_vect_source[i] += 2.0 * conformal_shift_vect_der[i][r] * conformal_aux_vect[r];

        if (i == r) {
          conformal_aux_vect_source[i] -= 2.0 * conformal_shift_vect_der_trace * conformal_aux_vect[r];
        }

        for (int s = 0; s < 3; s++) {
          conformal_aux_vect_source[i] += 2.0 * conformal_spatial_metric_der_raised3[r][i][s] * conformal_shift_vect_der_switched[r][s];

          if (i == s) {
            for (int j = 0; j < 3; j++) {
              conformal_aux_vect_source[i] -= 2.0 * conformal_spatial_metric_der_raised1[j][j][r] * conformal_shift_vect_der_switched[r][s];
            }
          }
        }

        conformal_aux_vect_source[i] += 4.0 * conformal_shift_vect_der[i][r] * conformal_fact_der[r];
        conformal_aux_vect_source[i] -= 4.0 * conformal_shift_vect_der_trace * conformal_fact_der[i];
      }
    }

    double bssn_conformal_fact_source = (1.0 / 3.0) * conformal_lapse * conformal_extrinsic_curvature_trace * bssn_conformal_fact;
    bssn_conformal_fact_source -= (8.0 / 3.0) * conformal_shift_vect_der_trace * bssn_conformal_fact;

    sout[0] = conformal_spatial_metric_source[0][0]; sout[1] = conformal_spatial_metric_source[0][1]; sout[2] = conformal_spatial_metric_source[0][2];
    sout[3] = conformal_spatial_metric_source[1][0]; sout[4] = conformal_spatial_metric_source[1][1]; sout[5] = conformal_spatial_metric_source[1][2];
    sout[6] = conformal_spatial_metric_source[2][0]; sout[7] = conformal_spatial_metric_source[2][1]; sout[8] = conformal_spatial_metric_source[2][2];

    sout[9] = conformal_lapse_source;

    sout[10] = conformal_extrinsic_curvature_source[0][0]; sout[11] = conformal_extrinsic_curvature_source[0][1]; sout[12] = conformal_extrinsic_curvature_source[0][2];
    sout[13] = conformal_extrinsic_curvature_source[1][0]; sout[14] = conformal_extrinsic_curvature_source[1][1]; sout[15] = conformal_extrinsic_curvature_source[1][2];
    sout[16] = conformal_extrinsic_curvature_source[2][0]; sout[17] = conformal_extrinsic_curvature_source[2][1]; sout[18] = conformal_extrinsic_curvature_source[2][2];

    for (int i = 19; i < 49; i++) {
      sout[i] = 0.0;
    }

    sout[49] = conformal_aux_vect_source[0];
    sout[50] = conformal_aux_vect_source[1];
    sout[51] = conformal_aux_vect_source[2];

    for (int i = 52; i < 77; i++) {
      sout[i] = 0.0;
    }

    for (int i = 0; i < 3; i++) {
      gkyl_free(inv_conformal_spatial_metric[i]);
    }
    gkyl_free(inv_conformal_spatial_metric);
  }
  else {
    for (int i = 0; i < 77; i++) {
      sout[i] = 0.0;
    }
  }
}

void
gkyl_vacuum_einstein_conformal_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_vacuum_einstein_conformal *vacuum_einstein_conformal = container_of(base->on_dev, struct wv_vacuum_einstein_conformal, eqn);
    gkyl_cu_free(vacuum_einstein_conformal);
  }

  struct wv_vacuum_einstein_conformal *vacuum_einstein_conformal = container_of(base, struct wv_vacuum_einstein_conformal, eqn);
  gkyl_free(vacuum_einstein_conformal);
}

struct gkyl_wv_eqn*
gkyl_wv_vacuum_einstein_conformal_new(double excision_threshold, enum gkyl_spacetime_slicing spacetime_slicing, enum gkyl_spacetime_evolution spacetime_evolution, bool use_gpu)
{
  return gkyl_wv_vacuum_einstein_conformal_inew(&(struct gkyl_wv_vacuum_einstein_conformal_inp) {
      .excision_threshold = excision_threshold,
      .spacetime_slicing = spacetime_slicing,
      .spacetime_evolution = spacetime_evolution,
      .rp_type = WV_VACUUM_EINSTEIN_CONFORMAL_RP_HLL,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_vacuum_einstein_conformal_inew(const struct gkyl_wv_vacuum_einstein_conformal_inp* inp)
{
  struct wv_vacuum_einstein_conformal *vacuum_einstein_conformal = gkyl_malloc(sizeof(struct wv_vacuum_einstein_conformal));

  vacuum_einstein_conformal->eqn.type = GKYL_EQN_VACUUM_EINSTEIN_CONFORMAL;
  vacuum_einstein_conformal->eqn.num_equations = 77;
  vacuum_einstein_conformal->eqn.num_diag = 1;

  vacuum_einstein_conformal->excision_threshold = inp->excision_threshold;
  vacuum_einstein_conformal->spacetime_slicing = inp->spacetime_slicing;
  vacuum_einstein_conformal->spacetime_evolution = inp->spacetime_evolution;

  if (inp->rp_type == WV_VACUUM_EINSTEIN_CONFORMAL_RP_LAX) {
    vacuum_einstein_conformal->eqn.num_waves = 2;
    vacuum_einstein_conformal->eqn.waves_func = wave_lax_l;
    vacuum_einstein_conformal->eqn.qfluct_func = qfluct_lax_l;
  }
  else if (inp->rp_type == WV_VACUUM_EINSTEIN_CONFORMAL_RP_HLL) {
    vacuum_einstein_conformal->eqn.num_waves = 2;
    vacuum_einstein_conformal->eqn.waves_func = wave_hll_l;
    vacuum_einstein_conformal->eqn.qfluct_func = qfluct_hll_l;
  }
  
  vacuum_einstein_conformal->eqn.flux_jump = flux_jump;
  vacuum_einstein_conformal->eqn.check_inv_func = check_inv;
  vacuum_einstein_conformal->eqn.max_speed_func = max_speed;
  vacuum_einstein_conformal->eqn.rotate_to_local_func = rot_to_local;
  vacuum_einstein_conformal->eqn.rotate_to_global_func = rot_to_global;

  vacuum_einstein_conformal->eqn.wall_bc_func = vacuum_einstein_conformal_wall;

  vacuum_einstein_conformal->eqn.cons_to_riem = cons_to_riem;
  vacuum_einstein_conformal->eqn.riem_to_cons = riem_to_cons;

  vacuum_einstein_conformal->eqn.cons_to_diag = vacuum_einstein_conformal_cons_to_diag;

  vacuum_einstein_conformal->eqn.source_func = vacuum_einstein_conformal_source;

  vacuum_einstein_conformal->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(vacuum_einstein_conformal->eqn.flags);
  vacuum_einstein_conformal->eqn.ref_count = gkyl_ref_count_init(gkyl_vacuum_einstein_conformal_free);
  vacuum_einstein_conformal->eqn.on_dev = &vacuum_einstein_conformal->eqn; // On the CPU, the equation object points to itself.

  return &vacuum_einstein_conformal->eqn;
}

double
gkyl_wv_vacuum_einstein_conformal_excision_threshold(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_vacuum_einstein_conformal *vacuum_einstein_conformal = container_of(eqn, struct wv_vacuum_einstein_conformal, eqn);
  double excision_threshold = vacuum_einstein_conformal->excision_threshold;

  return excision_threshold;
}

enum gkyl_spacetime_slicing
gkyl_wv_vacuum_einstein_conformal_spacetime_slicing(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_vacuum_einstein_conformal *vacuum_einstein_conformal = container_of(eqn, struct wv_vacuum_einstein_conformal, eqn);
  enum gkyl_spacetime_slicing spacetime_slicing = vacuum_einstein_conformal->spacetime_slicing;

  return spacetime_slicing;
}

enum gkyl_spacetime_evolution
gkyl_wv_vacuum_einstein_conformal_spacetime_evolution(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_vacuum_einstein_conformal *vacuum_einstein_conformal = container_of(eqn, struct wv_vacuum_einstein_conformal, eqn);
  enum gkyl_spacetime_evolution spacetime_evolution = vacuum_einstein_conformal->spacetime_evolution;

  return spacetime_evolution;
}