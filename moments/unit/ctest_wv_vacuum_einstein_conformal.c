#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_vacuum_einstein_conformal.h>
#include <gkyl_wv_vacuum_einstein_conformal_priv.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_gr_blackhole.h>

void
test_vacuum_einstein_conformal_basic_minkowski()
{
  double excision_threshold = 0.3;
  enum gkyl_spacetime_slicing spacetime_slicing = GKYL_HARMONIC_SLICING;
  enum gkyl_spacetime_evolution spacetime_evolution = GKYL_EINSTEIN_EVOLUTION;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  struct gkyl_wv_eqn *vacuum_einstein_conformal = gkyl_wv_vacuum_einstein_conformal_new(excision_threshold, spacetime_slicing, spacetime_evolution, false);

  TEST_CHECK( vacuum_einstein_conformal->num_equations == 77 );
  TEST_CHECK( vacuum_einstein_conformal->num_waves == 2 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double conformal_spatial_det, conformal_lapse;
      double *conformal_shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **conformal_spatial_metric = gkyl_malloc(sizeof(double*[3]));
      double **inv_conformal_spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        conformal_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
        inv_conformal_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **conformal_extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        conformal_extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *conformal_lapse_der = gkyl_malloc(sizeof(double[3]));
      double **conformal_shift_der = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        conformal_shift_der[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***conformal_spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        conformal_spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          conformal_spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      double conformal_fact, bssn_conformal_fact;
      double *conformal_fact_der = gkyl_malloc(sizeof(double[3]));
      double *bssn_conformal_fact_der = gkyl_malloc(sizeof(double[3]));
      double **bssn_conformal_fact_der2 = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        bssn_conformal_fact_der2[i] = gkyl_malloc(sizeof(double[3]));
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &conformal_spatial_det);
      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &conformal_lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &conformal_shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
      
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &conformal_spatial_metric);
      spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_conformal_spatial_metric);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &conformal_extrinsic_curvature);

      spacetime->conformal_factor_func(spacetime, 0.0, x, y, 0.0, &conformal_fact);
      spacetime->bssn_conformal_factor_func(spacetime, 0.0, x, y, 0.0, &bssn_conformal_fact);

      spacetime->conformal_factor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &conformal_fact_der);
      spacetime->bssn_conformal_factor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &bssn_conformal_fact_der);
      spacetime->bssn_conformal_factor_der2_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &bssn_conformal_fact_der2);

      spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &conformal_lapse_der);
      spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &conformal_shift_der);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &conformal_spatial_metric_der);

      for (int i = 0; i < 3; i++) {
        conformal_fact_der[i] /= conformal_fact;
      }

      // Set first and second conformal derivatives to zero, to improve stability.
      for (int i = 0; i < 3; i++) {
        conformal_fact_der[i] = 0.0;
        bssn_conformal_fact_der[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          bssn_conformal_fact_der2[i][j] = 0.0;
        }
      }

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          conformal_spatial_metric[i][j] /= (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
          inv_conformal_spatial_metric[i][j] *= (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
        }
      }

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            conformal_spatial_metric_der[i][j][k] = 0.5 * conformal_spatial_metric_der[i][j][k];
            conformal_spatial_metric_der[i][j][k] /= (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
            conformal_spatial_metric_der[i][j][k] -= 2.0 * conformal_fact_der[i] * conformal_spatial_metric[j][k];
          }

          conformal_shift_der[i][j] = 0.5 * conformal_shift_der[i][j];
        }
      }

      for (int i = 0; i < 3; i++) {
        conformal_lapse_der[i] = conformal_lapse_der[i] / conformal_lapse;
      }

      double conformal_extrinsic_curvature_trace = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          conformal_extrinsic_curvature_trace += inv_conformal_spatial_metric[i][j] * conformal_extrinsic_curvature[i][j];
        }
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

      double conformal_aux_vect[3];
      for (int i = 0; i < 3; i++) {
        conformal_aux_vect[i] = 0.0;

        for (int s = 0; s < 3; s++) {
          conformal_aux_vect[i] += conformal_spatial_metric_der_raised3[i][s][s];
          conformal_aux_vect[i] -= conformal_spatial_metric_der_raised1[s][s][i];
        }

        conformal_aux_vect[i] -= 4.0 * conformal_fact_der[i];
      }

      double conformal_aux_vect_raised[3];
      for (int k = 0; k < 3; k++) {
        conformal_aux_vect_raised[k] = 0.0;
            
        for (int l = 0; l < 3; l++) {
          conformal_aux_vect_raised[k] += inv_conformal_spatial_metric[k][l] * conformal_aux_vect[l];
        }
      }

      double conformal_shift_der_lowered[3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          conformal_shift_der_lowered[i][j] = 0.0;

          for (int k = 0; k < 3; k++) {
            conformal_shift_der_lowered[i][j] += conformal_spatial_metric[k][j] * conformal_shift_der[i][k];
          }
        }
      }

      double conformal_shift_der_switched[3][3];
      for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
          conformal_shift_der_switched[i][k] = 0.0;

          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              conformal_shift_der_switched[i][k] += inv_conformal_spatial_metric[i][l] * conformal_spatial_metric[m][k] * conformal_shift_der[l][m];
            }
          }
        }
      }

      double symmetrized_conformal_shift[3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          symmetrized_conformal_shift[i][j] = (1.0 / conformal_lapse) * (conformal_shift_der_lowered[i][j] + conformal_shift_der_lowered[j][i]);
        }
      }

      double q[77];
      q[0] = conformal_spatial_metric[0][0]; q[1] = conformal_spatial_metric[0][1]; q[2] = conformal_spatial_metric[0][2];
      q[3] = conformal_spatial_metric[1][0]; q[4] = conformal_spatial_metric[1][1]; q[5] = conformal_spatial_metric[1][2];
      q[6] = conformal_spatial_metric[2][0]; q[7] = conformal_spatial_metric[2][1]; q[8] = conformal_spatial_metric[2][2];

      q[9] = conformal_lapse;

      q[10] = conformal_extrinsic_curvature[0][0]; q[11] = conformal_extrinsic_curvature[0][1]; q[12] = conformal_extrinsic_curvature[0][2];
      q[13] = conformal_extrinsic_curvature[1][0]; q[14] = conformal_extrinsic_curvature[1][1]; q[15] = conformal_extrinsic_curvature[1][2];
      q[16] = conformal_extrinsic_curvature[2][0]; q[17] = conformal_extrinsic_curvature[2][1]; q[18] = conformal_extrinsic_curvature[2][2];

      q[19] = conformal_spatial_metric_der[0][0][0]; q[20] = conformal_spatial_metric_der[0][0][1]; q[21] = conformal_spatial_metric_der[0][0][2];
      q[22] = conformal_spatial_metric_der[0][1][0]; q[23] = conformal_spatial_metric_der[0][1][1]; q[24] = conformal_spatial_metric_der[0][1][2];
      q[25] = conformal_spatial_metric_der[0][2][0]; q[26] = conformal_spatial_metric_der[0][2][1]; q[27] = conformal_spatial_metric_der[0][2][2];

      q[28] = conformal_spatial_metric_der[1][0][0]; q[29] = conformal_spatial_metric_der[1][0][1]; q[30] = conformal_spatial_metric_der[1][0][2];
      q[31] = conformal_spatial_metric_der[1][1][0]; q[32] = conformal_spatial_metric_der[1][1][1]; q[33] = conformal_spatial_metric_der[1][1][2];
      q[34] = conformal_spatial_metric_der[1][2][0]; q[35] = conformal_spatial_metric_der[1][2][1]; q[36] = conformal_spatial_metric_der[1][2][2];

      q[37] = conformal_spatial_metric_der[2][0][0]; q[38] = conformal_spatial_metric_der[2][0][1]; q[39] = conformal_spatial_metric_der[2][0][2];
      q[40] = conformal_spatial_metric_der[2][1][0]; q[41] = conformal_spatial_metric_der[2][1][1]; q[42] = conformal_spatial_metric_der[2][1][2];
      q[43] = conformal_spatial_metric_der[2][2][0]; q[44] = conformal_spatial_metric_der[2][2][1]; q[45] = conformal_spatial_metric_der[2][2][2];

      q[46] = conformal_lapse_der[0]; q[47] = conformal_lapse_der[1]; q[48] = conformal_lapse_der[2];

      q[49] = conformal_aux_vect[0]; q[50] = conformal_aux_vect[1]; q[51] = conformal_aux_vect[2];

      q[52] = conformal_shift[0]; q[53] = conformal_shift[1]; q[54] = conformal_shift[2];
      
      q[55] = conformal_shift_der[0][0]; q[56] = conformal_shift_der[0][1]; q[57] = conformal_shift_der[0][2];
      q[58] = conformal_shift_der[1][0]; q[59] = conformal_shift_der[1][1]; q[60] = conformal_shift_der[1][2];
      q[61] = conformal_shift_der[2][0]; q[62] = conformal_shift_der[2][1]; q[63] = conformal_shift_der[2][2];

      q[64] = bssn_conformal_fact;

      q[65] = bssn_conformal_fact_der[0]; q[66] = bssn_conformal_fact_der[1]; q[67] = bssn_conformal_fact_der[2];

      q[68] = bssn_conformal_fact_der2[0][0]; q[69] = bssn_conformal_fact_der2[0][1]; q[70] = bssn_conformal_fact_der2[0][2];
      q[71] = bssn_conformal_fact_der2[1][0]; q[72] = bssn_conformal_fact_der2[1][1]; q[73] = bssn_conformal_fact_der2[1][2];
      q[74] = bssn_conformal_fact_der2[2][0]; q[75] = bssn_conformal_fact_der2[2][1]; q[76] = bssn_conformal_fact_der2[2][2];

      double evolution_func = 1.0;
      double slicing_func = conformal_extrinsic_curvature_trace / (conformal_fact * conformal_fact * conformal_fact * conformal_fact);

      double conformal_extrinsic_curvature_flux[3][3][3];
      for (int d = 0; d < 3; d++) {
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            conformal_extrinsic_curvature_flux[d][i][j] = -conformal_shift[d] * conformal_extrinsic_curvature[i][j];
            conformal_extrinsic_curvature_flux[d][i][j] += conformal_lapse * conformal_spatial_metric_der_raised1[d][i][j];
            conformal_extrinsic_curvature_flux[d][i][j] -= conformal_lapse * (0.5 * evolution_func) * conformal_aux_vect_raised[d] * conformal_spatial_metric[i][j];

            if (i == d) {
              conformal_extrinsic_curvature_flux[d][i][j] += 0.5 * conformal_lapse * conformal_lapse_der[j];
              conformal_extrinsic_curvature_flux[d][i][j] += conformal_lapse * conformal_aux_vect[j];
              for (int r = 0; r < 3; r++) {
                conformal_extrinsic_curvature_flux[d][i][j] -= 0.5 * conformal_lapse * conformal_spatial_metric_der_raised3[j][r][r];
              }
            }

            if (j == d) {
              conformal_extrinsic_curvature_flux[d][i][j] += 0.5 * conformal_lapse * conformal_lapse_der[i];
              conformal_extrinsic_curvature_flux[d][i][j] += conformal_lapse * conformal_aux_vect[i];
              for (int r = 0; r < 3; r++) {
                conformal_extrinsic_curvature_flux[d][i][j] -= 0.5 * conformal_lapse * conformal_spatial_metric_der_raised3[i][r][r];
              }
            }
          }
        }
      }

      double conformal_spatial_metric_der_flux[3][3][3][3];
      for (int d = 0; d < 3; d++) {
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              conformal_spatial_metric_der_flux[d][k][i][j] = 0.0;
            }

            for (int r = 0; r < 3; r++) {
              conformal_spatial_metric_der_flux[d][d][i][j] -= conformal_shift[r] * conformal_spatial_metric_der[r][i][j];
              conformal_spatial_metric_der_flux[d][d][i][j] -= 2.0 * conformal_shift[r] * conformal_fact_der[r] * conformal_spatial_metric[i][j];
            }

            conformal_spatial_metric_der_flux[d][d][i][j] += conformal_lapse * (conformal_extrinsic_curvature[i][j] /
              (conformal_fact * conformal_fact * conformal_fact * conformal_fact));
            conformal_spatial_metric_der_flux[d][d][i][j] -= conformal_lapse * symmetrized_conformal_shift[i][j];
          }
        }
      }

      double conformal_lapse_der_flux[3][3];
      for (int d = 0; d < 3; d++) {
        for (int i = 0; i < 3; i++) {
          conformal_lapse_der_flux[d][i] = 0.0;
        }
      }
      
      for (int d = 0; d < 3; d++) {
        for (int r = 0; r < 3; r++) {
          conformal_lapse_der_flux[d][d] -= conformal_shift[r] * conformal_lapse_der[r];
        }
        conformal_lapse_der_flux[d][d] += conformal_lapse * slicing_func;
      }

      double conformal_aux_vect_flux[3][3];
      for (int d = 0; d < 3; d++) {
        for (int i = 0; i < 3; i++) {
          conformal_aux_vect_flux[d][i] = -conformal_shift[d] * conformal_aux_vect[i];
          conformal_aux_vect_flux[d][i] -= 4.0 * conformal_shift[d] * conformal_fact_der[i];
          conformal_aux_vect_flux[d][i] += conformal_shift_der_switched[d][i];
          conformal_aux_vect_flux[d][i] -= conformal_shift_der[i][d];
        }
      }

      double fluxes[3][42] = {
        { conformal_extrinsic_curvature_flux[0][0][0], conformal_extrinsic_curvature_flux[0][0][1], conformal_extrinsic_curvature_flux[0][0][2],
          conformal_extrinsic_curvature_flux[0][1][0], conformal_extrinsic_curvature_flux[0][1][1], conformal_extrinsic_curvature_flux[0][1][2],
          conformal_extrinsic_curvature_flux[0][2][0], conformal_extrinsic_curvature_flux[0][2][1], conformal_extrinsic_curvature_flux[0][2][2],
          conformal_spatial_metric_der_flux[0][0][0][0], conformal_spatial_metric_der_flux[0][0][0][1], conformal_spatial_metric_der_flux[0][0][0][2],
          conformal_spatial_metric_der_flux[0][0][1][0], conformal_spatial_metric_der_flux[0][0][1][1], conformal_spatial_metric_der_flux[0][0][1][2],
          conformal_spatial_metric_der_flux[0][0][2][0], conformal_spatial_metric_der_flux[0][0][2][1], conformal_spatial_metric_der_flux[0][0][2][2],
          conformal_spatial_metric_der_flux[0][1][0][0], conformal_spatial_metric_der_flux[0][1][0][1], conformal_spatial_metric_der_flux[0][1][0][2],
          conformal_spatial_metric_der_flux[0][1][1][0], conformal_spatial_metric_der_flux[0][1][1][1], conformal_spatial_metric_der_flux[0][1][1][2],
          conformal_spatial_metric_der_flux[0][1][2][0], conformal_spatial_metric_der_flux[0][1][2][1], conformal_spatial_metric_der_flux[0][1][2][2],
          conformal_spatial_metric_der_flux[0][2][0][0], conformal_spatial_metric_der_flux[0][2][0][1], conformal_spatial_metric_der_flux[0][2][0][2],
          conformal_spatial_metric_der_flux[0][2][1][0], conformal_spatial_metric_der_flux[0][2][1][1], conformal_spatial_metric_der_flux[0][2][1][2],
          conformal_spatial_metric_der_flux[0][2][2][0], conformal_spatial_metric_der_flux[0][2][2][1], conformal_spatial_metric_der_flux[0][2][2][2],
          conformal_lapse_der_flux[0][0], conformal_lapse_der_flux[0][1], conformal_lapse_der_flux[0][2],
          conformal_aux_vect_flux[0][0], conformal_aux_vect_flux[0][1], conformal_aux_vect_flux[0][2] },
        { conformal_extrinsic_curvature_flux[1][0][0], conformal_extrinsic_curvature_flux[1][0][1], conformal_extrinsic_curvature_flux[1][0][2],
          conformal_extrinsic_curvature_flux[1][1][0], conformal_extrinsic_curvature_flux[1][1][1], conformal_extrinsic_curvature_flux[1][1][2],
          conformal_extrinsic_curvature_flux[1][2][0], conformal_extrinsic_curvature_flux[1][2][1], conformal_extrinsic_curvature_flux[1][2][2],
          conformal_spatial_metric_der_flux[1][0][0][0], conformal_spatial_metric_der_flux[1][0][0][1], conformal_spatial_metric_der_flux[1][0][0][2],
          conformal_spatial_metric_der_flux[1][0][1][0], conformal_spatial_metric_der_flux[1][0][1][1], conformal_spatial_metric_der_flux[1][0][1][2],
          conformal_spatial_metric_der_flux[1][0][2][0], conformal_spatial_metric_der_flux[1][0][2][1], conformal_spatial_metric_der_flux[1][0][2][2],
          conformal_spatial_metric_der_flux[1][1][0][0], conformal_spatial_metric_der_flux[1][1][0][1], conformal_spatial_metric_der_flux[1][1][0][2],
          conformal_spatial_metric_der_flux[1][1][1][0], conformal_spatial_metric_der_flux[1][1][1][1], conformal_spatial_metric_der_flux[1][1][1][2],
          conformal_spatial_metric_der_flux[1][1][2][0], conformal_spatial_metric_der_flux[1][1][2][1], conformal_spatial_metric_der_flux[1][1][2][2],
          conformal_spatial_metric_der_flux[1][2][0][0], conformal_spatial_metric_der_flux[1][2][0][1], conformal_spatial_metric_der_flux[1][2][0][2],
          conformal_spatial_metric_der_flux[1][2][1][0], conformal_spatial_metric_der_flux[1][2][1][1], conformal_spatial_metric_der_flux[1][2][1][2],
          conformal_spatial_metric_der_flux[1][2][2][0], conformal_spatial_metric_der_flux[1][2][2][1], conformal_spatial_metric_der_flux[1][2][2][2],
          conformal_lapse_der_flux[1][0], conformal_lapse_der_flux[1][1], conformal_lapse_der_flux[1][2],
          conformal_aux_vect_flux[1][0], conformal_aux_vect_flux[1][1], conformal_aux_vect_flux[1][2] },
        { conformal_extrinsic_curvature_flux[2][0][0], conformal_extrinsic_curvature_flux[2][0][1], conformal_extrinsic_curvature_flux[2][0][2],
          conformal_extrinsic_curvature_flux[2][1][0], conformal_extrinsic_curvature_flux[2][1][1], conformal_extrinsic_curvature_flux[2][1][2],
          conformal_extrinsic_curvature_flux[2][2][0], conformal_extrinsic_curvature_flux[2][2][1], conformal_extrinsic_curvature_flux[2][2][2],
          conformal_spatial_metric_der_flux[2][0][0][0], conformal_spatial_metric_der_flux[2][0][0][1], conformal_spatial_metric_der_flux[2][0][0][2],
          conformal_spatial_metric_der_flux[2][0][1][0], conformal_spatial_metric_der_flux[2][0][1][1], conformal_spatial_metric_der_flux[2][0][1][2],
          conformal_spatial_metric_der_flux[2][0][2][0], conformal_spatial_metric_der_flux[2][0][2][1], conformal_spatial_metric_der_flux[2][0][2][2],
          conformal_spatial_metric_der_flux[2][1][0][0], conformal_spatial_metric_der_flux[2][1][0][1], conformal_spatial_metric_der_flux[2][1][0][2],
          conformal_spatial_metric_der_flux[2][1][1][0], conformal_spatial_metric_der_flux[2][1][1][1], conformal_spatial_metric_der_flux[2][1][1][2],
          conformal_spatial_metric_der_flux[2][1][2][0], conformal_spatial_metric_der_flux[2][1][2][1], conformal_spatial_metric_der_flux[2][1][2][2],
          conformal_spatial_metric_der_flux[2][2][0][0], conformal_spatial_metric_der_flux[2][2][0][1], conformal_spatial_metric_der_flux[2][2][0][2],
          conformal_spatial_metric_der_flux[2][2][1][0], conformal_spatial_metric_der_flux[2][2][1][1], conformal_spatial_metric_der_flux[2][2][1][2],
          conformal_spatial_metric_der_flux[2][2][2][0], conformal_spatial_metric_der_flux[2][2][2][1], conformal_spatial_metric_der_flux[2][2][2][2],
          conformal_lapse_der_flux[2][0], conformal_lapse_der_flux[2][1], conformal_lapse_der_flux[2][2],
          conformal_aux_vect_flux[2][0], conformal_aux_vect_flux[2][1], conformal_aux_vect_flux[2][2] },
      };

      double norm[3][3] = {
        { 1.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 },
      };

      double tau1[3][3] = {
        { 0.0, 1.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0 },
      };

      double tau2[3][3] = {
        { 0.0, 0.0, 1.0 },
        { 0.0, 0.0, -1.0 },
        { 0.0, 1.0, 0.0 },
      };

      double q_local[77], flux_local[77], flux[77];
      for (int d = 0; d < 3; d++) {
        vacuum_einstein_conformal->rotate_to_local_func(vacuum_einstein_conformal, tau1[d], tau2[d], norm[d], q, q_local);
        gkyl_vacuum_einstein_conformal_flux(excision_threshold, spacetime_slicing, spacetime_evolution, q_local, flux_local);
        vacuum_einstein_conformal->rotate_to_global_func(vacuum_einstein_conformal, tau1[d], tau2[d], norm[d], flux_local, flux);

        for (int i = 0; i < 42; i++) {
          TEST_CHECK( gkyl_compare(flux[i + 10], fluxes[d][i], 1e-8) );
        }
      }
      
      double q_l[77], q_g[77];
      for (int d = 0; d < 3; d++) {
        gkyl_wv_eqn_rotate_to_local(vacuum_einstein_conformal, tau1[d], tau2[d], norm[d], q, q_l);
        gkyl_wv_eqn_rotate_to_global(vacuum_einstein_conformal, tau1[d], tau2[d], norm[d], q_l, q_g);

        for (int i = 0; i < 77; i++) {
          TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
        }

        double w1[77], q1[77];
        vacuum_einstein_conformal->cons_to_riem(vacuum_einstein_conformal, q_local, q_local, w1);
        vacuum_einstein_conformal->riem_to_cons(vacuum_einstein_conformal, q_local, w1, q1);

        for (int i = 0; i < 77; i++) {
          TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(conformal_spatial_metric[i]);
        gkyl_free(inv_conformal_spatial_metric[i]);
        gkyl_free(conformal_extrinsic_curvature[i]);
        gkyl_free(conformal_shift_der[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(conformal_spatial_metric_der[i][j]);
        }
        gkyl_free(conformal_spatial_metric_der[i]);
        gkyl_free(bssn_conformal_fact_der2[i]);
      }
      gkyl_free(conformal_spatial_metric);
      gkyl_free(inv_conformal_spatial_metric);
      gkyl_free(conformal_extrinsic_curvature);
      gkyl_free(conformal_shift);
      gkyl_free(conformal_lapse_der);
      gkyl_free(conformal_shift_der);
      gkyl_free(conformal_spatial_metric_der);
      gkyl_free(conformal_fact_der);
      gkyl_free(bssn_conformal_fact_der);
      gkyl_free(bssn_conformal_fact_der2);
    }
  }

  gkyl_wv_eqn_release(vacuum_einstein_conformal);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_vacuum_einstein_conformal_basic_schwarzschild()
{
  double excision_threshold = 0.3;
  enum gkyl_spacetime_slicing spacetime_slicing = GKYL_1PLUSLOG_SLICING;
  enum gkyl_spacetime_evolution spacetime_evolution = GKYL_EINSTEIN_EVOLUTION;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *vacuum_einstein_conformal = gkyl_wv_vacuum_einstein_conformal_new(excision_threshold, spacetime_slicing, spacetime_evolution, false);

  TEST_CHECK( vacuum_einstein_conformal->num_equations == 77 );
  TEST_CHECK( vacuum_einstein_conformal->num_waves == 2 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double conformal_spatial_det, conformal_lapse;
      double *conformal_shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **conformal_spatial_metric = gkyl_malloc(sizeof(double*[3]));
      double **inv_conformal_spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        conformal_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
        inv_conformal_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **conformal_extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        conformal_extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *conformal_lapse_der = gkyl_malloc(sizeof(double[3]));
      double **conformal_shift_der = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        conformal_shift_der[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***conformal_spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        conformal_spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          conformal_spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      double conformal_fact, bssn_conformal_fact;
      double *conformal_fact_der = gkyl_malloc(sizeof(double[3]));
      double *bssn_conformal_fact_der = gkyl_malloc(sizeof(double[3]));
      double **bssn_conformal_fact_der2 = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        bssn_conformal_fact_der2[i] = gkyl_malloc(sizeof(double[3]));
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &conformal_spatial_det);
      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &conformal_lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &conformal_shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
      
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &conformal_spatial_metric);
      spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_conformal_spatial_metric);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &conformal_extrinsic_curvature);

      spacetime->conformal_factor_func(spacetime, 0.0, x, y, 0.0, &conformal_fact);
      spacetime->bssn_conformal_factor_func(spacetime, 0.0, x, y, 0.0, &bssn_conformal_fact);

      spacetime->conformal_factor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &conformal_fact_der);
      spacetime->bssn_conformal_factor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &bssn_conformal_fact_der);
      spacetime->bssn_conformal_factor_der2_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &bssn_conformal_fact_der2);

      spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &conformal_lapse_der);
      spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &conformal_shift_der);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &conformal_spatial_metric_der);

      for (int i = 0; i < 3; i++) {
        conformal_fact_der[i] /= conformal_fact;
      }

      // Set first and second conformal derivatives to zero, to improve stability.
      for (int i = 0; i < 3; i++) {
        conformal_fact_der[i] = 0.0;
        bssn_conformal_fact_der[i] = 0.0;

        for (int j = 0; j < 3; j++) {
          bssn_conformal_fact_der2[i][j] = 0.0;
        }
      }

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          conformal_spatial_metric[i][j] /= (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
          inv_conformal_spatial_metric[i][j] *= (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
        }
      }

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            conformal_spatial_metric_der[i][j][k] = 0.5 * conformal_spatial_metric_der[i][j][k];
            conformal_spatial_metric_der[i][j][k] /= (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
            conformal_spatial_metric_der[i][j][k] -= 2.0 * conformal_fact_der[i] * conformal_spatial_metric[j][k];
          }

          conformal_shift_der[i][j] = 0.5 * conformal_shift_der[i][j];
        }
      }

      for (int i = 0; i < 3; i++) {
        conformal_lapse_der[i] = conformal_lapse_der[i] / conformal_lapse;
      }

      double conformal_extrinsic_curvature_trace = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          conformal_extrinsic_curvature_trace += inv_conformal_spatial_metric[i][j] * conformal_extrinsic_curvature[i][j];
        }
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

      double conformal_aux_vect[3];
      for (int i = 0; i < 3; i++) {
        conformal_aux_vect[i] = 0.0;

        for (int s = 0; s < 3; s++) {
          conformal_aux_vect[i] += conformal_spatial_metric_der_raised3[i][s][s];
          conformal_aux_vect[i] -= conformal_spatial_metric_der_raised1[s][s][i];
        }

        conformal_aux_vect[i] -= 4.0 * conformal_fact_der[i];
      }

      double conformal_aux_vect_raised[3];
      for (int k = 0; k < 3; k++) {
        conformal_aux_vect_raised[k] = 0.0;
            
        for (int l = 0; l < 3; l++) {
          conformal_aux_vect_raised[k] += inv_conformal_spatial_metric[k][l] * conformal_aux_vect[l];
        }
      }

      double conformal_shift_der_lowered[3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          conformal_shift_der_lowered[i][j] = 0.0;

          for (int k = 0; k < 3; k++) {
            conformal_shift_der_lowered[i][j] += conformal_spatial_metric[k][j] * conformal_shift_der[i][k];
          }
        }
      }

      double conformal_shift_der_switched[3][3];
      for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
          conformal_shift_der_switched[i][k] = 0.0;

          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              conformal_shift_der_switched[i][k] += inv_conformal_spatial_metric[i][l] * conformal_spatial_metric[m][k] * conformal_shift_der[l][m];
            }
          }
        }
      }

      double symmetrized_conformal_shift[3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          symmetrized_conformal_shift[i][j] = (1.0 / conformal_lapse) * (conformal_shift_der_lowered[i][j] + conformal_shift_der_lowered[j][i]);
        }
      }

      if (!in_excision_region) {
        double q[77];
        q[0] = conformal_spatial_metric[0][0]; q[1] = conformal_spatial_metric[0][1]; q[2] = conformal_spatial_metric[0][2];
        q[3] = conformal_spatial_metric[1][0]; q[4] = conformal_spatial_metric[1][1]; q[5] = conformal_spatial_metric[1][2];
        q[6] = conformal_spatial_metric[2][0]; q[7] = conformal_spatial_metric[2][1]; q[8] = conformal_spatial_metric[2][2];

        q[9] = conformal_lapse;

        q[10] = conformal_extrinsic_curvature[0][0]; q[11] = conformal_extrinsic_curvature[0][1]; q[12] = conformal_extrinsic_curvature[0][2];
        q[13] = conformal_extrinsic_curvature[1][0]; q[14] = conformal_extrinsic_curvature[1][1]; q[15] = conformal_extrinsic_curvature[1][2];
        q[16] = conformal_extrinsic_curvature[2][0]; q[17] = conformal_extrinsic_curvature[2][1]; q[18] = conformal_extrinsic_curvature[2][2];

        q[19] = conformal_spatial_metric_der[0][0][0]; q[20] = conformal_spatial_metric_der[0][0][1]; q[21] = conformal_spatial_metric_der[0][0][2];
        q[22] = conformal_spatial_metric_der[0][1][0]; q[23] = conformal_spatial_metric_der[0][1][1]; q[24] = conformal_spatial_metric_der[0][1][2];
        q[25] = conformal_spatial_metric_der[0][2][0]; q[26] = conformal_spatial_metric_der[0][2][1]; q[27] = conformal_spatial_metric_der[0][2][2];

        q[28] = conformal_spatial_metric_der[1][0][0]; q[29] = conformal_spatial_metric_der[1][0][1]; q[30] = conformal_spatial_metric_der[1][0][2];
        q[31] = conformal_spatial_metric_der[1][1][0]; q[32] = conformal_spatial_metric_der[1][1][1]; q[33] = conformal_spatial_metric_der[1][1][2];
        q[34] = conformal_spatial_metric_der[1][2][0]; q[35] = conformal_spatial_metric_der[1][2][1]; q[36] = conformal_spatial_metric_der[1][2][2];

        q[37] = conformal_spatial_metric_der[2][0][0]; q[38] = conformal_spatial_metric_der[2][0][1]; q[39] = conformal_spatial_metric_der[2][0][2];
        q[40] = conformal_spatial_metric_der[2][1][0]; q[41] = conformal_spatial_metric_der[2][1][1]; q[42] = conformal_spatial_metric_der[2][1][2];
        q[43] = conformal_spatial_metric_der[2][2][0]; q[44] = conformal_spatial_metric_der[2][2][1]; q[45] = conformal_spatial_metric_der[2][2][2];

        q[46] = conformal_lapse_der[0]; q[47] = conformal_lapse_der[1]; q[48] = conformal_lapse_der[2];

        q[49] = conformal_aux_vect[0]; q[50] = conformal_aux_vect[1]; q[51] = conformal_aux_vect[2];

        q[52] = conformal_shift[0]; q[53] = conformal_shift[1]; q[54] = conformal_shift[2];
        
        q[55] = conformal_shift_der[0][0]; q[56] = conformal_shift_der[0][1]; q[57] = conformal_shift_der[0][2];
        q[58] = conformal_shift_der[1][0]; q[59] = conformal_shift_der[1][1]; q[60] = conformal_shift_der[1][2];
        q[61] = conformal_shift_der[2][0]; q[62] = conformal_shift_der[2][1]; q[63] = conformal_shift_der[2][2];

        q[64] = bssn_conformal_fact;

        q[65] = bssn_conformal_fact_der[0]; q[66] = bssn_conformal_fact_der[1]; q[67] = bssn_conformal_fact_der[2];

        q[68] = bssn_conformal_fact_der2[0][0]; q[69] = bssn_conformal_fact_der2[0][1]; q[70] = bssn_conformal_fact_der2[0][2];
        q[71] = bssn_conformal_fact_der2[1][0]; q[72] = bssn_conformal_fact_der2[1][1]; q[73] = bssn_conformal_fact_der2[1][2];
        q[74] = bssn_conformal_fact_der2[2][0]; q[75] = bssn_conformal_fact_der2[2][1]; q[76] = bssn_conformal_fact_der2[2][2];

        double evolution_func = 1.0;
        double slicing_func = 2.0 * conformal_extrinsic_curvature_trace / (conformal_lapse * conformal_fact * conformal_fact * conformal_fact * conformal_fact);

        double conformal_extrinsic_curvature_flux[3][3][3];
        for (int d = 0; d < 3; d++) {
          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
              conformal_extrinsic_curvature_flux[d][i][j] = -conformal_shift[d] * conformal_extrinsic_curvature[i][j];
              conformal_extrinsic_curvature_flux[d][i][j] += conformal_lapse * conformal_spatial_metric_der_raised1[d][i][j];
              conformal_extrinsic_curvature_flux[d][i][j] -= conformal_lapse * (0.5 * evolution_func) * conformal_aux_vect_raised[d] * conformal_spatial_metric[i][j];

              if (i == d) {
                conformal_extrinsic_curvature_flux[d][i][j] += 0.5 * conformal_lapse * conformal_lapse_der[j];
                conformal_extrinsic_curvature_flux[d][i][j] += conformal_lapse * conformal_aux_vect[j];
                for (int r = 0; r < 3; r++) {
                  conformal_extrinsic_curvature_flux[d][i][j] -= 0.5 * conformal_lapse * conformal_spatial_metric_der_raised3[j][r][r];
                }
              }

              if (j == d) {
                conformal_extrinsic_curvature_flux[d][i][j] += 0.5 * conformal_lapse * conformal_lapse_der[i];
                conformal_extrinsic_curvature_flux[d][i][j] += conformal_lapse * conformal_aux_vect[i];
                for (int r = 0; r < 3; r++) {
                  conformal_extrinsic_curvature_flux[d][i][j] -= 0.5 * conformal_lapse * conformal_spatial_metric_der_raised3[i][r][r];
                }
              }
            }
          }
        }

        double conformal_spatial_metric_der_flux[3][3][3][3];
        for (int d = 0; d < 3; d++) {
          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
              for (int k = 0; k < 3; k++) {
                conformal_spatial_metric_der_flux[d][k][i][j] = 0.0;
              }

              for (int r = 0; r < 3; r++) {
                conformal_spatial_metric_der_flux[d][d][i][j] -= conformal_shift[r] * conformal_spatial_metric_der[r][i][j];
                conformal_spatial_metric_der_flux[d][d][i][j] -= 2.0 * conformal_shift[r] * conformal_fact_der[r] * conformal_spatial_metric[i][j];
              }

              conformal_spatial_metric_der_flux[d][d][i][j] += conformal_lapse * (conformal_extrinsic_curvature[i][j] /
                (conformal_fact * conformal_fact * conformal_fact * conformal_fact));
              conformal_spatial_metric_der_flux[d][d][i][j] -= conformal_lapse * symmetrized_conformal_shift[i][j];
            }
          }
        }

        double conformal_lapse_der_flux[3][3];
        for (int d = 0; d < 3; d++) {
          for (int i = 0; i < 3; i++) {
            conformal_lapse_der_flux[d][i] = 0.0;
          }
        }
        
        for (int d = 0; d < 3; d++) {
          for (int r = 0; r < 3; r++) {
            conformal_lapse_der_flux[d][d] -= conformal_shift[r] * conformal_lapse_der[r];
          }
          conformal_lapse_der_flux[d][d] += conformal_lapse * slicing_func;
        }

        double conformal_aux_vect_flux[3][3];
        for (int d = 0; d < 3; d++) {
          for (int i = 0; i < 3; i++) {
            conformal_aux_vect_flux[d][i] = -conformal_shift[d] * conformal_aux_vect[i];
            conformal_aux_vect_flux[d][i] -= 4.0 * conformal_shift[d] * conformal_fact_der[i];
            conformal_aux_vect_flux[d][i] += conformal_shift_der_switched[d][i];
            conformal_aux_vect_flux[d][i] -= conformal_shift_der[i][d];
          }
        }

        double fluxes[3][42] = {
          { conformal_extrinsic_curvature_flux[0][0][0], conformal_extrinsic_curvature_flux[0][0][1], conformal_extrinsic_curvature_flux[0][0][2],
            conformal_extrinsic_curvature_flux[0][1][0], conformal_extrinsic_curvature_flux[0][1][1], conformal_extrinsic_curvature_flux[0][1][2],
            conformal_extrinsic_curvature_flux[0][2][0], conformal_extrinsic_curvature_flux[0][2][1], conformal_extrinsic_curvature_flux[0][2][2],
            conformal_spatial_metric_der_flux[0][0][0][0], conformal_spatial_metric_der_flux[0][0][0][1], conformal_spatial_metric_der_flux[0][0][0][2],
            conformal_spatial_metric_der_flux[0][0][1][0], conformal_spatial_metric_der_flux[0][0][1][1], conformal_spatial_metric_der_flux[0][0][1][2],
            conformal_spatial_metric_der_flux[0][0][2][0], conformal_spatial_metric_der_flux[0][0][2][1], conformal_spatial_metric_der_flux[0][0][2][2],
            conformal_spatial_metric_der_flux[0][1][0][0], conformal_spatial_metric_der_flux[0][1][0][1], conformal_spatial_metric_der_flux[0][1][0][2],
            conformal_spatial_metric_der_flux[0][1][1][0], conformal_spatial_metric_der_flux[0][1][1][1], conformal_spatial_metric_der_flux[0][1][1][2],
            conformal_spatial_metric_der_flux[0][1][2][0], conformal_spatial_metric_der_flux[0][1][2][1], conformal_spatial_metric_der_flux[0][1][2][2],
            conformal_spatial_metric_der_flux[0][2][0][0], conformal_spatial_metric_der_flux[0][2][0][1], conformal_spatial_metric_der_flux[0][2][0][2],
            conformal_spatial_metric_der_flux[0][2][1][0], conformal_spatial_metric_der_flux[0][2][1][1], conformal_spatial_metric_der_flux[0][2][1][2],
            conformal_spatial_metric_der_flux[0][2][2][0], conformal_spatial_metric_der_flux[0][2][2][1], conformal_spatial_metric_der_flux[0][2][2][2],
            conformal_lapse_der_flux[0][0], conformal_lapse_der_flux[0][1], conformal_lapse_der_flux[0][2],
            conformal_aux_vect_flux[0][0], conformal_aux_vect_flux[0][1], conformal_aux_vect_flux[0][2] },
          { conformal_extrinsic_curvature_flux[1][0][0], conformal_extrinsic_curvature_flux[1][0][1], conformal_extrinsic_curvature_flux[1][0][2],
            conformal_extrinsic_curvature_flux[1][1][0], conformal_extrinsic_curvature_flux[1][1][1], conformal_extrinsic_curvature_flux[1][1][2],
            conformal_extrinsic_curvature_flux[1][2][0], conformal_extrinsic_curvature_flux[1][2][1], conformal_extrinsic_curvature_flux[1][2][2],
            conformal_spatial_metric_der_flux[1][0][0][0], conformal_spatial_metric_der_flux[1][0][0][1], conformal_spatial_metric_der_flux[1][0][0][2],
            conformal_spatial_metric_der_flux[1][0][1][0], conformal_spatial_metric_der_flux[1][0][1][1], conformal_spatial_metric_der_flux[1][0][1][2],
            conformal_spatial_metric_der_flux[1][0][2][0], conformal_spatial_metric_der_flux[1][0][2][1], conformal_spatial_metric_der_flux[1][0][2][2],
            conformal_spatial_metric_der_flux[1][1][0][0], conformal_spatial_metric_der_flux[1][1][0][1], conformal_spatial_metric_der_flux[1][1][0][2],
            conformal_spatial_metric_der_flux[1][1][1][0], conformal_spatial_metric_der_flux[1][1][1][1], conformal_spatial_metric_der_flux[1][1][1][2],
            conformal_spatial_metric_der_flux[1][1][2][0], conformal_spatial_metric_der_flux[1][1][2][1], conformal_spatial_metric_der_flux[1][1][2][2],
            conformal_spatial_metric_der_flux[1][2][0][0], conformal_spatial_metric_der_flux[1][2][0][1], conformal_spatial_metric_der_flux[1][2][0][2],
            conformal_spatial_metric_der_flux[1][2][1][0], conformal_spatial_metric_der_flux[1][2][1][1], conformal_spatial_metric_der_flux[1][2][1][2],
            conformal_spatial_metric_der_flux[1][2][2][0], conformal_spatial_metric_der_flux[1][2][2][1], conformal_spatial_metric_der_flux[1][2][2][2],
            conformal_lapse_der_flux[1][0], conformal_lapse_der_flux[1][1], conformal_lapse_der_flux[1][2],
            conformal_aux_vect_flux[1][0], conformal_aux_vect_flux[1][1], conformal_aux_vect_flux[1][2] },
          { conformal_extrinsic_curvature_flux[2][0][0], conformal_extrinsic_curvature_flux[2][0][1], conformal_extrinsic_curvature_flux[2][0][2],
            conformal_extrinsic_curvature_flux[2][1][0], conformal_extrinsic_curvature_flux[2][1][1], conformal_extrinsic_curvature_flux[2][1][2],
            conformal_extrinsic_curvature_flux[2][2][0], conformal_extrinsic_curvature_flux[2][2][1], conformal_extrinsic_curvature_flux[2][2][2],
            conformal_spatial_metric_der_flux[2][0][0][0], conformal_spatial_metric_der_flux[2][0][0][1], conformal_spatial_metric_der_flux[2][0][0][2],
            conformal_spatial_metric_der_flux[2][0][1][0], conformal_spatial_metric_der_flux[2][0][1][1], conformal_spatial_metric_der_flux[2][0][1][2],
            conformal_spatial_metric_der_flux[2][0][2][0], conformal_spatial_metric_der_flux[2][0][2][1], conformal_spatial_metric_der_flux[2][0][2][2],
            conformal_spatial_metric_der_flux[2][1][0][0], conformal_spatial_metric_der_flux[2][1][0][1], conformal_spatial_metric_der_flux[2][1][0][2],
            conformal_spatial_metric_der_flux[2][1][1][0], conformal_spatial_metric_der_flux[2][1][1][1], conformal_spatial_metric_der_flux[2][1][1][2],
            conformal_spatial_metric_der_flux[2][1][2][0], conformal_spatial_metric_der_flux[2][1][2][1], conformal_spatial_metric_der_flux[2][1][2][2],
            conformal_spatial_metric_der_flux[2][2][0][0], conformal_spatial_metric_der_flux[2][2][0][1], conformal_spatial_metric_der_flux[2][2][0][2],
            conformal_spatial_metric_der_flux[2][2][1][0], conformal_spatial_metric_der_flux[2][2][1][1], conformal_spatial_metric_der_flux[2][2][1][2],
            conformal_spatial_metric_der_flux[2][2][2][0], conformal_spatial_metric_der_flux[2][2][2][1], conformal_spatial_metric_der_flux[2][2][2][2],
            conformal_lapse_der_flux[2][0], conformal_lapse_der_flux[2][1], conformal_lapse_der_flux[2][2],
            conformal_aux_vect_flux[2][0], conformal_aux_vect_flux[2][1], conformal_aux_vect_flux[2][2] },
        };

        double norm[3][3] = {
          { 1.0, 0.0, 0.0 },
          { 0.0, 1.0, 0.0 },
          { 0.0, 0.0, 1.0 },
        };

        double tau1[3][3] = {
          { 0.0, 1.0, 0.0 },
          { 1.0, 0.0, 0.0 },
          { 1.0, 0.0, 0.0 },
        };

        double tau2[3][3] = {
          { 0.0, 0.0, 1.0 },
          { 0.0, 0.0, -1.0 },
          { 0.0, 1.0, 0.0 },
        };

        double q_local[77], flux_local[77], flux[77];
        for (int d = 0; d < 3; d++) {
          vacuum_einstein_conformal->rotate_to_local_func(vacuum_einstein_conformal, tau1[d], tau2[d], norm[d], q, q_local);
          gkyl_vacuum_einstein_conformal_flux(excision_threshold, spacetime_slicing, spacetime_evolution, q_local, flux_local);
          vacuum_einstein_conformal->rotate_to_global_func(vacuum_einstein_conformal, tau1[d], tau2[d], norm[d], flux_local, flux);

          for (int i = 0; i < 42; i++) {
            TEST_CHECK( gkyl_compare(flux[i + 10], fluxes[d][i], 1e-6) );
          }
        }
        
        double q_l[77], q_g[77];
        for (int d = 0; d < 3; d++) {
          gkyl_wv_eqn_rotate_to_local(vacuum_einstein_conformal, tau1[d], tau2[d], norm[d], q, q_l);
          gkyl_wv_eqn_rotate_to_global(vacuum_einstein_conformal, tau1[d], tau2[d], norm[d], q_l, q_g);

          for (int i = 0; i < 77; i++) {
            TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
          }

          double w1[77], q1[77];
          vacuum_einstein_conformal->cons_to_riem(vacuum_einstein_conformal, q_local, q_local, w1);
          vacuum_einstein_conformal->riem_to_cons(vacuum_einstein_conformal, q_local, w1, q1);

          for (int i = 0; i < 77; i++) {
            TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(conformal_spatial_metric[i]);
        gkyl_free(inv_conformal_spatial_metric[i]);
        gkyl_free(conformal_extrinsic_curvature[i]);
        gkyl_free(conformal_shift_der[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(conformal_spatial_metric_der[i][j]);
        }
        gkyl_free(conformal_spatial_metric_der[i]);
        gkyl_free(bssn_conformal_fact_der2[i]);
      }
      gkyl_free(conformal_spatial_metric);
      gkyl_free(inv_conformal_spatial_metric);
      gkyl_free(conformal_extrinsic_curvature);
      gkyl_free(conformal_shift);
      gkyl_free(conformal_lapse_der);
      gkyl_free(conformal_shift_der);
      gkyl_free(conformal_spatial_metric_der);
      gkyl_free(conformal_fact_der);
      gkyl_free(bssn_conformal_fact_der);
      gkyl_free(bssn_conformal_fact_der2);
    }
  }

  gkyl_wv_eqn_release(vacuum_einstein_conformal);
  gkyl_gr_spacetime_release(spacetime);
}

TEST_LIST = {
  { "vacuum_einstein_conformal_basic_minkowski", test_vacuum_einstein_conformal_basic_minkowski },
  { "vacuum_einstein_conformal_basic_schwarzschild", test_vacuum_einstein_conformal_basic_schwarzschild },
  { NULL, NULL },
};