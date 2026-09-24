#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_vacuum_einstein.h>
#include <gkyl_wv_vacuum_einstein_priv.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_gr_blackhole.h>

void
test_vacuum_einstein_basic_minkowski()
{
  double excision_threshold = 0.3;
  enum gkyl_spacetime_slicing spacetime_slicing = GKYL_HARMONIC_SLICING;
  enum gkyl_spacetime_evolution spacetime_evolution = GKYL_EINSTEIN_EVOLUTION;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  struct gkyl_wv_eqn *vacuum_einstein = gkyl_wv_vacuum_einstein_new(excision_threshold, spacetime_slicing, spacetime_evolution, false);

  TEST_CHECK( vacuum_einstein->num_equations == 64 );
  TEST_CHECK( vacuum_einstein->num_waves == 2 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double spatial_det, lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
        inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
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

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
      
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
      spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

      spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
      spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            spatial_metric_der[i][j][k] = 0.5 * spatial_metric_der[i][j][k];
          }

          shift_der[i][j] = 0.5 * shift_der[i][j];
        }
      }

      for (int i = 0; i < 3; i++) {
        lapse_der[i] = lapse_der[i] / lapse;
      }

      double extrinsic_curvature_trace = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          extrinsic_curvature_trace += inv_spatial_metric[i][j] * extrinsic_curvature[i][j];
        }
      }

      double spatial_metric_der_raised1[3][3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            spatial_metric_der_raised1[k][i][j] = 0.0;
            
            for (int l = 0; l < 3; l++) {
              spatial_metric_der_raised1[k][i][j] += inv_spatial_metric[k][l] * spatial_metric_der[l][i][j];
            }
          }
        }
      }
    
      double spatial_metric_der_raised3[3][3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            spatial_metric_der_raised3[i][j][k] = 0.0;

            for (int l = 0; l < 3; l++) {
              spatial_metric_der_raised3[i][j][k] += inv_spatial_metric[l][k] * spatial_metric_der[i][j][l];
            }
          }
        }
      }

      double aux_vect[3];
      for (int i = 0; i < 3; i++) {
        aux_vect[i] = 0.0;

        for (int s = 0; s < 3; s++) {
          aux_vect[i] += spatial_metric_der_raised3[i][s][s];
          aux_vect[i] -= spatial_metric_der_raised1[s][s][i];
        }
      }

      double aux_vect_raised[3];
      for (int k = 0; k < 3; k++) {
        aux_vect_raised[k] = 0.0;
            
        for (int l = 0; l < 3; l++) {
          aux_vect_raised[k] += inv_spatial_metric[k][l] * aux_vect[l];
        }
      }

      double shift_der_lowered[3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          shift_der_lowered[i][j] = 0.0;

          for (int k = 0; k < 3; k++) {
            shift_der_lowered[i][j] += spatial_metric[k][j] * shift_der[i][k];
          }
        }
      }

      double shift_der_switched[3][3];
      for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
          shift_der_switched[i][k] = 0.0;

          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              shift_der_switched[i][k] += inv_spatial_metric[i][l] * spatial_metric[m][k] * shift_der[l][m];
            }
          }
        }
      }

      double symmetrized_shift[3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          symmetrized_shift[i][j] = (1.0 / lapse) * (shift_der_lowered[i][j] + shift_der_lowered[j][i]);
        }
      }

      double q[64];
      q[0] = spatial_metric[0][0]; q[1] = spatial_metric[0][1]; q[2] = spatial_metric[0][2];
      q[3] = spatial_metric[1][0]; q[4] = spatial_metric[1][1]; q[5] = spatial_metric[1][2];
      q[6] = spatial_metric[2][0]; q[7] = spatial_metric[2][1]; q[8] = spatial_metric[2][2];

      q[9] = lapse;

      q[10] = extrinsic_curvature[0][0]; q[11] = extrinsic_curvature[0][1]; q[12] = extrinsic_curvature[0][2];
      q[13] = extrinsic_curvature[1][0]; q[14] = extrinsic_curvature[1][1]; q[15] = extrinsic_curvature[1][2];
      q[16] = extrinsic_curvature[2][0]; q[17] = extrinsic_curvature[2][1]; q[18] = extrinsic_curvature[2][2];

      q[19] = spatial_metric_der[0][0][0]; q[20] = spatial_metric_der[0][0][1]; q[21] = spatial_metric_der[0][0][2];
      q[22] = spatial_metric_der[0][1][0]; q[23] = spatial_metric_der[0][1][1]; q[24] = spatial_metric_der[0][1][2];
      q[25] = spatial_metric_der[0][2][0]; q[26] = spatial_metric_der[0][2][1]; q[27] = spatial_metric_der[0][2][2];

      q[28] = spatial_metric_der[1][0][0]; q[29] = spatial_metric_der[1][0][1]; q[30] = spatial_metric_der[1][0][2];
      q[31] = spatial_metric_der[1][1][0]; q[32] = spatial_metric_der[1][1][1]; q[33] = spatial_metric_der[1][1][2];
      q[34] = spatial_metric_der[1][2][0]; q[35] = spatial_metric_der[1][2][1]; q[36] = spatial_metric_der[1][2][2];

      q[37] = spatial_metric_der[2][0][0]; q[38] = spatial_metric_der[2][0][1]; q[39] = spatial_metric_der[2][0][2];
      q[40] = spatial_metric_der[2][1][0]; q[41] = spatial_metric_der[2][1][1]; q[42] = spatial_metric_der[2][1][2];
      q[43] = spatial_metric_der[2][2][0]; q[44] = spatial_metric_der[2][2][1]; q[45] = spatial_metric_der[2][2][2];

      q[46] = lapse_der[0]; q[47] = lapse_der[1]; q[48] = lapse_der[2];

      q[49] = aux_vect[0]; q[50] = aux_vect[1]; q[51] = aux_vect[2];

      q[52] = shift[0]; q[53] = shift[1]; q[54] = shift[2];
      
      q[55] = shift_der[0][0]; q[56] = shift_der[0][1]; q[57] = shift_der[0][2];
      q[58] = shift_der[1][0]; q[59] = shift_der[1][1]; q[60] = shift_der[1][2];
      q[61] = shift_der[2][0]; q[62] = shift_der[2][1]; q[63] = shift_der[2][2];

      double evolution_func = 1.0;
      double slicing_func = extrinsic_curvature_trace;

      double extrinsic_curvature_flux[3][3][3];
      for (int d = 0; d < 3; d++) {
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            extrinsic_curvature_flux[d][i][j] = -shift[d] * extrinsic_curvature[i][j];
            extrinsic_curvature_flux[d][i][j] += lapse * spatial_metric_der_raised1[d][i][j];
            extrinsic_curvature_flux[d][i][j] -= lapse * (0.5 * evolution_func) * aux_vect_raised[d] * spatial_metric[i][j];

            if (i == d) {
              extrinsic_curvature_flux[d][i][j] += 0.5 * lapse * lapse_der[j];
              extrinsic_curvature_flux[d][i][j] += lapse * aux_vect[j];
              for (int r = 0; r < 3; r++) {
                extrinsic_curvature_flux[d][i][j] -= 0.5 * lapse * spatial_metric_der_raised3[j][r][r];
              }
            }

            if (j == d) {
              extrinsic_curvature_flux[d][i][j] += 0.5 * lapse * lapse_der[i];
              extrinsic_curvature_flux[d][i][j] += lapse * aux_vect[i];
              for (int r = 0; r < 3; r++) {
                extrinsic_curvature_flux[d][i][j] -= 0.5 * lapse * spatial_metric_der_raised3[i][r][r];
              }
            }
          }
        }
      }
      
      double spatial_metric_der_flux[3][3][3][3];
      for (int d = 0; d < 3; d++) {
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              spatial_metric_der_flux[d][k][i][j] = 0.0;
            }

            for (int r = 0; r < 3; r++) {
              spatial_metric_der_flux[d][d][i][j] -= shift[r] * spatial_metric_der[r][i][j];
            }

            spatial_metric_der_flux[d][d][i][j] += lapse * (extrinsic_curvature[i][j] - symmetrized_shift[i][j]);
          }
        }
      }

      double lapse_der_flux[3][3];
      for (int d = 0; d < 3; d++) {
        for (int i = 0; i < 3; i++) {
          lapse_der_flux[d][i] = 0.0;
        }
      }
      
      for (int d = 0; d < 3; d++) {
        for (int r = 0; r < 3; r++) {
          lapse_der_flux[d][d] -= shift[r] * lapse_der[r];
        }
        lapse_der_flux[d][d] += lapse * slicing_func;
      }

      double aux_vect_flux[3][3];
      for (int d = 0; d < 3; d++) {
        for (int i = 0; i < 3; i++) {
          aux_vect_flux[d][i] = -shift[d] * aux_vect[i];
          aux_vect_flux[d][i] += shift_der_switched[d][i];
          aux_vect_flux[d][i] -= shift_der[i][d];
        }
      }

      double fluxes[3][42] = {
        { extrinsic_curvature_flux[0][0][0], extrinsic_curvature_flux[0][0][1], extrinsic_curvature_flux[0][0][2],
          extrinsic_curvature_flux[0][1][0], extrinsic_curvature_flux[0][1][1], extrinsic_curvature_flux[0][1][2],
          extrinsic_curvature_flux[0][2][0], extrinsic_curvature_flux[0][2][1], extrinsic_curvature_flux[0][2][2],
          spatial_metric_der_flux[0][0][0][0], spatial_metric_der_flux[0][0][0][1], spatial_metric_der_flux[0][0][0][2],
          spatial_metric_der_flux[0][0][1][0], spatial_metric_der_flux[0][0][1][1], spatial_metric_der_flux[0][0][1][2],
          spatial_metric_der_flux[0][0][2][0], spatial_metric_der_flux[0][0][2][1], spatial_metric_der_flux[0][0][2][2],
          spatial_metric_der_flux[0][1][0][0], spatial_metric_der_flux[0][1][0][1], spatial_metric_der_flux[0][1][0][2],
          spatial_metric_der_flux[0][1][1][0], spatial_metric_der_flux[0][1][1][1], spatial_metric_der_flux[0][1][1][2],
          spatial_metric_der_flux[0][1][2][0], spatial_metric_der_flux[0][1][2][1], spatial_metric_der_flux[0][1][2][2],
          spatial_metric_der_flux[0][2][0][0], spatial_metric_der_flux[0][2][0][1], spatial_metric_der_flux[0][2][0][2],
          spatial_metric_der_flux[0][2][1][0], spatial_metric_der_flux[0][2][1][1], spatial_metric_der_flux[0][2][1][2],
          spatial_metric_der_flux[0][2][2][0], spatial_metric_der_flux[0][2][2][1], spatial_metric_der_flux[0][2][2][2],
          lapse_der_flux[0][0], lapse_der_flux[0][1], lapse_der_flux[0][2],
          aux_vect_flux[0][0], aux_vect_flux[0][1], aux_vect_flux[0][2] },
        { extrinsic_curvature_flux[1][0][0], extrinsic_curvature_flux[1][0][1], extrinsic_curvature_flux[1][0][2],
          extrinsic_curvature_flux[1][1][0], extrinsic_curvature_flux[1][1][1], extrinsic_curvature_flux[1][1][2],
          extrinsic_curvature_flux[1][2][0], extrinsic_curvature_flux[1][2][1], extrinsic_curvature_flux[1][2][2],
          spatial_metric_der_flux[1][0][0][0], spatial_metric_der_flux[1][0][0][1], spatial_metric_der_flux[1][0][0][2],
          spatial_metric_der_flux[1][0][1][0], spatial_metric_der_flux[1][0][1][1], spatial_metric_der_flux[1][0][1][2],
          spatial_metric_der_flux[1][0][2][0], spatial_metric_der_flux[1][0][2][1], spatial_metric_der_flux[1][0][2][2],
          spatial_metric_der_flux[1][1][0][0], spatial_metric_der_flux[1][1][0][1], spatial_metric_der_flux[1][1][0][2],
          spatial_metric_der_flux[1][1][1][0], spatial_metric_der_flux[1][1][1][1], spatial_metric_der_flux[1][1][1][2],
          spatial_metric_der_flux[1][1][2][0], spatial_metric_der_flux[1][1][2][1], spatial_metric_der_flux[1][1][2][2],
          spatial_metric_der_flux[1][2][0][0], spatial_metric_der_flux[1][2][0][1], spatial_metric_der_flux[1][2][0][2],
          spatial_metric_der_flux[1][2][1][0], spatial_metric_der_flux[1][2][1][1], spatial_metric_der_flux[1][2][1][2],
          spatial_metric_der_flux[1][2][2][0], spatial_metric_der_flux[1][2][2][1], spatial_metric_der_flux[1][2][2][2],
          lapse_der_flux[1][0], lapse_der_flux[1][1], lapse_der_flux[1][2],
          aux_vect_flux[1][0], aux_vect_flux[1][1], aux_vect_flux[1][2] },
        { extrinsic_curvature_flux[2][0][0], extrinsic_curvature_flux[2][0][1], extrinsic_curvature_flux[2][0][2],
          extrinsic_curvature_flux[2][1][0], extrinsic_curvature_flux[2][1][1], extrinsic_curvature_flux[2][1][2],
          extrinsic_curvature_flux[2][2][0], extrinsic_curvature_flux[2][2][1], extrinsic_curvature_flux[2][2][2],
          spatial_metric_der_flux[2][0][0][0], spatial_metric_der_flux[2][0][0][1], spatial_metric_der_flux[2][0][0][2],
          spatial_metric_der_flux[2][0][1][0], spatial_metric_der_flux[2][0][1][1], spatial_metric_der_flux[2][0][1][2],
          spatial_metric_der_flux[2][0][2][0], spatial_metric_der_flux[2][0][2][1], spatial_metric_der_flux[2][0][2][2],
          spatial_metric_der_flux[2][1][0][0], spatial_metric_der_flux[2][1][0][1], spatial_metric_der_flux[2][1][0][2],
          spatial_metric_der_flux[2][1][1][0], spatial_metric_der_flux[2][1][1][1], spatial_metric_der_flux[2][1][1][2],
          spatial_metric_der_flux[2][1][2][0], spatial_metric_der_flux[2][1][2][1], spatial_metric_der_flux[2][1][2][2],
          spatial_metric_der_flux[2][2][0][0], spatial_metric_der_flux[2][2][0][1], spatial_metric_der_flux[2][2][0][2],
          spatial_metric_der_flux[2][2][1][0], spatial_metric_der_flux[2][2][1][1], spatial_metric_der_flux[2][2][1][2],
          spatial_metric_der_flux[2][2][2][0], spatial_metric_der_flux[2][2][2][1], spatial_metric_der_flux[2][2][2][2],
          lapse_der_flux[2][0], lapse_der_flux[2][1], lapse_der_flux[2][2],
          aux_vect_flux[2][0], aux_vect_flux[2][1], aux_vect_flux[2][2] },
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

      double q_local[64], flux_local[64], flux[64];
      for (int d = 0; d < 3; d++) {
        vacuum_einstein->rotate_to_local_func(vacuum_einstein, tau1[d], tau2[d], norm[d], q, q_local);
        gkyl_vacuum_einstein_flux(excision_threshold, spacetime_slicing, spacetime_evolution, q_local, flux_local);
        vacuum_einstein->rotate_to_global_func(vacuum_einstein, tau1[d], tau2[d], norm[d], flux_local, flux);

        for (int i = 0; i < 42; i++) {
          TEST_CHECK( gkyl_compare(flux[i + 10], fluxes[d][i], 1e-8) );
        }
      }
      
      double q_l[64], q_g[64];
      for (int d = 0; d < 3; d++) {
        gkyl_wv_eqn_rotate_to_local(vacuum_einstein, tau1[d], tau2[d], norm[d], q, q_l);
        gkyl_wv_eqn_rotate_to_global(vacuum_einstein, tau1[d], tau2[d], norm[d], q_l, q_g);

        for (int i = 0; i < 64; i++) {
          TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
        }

        double w1[64], q1[64];
        vacuum_einstein->cons_to_riem(vacuum_einstein, q_local, q_local, w1);
        vacuum_einstein->riem_to_cons(vacuum_einstein, q_local, w1, q1);

        for (int i = 0; i < 64; i++) {
          TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(inv_spatial_metric[i]);
        gkyl_free(extrinsic_curvature[i]);
        gkyl_free(shift_der[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der[i][j]);
        }
        gkyl_free(spatial_metric_der[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(inv_spatial_metric);
      gkyl_free(extrinsic_curvature);
      gkyl_free(shift);
      gkyl_free(lapse_der);
      gkyl_free(shift_der);
      gkyl_free(spatial_metric_der);
    }
  }

  gkyl_wv_eqn_release(vacuum_einstein);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_vacuum_einstein_basic_schwarzschild()
{
  double excision_threshold = 0.3;
  enum gkyl_spacetime_slicing spacetime_slicing = GKYL_1PLUSLOG_SLICING;
  enum gkyl_spacetime_evolution spacetime_evolution = GKYL_EINSTEIN_EVOLUTION;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *vacuum_einstein = gkyl_wv_vacuum_einstein_new(excision_threshold, spacetime_slicing, spacetime_evolution, false);

  TEST_CHECK( vacuum_einstein->num_equations == 64 );
  TEST_CHECK( vacuum_einstein->num_waves == 2 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double spatial_det, lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
        inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
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

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
      
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
      spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

      spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
      spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            spatial_metric_der[i][j][k] = 0.5 * spatial_metric_der[i][j][k];
          }

          shift_der[i][j] = 0.5 * shift_der[i][j];
        }
      }

      for (int i = 0; i < 3; i++) {
        lapse_der[i] = lapse_der[i] / lapse;
      }

      double extrinsic_curvature_trace = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          extrinsic_curvature_trace += inv_spatial_metric[i][j] * extrinsic_curvature[i][j];
        }
      }

      double spatial_metric_der_raised1[3][3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            spatial_metric_der_raised1[k][i][j] = 0.0;
            
            for (int l = 0; l < 3; l++) {
              spatial_metric_der_raised1[k][i][j] += inv_spatial_metric[k][l] * spatial_metric_der[l][i][j];
            }
          }
        }
      }
    
      double spatial_metric_der_raised3[3][3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            spatial_metric_der_raised3[i][j][k] = 0.0;

            for (int l = 0; l < 3; l++) {
              spatial_metric_der_raised3[i][j][k] += inv_spatial_metric[l][k] * spatial_metric_der[i][j][l];
            }
          }
        }
      }

      double aux_vect[3];
      for (int i = 0; i < 3; i++) {
        aux_vect[i] = 0.0;

        for (int s = 0; s < 3; s++) {
          aux_vect[i] += spatial_metric_der_raised3[i][s][s];
          aux_vect[i] -= spatial_metric_der_raised1[s][s][i];
        }
      }

      double aux_vect_raised[3];
      for (int k = 0; k < 3; k++) {
        aux_vect_raised[k] = 0.0;
            
        for (int l = 0; l < 3; l++) {
          aux_vect_raised[k] += inv_spatial_metric[k][l] * aux_vect[l];
        }
      }

      double shift_der_lowered[3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          shift_der_lowered[i][j] = 0.0;

          for (int k = 0; k < 3; k++) {
            shift_der_lowered[i][j] += spatial_metric[k][j] * shift_der[i][k];
          }
        }
      }

      double shift_der_switched[3][3];
      for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
          shift_der_switched[i][k] = 0.0;

          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              shift_der_switched[i][k] += inv_spatial_metric[i][l] * spatial_metric[m][k] * shift_der[l][m];
            }
          }
        }
      }

      double symmetrized_shift[3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          symmetrized_shift[i][j] = (1.0 / lapse) * (shift_der_lowered[i][j] + shift_der_lowered[j][i]);
        }
      }

      if (!in_excision_region) {
        double q[64];
        q[0] = spatial_metric[0][0]; q[1] = spatial_metric[0][1]; q[2] = spatial_metric[0][2];
        q[3] = spatial_metric[1][0]; q[4] = spatial_metric[1][1]; q[5] = spatial_metric[1][2];
        q[6] = spatial_metric[2][0]; q[7] = spatial_metric[2][1]; q[8] = spatial_metric[2][2];

        q[9] = lapse;

        q[10] = extrinsic_curvature[0][0]; q[11] = extrinsic_curvature[0][1]; q[12] = extrinsic_curvature[0][2];
        q[13] = extrinsic_curvature[1][0]; q[14] = extrinsic_curvature[1][1]; q[15] = extrinsic_curvature[1][2];
        q[16] = extrinsic_curvature[2][0]; q[17] = extrinsic_curvature[2][1]; q[18] = extrinsic_curvature[2][2];

        q[19] = spatial_metric_der[0][0][0]; q[20] = spatial_metric_der[0][0][1]; q[21] = spatial_metric_der[0][0][2];
        q[22] = spatial_metric_der[0][1][0]; q[23] = spatial_metric_der[0][1][1]; q[24] = spatial_metric_der[0][1][2];
        q[25] = spatial_metric_der[0][2][0]; q[26] = spatial_metric_der[0][2][1]; q[27] = spatial_metric_der[0][2][2];

        q[28] = spatial_metric_der[1][0][0]; q[29] = spatial_metric_der[1][0][1]; q[30] = spatial_metric_der[1][0][2];
        q[31] = spatial_metric_der[1][1][0]; q[32] = spatial_metric_der[1][1][1]; q[33] = spatial_metric_der[1][1][2];
        q[34] = spatial_metric_der[1][2][0]; q[35] = spatial_metric_der[1][2][1]; q[36] = spatial_metric_der[1][2][2];

        q[37] = spatial_metric_der[2][0][0]; q[38] = spatial_metric_der[2][0][1]; q[39] = spatial_metric_der[2][0][2];
        q[40] = spatial_metric_der[2][1][0]; q[41] = spatial_metric_der[2][1][1]; q[42] = spatial_metric_der[2][1][2];
        q[43] = spatial_metric_der[2][2][0]; q[44] = spatial_metric_der[2][2][1]; q[45] = spatial_metric_der[2][2][2];

        q[46] = lapse_der[0]; q[47] = lapse_der[1]; q[48] = lapse_der[2];

        q[49] = aux_vect[0]; q[50] = aux_vect[1]; q[51] = aux_vect[2];

        q[52] = shift[0]; q[53] = shift[1]; q[54] = shift[2];
        
        q[55] = shift_der[0][0]; q[56] = shift_der[0][1]; q[57] = shift_der[0][2];
        q[58] = shift_der[1][0]; q[59] = shift_der[1][1]; q[60] = shift_der[1][2];
        q[61] = shift_der[2][0]; q[62] = shift_der[2][1]; q[63] = shift_der[2][2];

        double evolution_func = 1.0;
        double slicing_func = 2.0 * extrinsic_curvature_trace / lapse;

        double extrinsic_curvature_flux[3][3][3];
        for (int d = 0; d < 3; d++) {
          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
              extrinsic_curvature_flux[d][i][j] = -shift[d] * extrinsic_curvature[i][j];
              extrinsic_curvature_flux[d][i][j] += lapse * spatial_metric_der_raised1[d][i][j];
              extrinsic_curvature_flux[d][i][j] -= lapse * (0.5 * evolution_func) * aux_vect_raised[d] * spatial_metric[i][j];

              if (i == d) {
                extrinsic_curvature_flux[d][i][j] += 0.5 * lapse * lapse_der[j];
                extrinsic_curvature_flux[d][i][j] += lapse * aux_vect[j];
                for (int r = 0; r < 3; r++) {
                  extrinsic_curvature_flux[d][i][j] -= 0.5 * lapse * spatial_metric_der_raised3[j][r][r];
                }
              }

              if (j == d) {
                extrinsic_curvature_flux[d][i][j] += 0.5 * lapse * lapse_der[i];
                extrinsic_curvature_flux[d][i][j] += lapse * aux_vect[i];
                for (int r = 0; r < 3; r++) {
                  extrinsic_curvature_flux[d][i][j] -= 0.5 * lapse * spatial_metric_der_raised3[i][r][r];
                }
              }
            }
          }
        }

        double spatial_metric_der_flux[3][3][3][3];
        for (int d = 0; d < 3; d++) {
          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
              for (int k = 0; k < 3; k++) {
                spatial_metric_der_flux[d][k][i][j] = 0.0;
              }

              for (int r = 0; r < 3; r++) {
                spatial_metric_der_flux[d][d][i][j] -= shift[r] * spatial_metric_der[r][i][j];
              }

              spatial_metric_der_flux[d][d][i][j] += lapse * (extrinsic_curvature[i][j] - symmetrized_shift[i][j]);
            }
          }
        }

        double lapse_der_flux[3][3];
        for (int d = 0; d < 3; d++) {
          for (int i = 0; i < 3; i++) {
            lapse_der_flux[d][i] = 0.0;
          }
        }
        
        for (int d = 0; d < 3; d++) {
          for (int r = 0; r < 3; r++) {
            lapse_der_flux[d][d] -= shift[r] * lapse_der[r];
          }
          lapse_der_flux[d][d] += lapse * slicing_func;
        }

        double aux_vect_flux[3][3];
        for (int d = 0; d < 3; d++) {
          for (int i = 0; i < 3; i++) {
            aux_vect_flux[d][i] = -shift[d] * aux_vect[i];
            aux_vect_flux[d][i] += shift_der_switched[d][i];
            aux_vect_flux[d][i] -= shift_der[i][d];
          }
        }

        double fluxes[3][42] = {
          { extrinsic_curvature_flux[0][0][0], extrinsic_curvature_flux[0][0][1], extrinsic_curvature_flux[0][0][2],
            extrinsic_curvature_flux[0][1][0], extrinsic_curvature_flux[0][1][1], extrinsic_curvature_flux[0][1][2],
            extrinsic_curvature_flux[0][2][0], extrinsic_curvature_flux[0][2][1], extrinsic_curvature_flux[0][2][2],
            spatial_metric_der_flux[0][0][0][0], spatial_metric_der_flux[0][0][0][1], spatial_metric_der_flux[0][0][0][2],
            spatial_metric_der_flux[0][0][1][0], spatial_metric_der_flux[0][0][1][1], spatial_metric_der_flux[0][0][1][2],
            spatial_metric_der_flux[0][0][2][0], spatial_metric_der_flux[0][0][2][1], spatial_metric_der_flux[0][0][2][2],
            spatial_metric_der_flux[0][1][0][0], spatial_metric_der_flux[0][1][0][1], spatial_metric_der_flux[0][1][0][2],
            spatial_metric_der_flux[0][1][1][0], spatial_metric_der_flux[0][1][1][1], spatial_metric_der_flux[0][1][1][2],
            spatial_metric_der_flux[0][1][2][0], spatial_metric_der_flux[0][1][2][1], spatial_metric_der_flux[0][1][2][2],
            spatial_metric_der_flux[0][2][0][0], spatial_metric_der_flux[0][2][0][1], spatial_metric_der_flux[0][2][0][2],
            spatial_metric_der_flux[0][2][1][0], spatial_metric_der_flux[0][2][1][1], spatial_metric_der_flux[0][2][1][2],
            spatial_metric_der_flux[0][2][2][0], spatial_metric_der_flux[0][2][2][1], spatial_metric_der_flux[0][2][2][2],
            lapse_der_flux[0][0], lapse_der_flux[0][1], lapse_der_flux[0][2],
            aux_vect_flux[0][0], aux_vect_flux[0][1], aux_vect_flux[0][2] },
          { extrinsic_curvature_flux[1][0][0], extrinsic_curvature_flux[1][0][1], extrinsic_curvature_flux[1][0][2],
            extrinsic_curvature_flux[1][1][0], extrinsic_curvature_flux[1][1][1], extrinsic_curvature_flux[1][1][2],
            extrinsic_curvature_flux[1][2][0], extrinsic_curvature_flux[1][2][1], extrinsic_curvature_flux[1][2][2],
            spatial_metric_der_flux[1][0][0][0], spatial_metric_der_flux[1][0][0][1], spatial_metric_der_flux[1][0][0][2],
            spatial_metric_der_flux[1][0][1][0], spatial_metric_der_flux[1][0][1][1], spatial_metric_der_flux[1][0][1][2],
            spatial_metric_der_flux[1][0][2][0], spatial_metric_der_flux[1][0][2][1], spatial_metric_der_flux[1][0][2][2],
            spatial_metric_der_flux[1][1][0][0], spatial_metric_der_flux[1][1][0][1], spatial_metric_der_flux[1][1][0][2],
            spatial_metric_der_flux[1][1][1][0], spatial_metric_der_flux[1][1][1][1], spatial_metric_der_flux[1][1][1][2],
            spatial_metric_der_flux[1][1][2][0], spatial_metric_der_flux[1][1][2][1], spatial_metric_der_flux[1][1][2][2],
            spatial_metric_der_flux[1][2][0][0], spatial_metric_der_flux[1][2][0][1], spatial_metric_der_flux[1][2][0][2],
            spatial_metric_der_flux[1][2][1][0], spatial_metric_der_flux[1][2][1][1], spatial_metric_der_flux[1][2][1][2],
            spatial_metric_der_flux[1][2][2][0], spatial_metric_der_flux[1][2][2][1], spatial_metric_der_flux[1][2][2][2],
            lapse_der_flux[1][0], lapse_der_flux[1][1], lapse_der_flux[1][2],
            aux_vect_flux[1][0], aux_vect_flux[1][1], aux_vect_flux[1][2] },
          { extrinsic_curvature_flux[2][0][0], extrinsic_curvature_flux[2][0][1], extrinsic_curvature_flux[2][0][2],
            extrinsic_curvature_flux[2][1][0], extrinsic_curvature_flux[2][1][1], extrinsic_curvature_flux[2][1][2],
            extrinsic_curvature_flux[2][2][0], extrinsic_curvature_flux[2][2][1], extrinsic_curvature_flux[2][2][2],
            spatial_metric_der_flux[2][0][0][0], spatial_metric_der_flux[2][0][0][1], spatial_metric_der_flux[2][0][0][2],
            spatial_metric_der_flux[2][0][1][0], spatial_metric_der_flux[2][0][1][1], spatial_metric_der_flux[2][0][1][2],
            spatial_metric_der_flux[2][0][2][0], spatial_metric_der_flux[2][0][2][1], spatial_metric_der_flux[2][0][2][2],
            spatial_metric_der_flux[2][1][0][0], spatial_metric_der_flux[2][1][0][1], spatial_metric_der_flux[2][1][0][2],
            spatial_metric_der_flux[2][1][1][0], spatial_metric_der_flux[2][1][1][1], spatial_metric_der_flux[2][1][1][2],
            spatial_metric_der_flux[2][1][2][0], spatial_metric_der_flux[2][1][2][1], spatial_metric_der_flux[2][1][2][2],
            spatial_metric_der_flux[2][2][0][0], spatial_metric_der_flux[2][2][0][1], spatial_metric_der_flux[2][2][0][2],
            spatial_metric_der_flux[2][2][1][0], spatial_metric_der_flux[2][2][1][1], spatial_metric_der_flux[2][2][1][2],
            spatial_metric_der_flux[2][2][2][0], spatial_metric_der_flux[2][2][2][1], spatial_metric_der_flux[2][2][2][2],
            lapse_der_flux[2][0], lapse_der_flux[2][1], lapse_der_flux[2][2],
            aux_vect_flux[2][0], aux_vect_flux[2][1], aux_vect_flux[2][2] },
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

        double q_local[64], flux_local[64], flux[64];
        for (int d = 0; d < 3; d++) {
          vacuum_einstein->rotate_to_local_func(vacuum_einstein, tau1[d], tau2[d], norm[d], q, q_local);
          gkyl_vacuum_einstein_flux(excision_threshold, spacetime_slicing, spacetime_evolution, q_local, flux_local);
          vacuum_einstein->rotate_to_global_func(vacuum_einstein, tau1[d], tau2[d], norm[d], flux_local, flux);

          for (int i = 0; i < 42; i++) {
            TEST_CHECK( gkyl_compare(flux[i + 10], fluxes[d][i], 1e-6) );
          }
        }
        
        double q_l[64], q_g[64];
        for (int d = 0; d < 3; d++) {
          gkyl_wv_eqn_rotate_to_local(vacuum_einstein, tau1[d], tau2[d], norm[d], q, q_l);
          gkyl_wv_eqn_rotate_to_global(vacuum_einstein, tau1[d], tau2[d], norm[d], q_l, q_g);

          for (int i = 0; i < 64; i++) {
            TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
          }

          double w1[64], q1[64];
          vacuum_einstein->cons_to_riem(vacuum_einstein, q_local, q_local, w1);
          vacuum_einstein->riem_to_cons(vacuum_einstein, q_local, w1, q1);

          for (int i = 0; i < 64; i++) {
            TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(inv_spatial_metric[i]);
        gkyl_free(extrinsic_curvature[i]);
        gkyl_free(shift_der[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der[i][j]);
        }
        gkyl_free(spatial_metric_der[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(inv_spatial_metric);
      gkyl_free(extrinsic_curvature);
      gkyl_free(shift);
      gkyl_free(lapse_der);
      gkyl_free(shift_der);
      gkyl_free(spatial_metric_der);
    }
  }

  gkyl_wv_eqn_release(vacuum_einstein);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_vacuum_einstein_waves_schwarzschild()
{
  double excision_threshold = 0.3;
  enum gkyl_spacetime_slicing spacetime_slicing = GKYL_1PLUSLOG_SLICING;
  enum gkyl_spacetime_evolution spacetime_evolution = GKYL_EINSTEIN_EVOLUTION;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *vacuum_einstein = gkyl_wv_vacuum_einstein_new(excision_threshold, spacetime_slicing, spacetime_evolution, false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double spatial_det_l, spatial_det_r;
      double lapse_l, lapse_r;
      double *shift_l = gkyl_malloc(sizeof(double[3]));
      double *shift_r = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region_l, in_excision_region_r;

      double **spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
      double **spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
      double **inv_spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
      double **inv_spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
        spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
        inv_spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
        inv_spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **extrinsic_curvature_l = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        extrinsic_curvature_l[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der_l = gkyl_malloc(sizeof(double[3]));
      double *lapse_der_r = gkyl_malloc(sizeof(double[3]));
      double **shift_der_l = gkyl_malloc(sizeof(double*[3]));
      double **shift_der_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der_l[i] = gkyl_malloc(sizeof(double[3]));
        shift_der_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der_l = gkyl_malloc(sizeof(double**[3]));
      double ***spatial_metric_der_r = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der_l[i] = gkyl_malloc(sizeof(double*[3]));
        spatial_metric_der_r[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der_l[i][j] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_der_r[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_det_l);
      spacetime->spatial_metric_det_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_det_r);
      spacetime->lapse_function_func(spacetime, 0.0, x - 0.1, y, 0.0, &lapse_l);
      spacetime->lapse_function_func(spacetime, 0.0, x + 0.1, y, 0.0, &lapse_r);
      spacetime->shift_vector_func(spacetime, 0.0, x - 0.1, y, 0.0, &shift_l);
      spacetime->shift_vector_func(spacetime, 0.0, x + 0.1, y, 0.0, &shift_r);
      spacetime->excision_region_func(spacetime, 0.0, x - 0.1, y, 0.0, &in_excision_region_l);
      spacetime->excision_region_func(spacetime, 0.0, x + 0.1, y, 0.0, &in_excision_region_r);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_metric_l);
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_metric_r);
      spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, &inv_spatial_metric_l);
      spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, &inv_spatial_metric_r);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_l);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_r);

      spacetime->lapse_function_der_func(spacetime, 0.0, x - 0.1, y, 0.0, 0.1, 0.1, 0.1, &lapse_der_l);
      spacetime->lapse_function_der_func(spacetime, 0.0, x + 0.1, y, 0.0, 0.1, 0.1, 0.1, &lapse_der_r);
      spacetime->shift_vector_der_func(spacetime, 0.0, x - 0.1, y, 0.0, 0.1, 0.1, 0.1, &shift_der_l);
      spacetime->shift_vector_der_func(spacetime, 0.0, x + 0.1, y, 0.0, 0.1, 0.1, 0.1, &shift_der_r);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x - 0.1, y, 0.0, 0.1, 0.1, 0.1, &spatial_metric_der_l);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x + 0.1, y, 0.0, 0.1, 0.1, 0.1, &spatial_metric_der_r);

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            spatial_metric_der_l[i][j][k] = 0.5 * spatial_metric_der_l[i][j][k];
            spatial_metric_der_r[i][j][k] = 0.5 * spatial_metric_der_r[i][j][k];
          }

          shift_der_l[i][j] = 0.5 * shift_der_l[i][j];
          shift_der_r[i][j] = 0.5 * shift_der_r[i][j];
        }
      }

      for (int i = 0; i < 3; i++) {
        lapse_der_l[i] = lapse_der_l[i] / lapse_l;
        lapse_der_r[i] = lapse_der_r[i] / lapse_r;
      }

      double extrinsic_curvature_trace_l = 0.0;
      double extrinsic_curvature_trace_r = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          extrinsic_curvature_trace_l += inv_spatial_metric_l[i][j] * extrinsic_curvature_l[i][j];
          extrinsic_curvature_trace_r += inv_spatial_metric_r[i][j] * extrinsic_curvature_r[i][j];
        }
      }

      double spatial_metric_der_raised1_l[3][3][3];
      double spatial_metric_der_raised1_r[3][3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            spatial_metric_der_raised1_l[k][i][j] = 0.0;
            spatial_metric_der_raised1_r[k][i][j] = 0.0;
            
            for (int l = 0; l < 3; l++) {
              spatial_metric_der_raised1_l[k][i][j] += inv_spatial_metric_l[k][l] * spatial_metric_der_l[l][i][j];
              spatial_metric_der_raised1_r[k][i][j] += inv_spatial_metric_r[k][l] * spatial_metric_der_r[l][i][j];
            }
          }
        }
      }
    
      double spatial_metric_der_raised3_l[3][3][3];
      double spatial_metric_der_raised3_r[3][3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            spatial_metric_der_raised3_l[i][j][k] = 0.0;
            spatial_metric_der_raised3_r[i][j][k] = 0.0;

            for (int l = 0; l < 3; l++) {
              spatial_metric_der_raised3_l[i][j][k] += inv_spatial_metric_l[l][k] * spatial_metric_der_l[i][j][l];
              spatial_metric_der_raised3_r[i][j][k] += inv_spatial_metric_r[l][k] * spatial_metric_der_r[i][j][l];
            }
          }
        }
      }

      double aux_vect_l[3];
      double aux_vect_r[3];
      for (int i = 0; i < 3; i++) {
        aux_vect_l[i] = 0.0;
        aux_vect_r[i] = 0.0;

        for (int s = 0; s < 3; s++) {
          aux_vect_l[i] += spatial_metric_der_raised3_l[i][s][s];
          aux_vect_l[i] -= spatial_metric_der_raised1_l[s][s][i];

          aux_vect_r[i] += spatial_metric_der_raised3_r[i][s][s];
          aux_vect_r[i] -= spatial_metric_der_raised1_r[s][s][i];
        }
      }

      double aux_vect_raised_l[3];
      double aux_vect_raised_r[3];
      for (int k = 0; k < 3; k++) {
        aux_vect_raised_l[k] = 0.0;
        aux_vect_raised_r[k] = 0.0;

        for (int l = 0; l < 3; l++) {
          aux_vect_raised_l[k] += inv_spatial_metric_l[k][l] * aux_vect_l[l];
          aux_vect_raised_r[k] += inv_spatial_metric_r[k][l] * aux_vect_r[l];
        }
      }

      double shift_der_lowered_l[3][3];
      double shift_der_lowered_r[3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          shift_der_lowered_l[i][j] = 0.0;
          shift_der_lowered_r[i][j] = 0.0;

          for (int k = 0; k < 3; k++) {
            shift_der_lowered_l[i][j] += spatial_metric_l[k][j] * shift_der_l[i][k];
            shift_der_lowered_r[i][j] += spatial_metric_r[k][j] * shift_der_r[i][k];
          }
        }
      }

      double shift_der_switched_l[3][3];
      double shift_der_switched_r[3][3];
      for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
          shift_der_switched_l[i][k] = 0.0;
          shift_der_switched_r[i][k] = 0.0;

          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              shift_der_switched_l[i][k] += inv_spatial_metric_l[i][l] * spatial_metric_l[m][k] * shift_der_l[l][m];
              shift_der_switched_r[i][k] += inv_spatial_metric_r[i][l] * spatial_metric_r[m][k] * shift_der_r[l][m];
            }
          }
        }
      }

      double symmetrized_shift_l[3][3];
      double symmetrized_shift_r[3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          symmetrized_shift_l[i][j] = (1.0 / lapse_l) * (shift_der_lowered_l[i][j] + shift_der_lowered_l[j][i]);
          symmetrized_shift_r[i][j] = (1.0 / lapse_r) * (shift_der_lowered_r[i][j] + shift_der_lowered_r[j][i]);
        }
      }

      if (!in_excision_region_l && !in_excision_region_r) {
        double ql[64], qr[64];
        ql[0] = spatial_metric_l[0][0]; ql[1] = spatial_metric_l[0][1]; ql[2] = spatial_metric_l[0][2];
        ql[3] = spatial_metric_l[1][0]; ql[4] = spatial_metric_l[1][1]; ql[5] = spatial_metric_l[1][2];
        ql[6] = spatial_metric_l[2][0]; ql[7] = spatial_metric_l[2][1]; ql[8] = spatial_metric_l[2][2];

        ql[9] = lapse_l;

        ql[10] = extrinsic_curvature_l[0][0]; ql[11] = extrinsic_curvature_l[0][1]; ql[12] = extrinsic_curvature_l[0][2];
        ql[13] = extrinsic_curvature_l[1][0]; ql[14] = extrinsic_curvature_l[1][1]; ql[15] = extrinsic_curvature_l[1][2];
        ql[16] = extrinsic_curvature_l[2][0]; ql[17] = extrinsic_curvature_l[2][1]; ql[18] = extrinsic_curvature_l[2][2];

        ql[19] = spatial_metric_der_l[0][0][0]; ql[20] = spatial_metric_der_l[0][0][1]; ql[21] = spatial_metric_der_l[0][0][2];
        ql[22] = spatial_metric_der_l[0][1][0]; ql[23] = spatial_metric_der_l[0][1][1]; ql[24] = spatial_metric_der_l[0][1][2];
        ql[25] = spatial_metric_der_l[0][2][0]; ql[26] = spatial_metric_der_l[0][2][1]; ql[27] = spatial_metric_der_l[0][2][2];

        ql[28] = spatial_metric_der_l[1][0][0]; ql[29] = spatial_metric_der_l[1][0][1]; ql[30] = spatial_metric_der_l[1][0][2];
        ql[31] = spatial_metric_der_l[1][1][0]; ql[32] = spatial_metric_der_l[1][1][1]; ql[33] = spatial_metric_der_l[1][1][2];
        ql[34] = spatial_metric_der_l[1][2][0]; ql[35] = spatial_metric_der_l[1][2][1]; ql[36] = spatial_metric_der_l[1][2][2];

        ql[37] = spatial_metric_der_l[2][0][0]; ql[38] = spatial_metric_der_l[2][0][1]; ql[39] = spatial_metric_der_l[2][0][2];
        ql[40] = spatial_metric_der_l[2][1][0]; ql[41] = spatial_metric_der_l[2][1][1]; ql[42] = spatial_metric_der_l[2][1][2];
        ql[43] = spatial_metric_der_l[2][2][0]; ql[44] = spatial_metric_der_l[2][2][1]; ql[45] = spatial_metric_der_l[2][2][2];

        ql[46] = lapse_der_l[0]; ql[47] = lapse_der_l[1]; ql[48] = lapse_der_l[2];

        ql[49] = aux_vect_l[0]; ql[50] = aux_vect_l[1]; ql[51] = aux_vect_l[2];

        ql[52] = shift_l[0]; ql[53] = shift_l[1]; ql[54] = shift_l[2];
        
        ql[55] = shift_der_l[0][0]; ql[56] = shift_der_l[0][1]; ql[57] = shift_der_l[0][2];
        ql[58] = shift_der_l[1][0]; ql[59] = shift_der_l[1][1]; ql[60] = shift_der_l[1][2];
        ql[61] = shift_der_l[2][0]; ql[62] = shift_der_l[2][1]; ql[63] = shift_der_l[2][2];

        qr[0] = spatial_metric_r[0][0]; qr[1] = spatial_metric_r[0][1]; qr[2] = spatial_metric_r[0][2];
        qr[3] = spatial_metric_r[1][0]; qr[4] = spatial_metric_r[1][1]; qr[5] = spatial_metric_r[1][2];
        qr[6] = spatial_metric_r[2][0]; qr[7] = spatial_metric_r[2][1]; qr[8] = spatial_metric_r[2][2];

        qr[9] = lapse_r;

        qr[10] = extrinsic_curvature_r[0][0]; qr[11] = extrinsic_curvature_r[0][1]; qr[12] = extrinsic_curvature_r[0][2];
        qr[13] = extrinsic_curvature_r[1][0]; qr[14] = extrinsic_curvature_r[1][1]; qr[15] = extrinsic_curvature_r[1][2];
        qr[16] = extrinsic_curvature_r[2][0]; qr[17] = extrinsic_curvature_r[2][1]; qr[18] = extrinsic_curvature_r[2][2];

        qr[19] = spatial_metric_der_r[0][0][0]; qr[20] = spatial_metric_der_r[0][0][1]; qr[21] = spatial_metric_der_r[0][0][2];
        qr[22] = spatial_metric_der_r[0][1][0]; qr[23] = spatial_metric_der_r[0][1][1]; qr[24] = spatial_metric_der_r[0][1][2];
        qr[25] = spatial_metric_der_r[0][2][0]; qr[26] = spatial_metric_der_r[0][2][1]; qr[27] = spatial_metric_der_r[0][2][2];

        qr[28] = spatial_metric_der_r[1][0][0]; qr[29] = spatial_metric_der_r[1][0][1]; qr[30] = spatial_metric_der_r[1][0][2];
        qr[31] = spatial_metric_der_r[1][1][0]; qr[32] = spatial_metric_der_r[1][1][1]; qr[33] = spatial_metric_der_r[1][1][2];
        qr[34] = spatial_metric_der_r[1][2][0]; qr[35] = spatial_metric_der_r[1][2][1]; qr[36] = spatial_metric_der_r[1][2][2];

        qr[37] = spatial_metric_der_r[2][0][0]; qr[38] = spatial_metric_der_r[2][0][1]; qr[39] = spatial_metric_der_r[2][0][2];
        qr[40] = spatial_metric_der_r[2][1][0]; qr[41] = spatial_metric_der_r[2][1][1]; qr[42] = spatial_metric_der_r[2][1][2];
        qr[43] = spatial_metric_der_r[2][2][0]; qr[44] = spatial_metric_der_r[2][2][1]; qr[45] = spatial_metric_der_r[2][2][2];

        qr[46] = lapse_der_r[0]; qr[47] = lapse_der_r[1]; qr[48] = lapse_der_r[2];

        qr[49] = aux_vect_r[0]; qr[50] = aux_vect_r[1]; qr[51] = aux_vect_r[2];

        qr[52] = shift_r[0]; qr[53] = shift_r[1]; qr[54] = shift_r[2];
        
        qr[55] = shift_der_r[0][0]; qr[56] = shift_der_r[0][1]; qr[57] = shift_der_r[0][2];
        qr[58] = shift_der_r[1][0]; qr[59] = shift_der_r[1][1]; qr[60] = shift_der_r[1][2];
        qr[61] = shift_der_r[2][0]; qr[62] = shift_der_r[2][1]; qr[63] = shift_der_r[2][2];

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

        for (int d = 0; d < 3; d++) {
          double speeds[2], waves[2 * 64], waves_local[2 * 64];

          double ql_local[64], qr_local[64];
          gkyl_wv_eqn_rotate_to_local(vacuum_einstein, tau1[d], tau2[d], norm[d], ql, ql_local);
          gkyl_wv_eqn_rotate_to_local(vacuum_einstein, tau1[d], tau2[d], norm[d], qr, qr_local);

          double delta[64];
          for (int i = 0; i < 64; i++) {
            delta[i] = qr_local[i] - ql_local[i];
          }

          gkyl_wv_eqn_waves(vacuum_einstein, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

          double apdq_local[64], amdq_local[64];
          gkyl_wv_eqn_qfluct(vacuum_einstein, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

          for (int i = 0; i < 2; i++) {
            gkyl_wv_eqn_rotate_to_global(vacuum_einstein, tau1[d], tau2[d], norm[d], &waves_local[i * 64], &waves[i * 64]);
          }

          double apdq[64], amdq[64];
          gkyl_wv_eqn_rotate_to_global(vacuum_einstein, tau1[d], tau2[d], norm[d], apdq_local, apdq);
          gkyl_wv_eqn_rotate_to_global(vacuum_einstein, tau1[d], tau2[d], norm[d], amdq_local, amdq);

          double fl_local[64], fr_local[64];
          gkyl_vacuum_einstein_flux(excision_threshold, spacetime_slicing, spacetime_evolution, ql_local, fl_local);
          gkyl_vacuum_einstein_flux(excision_threshold, spacetime_slicing, spacetime_evolution, qr_local, fr_local);

          double fl[64], fr[64];
          gkyl_wv_eqn_rotate_to_global(vacuum_einstein, tau1[d], tau2[d], norm[d], fl_local, fl);
          gkyl_wv_eqn_rotate_to_global(vacuum_einstein, tau1[d], tau2[d], norm[d], fr_local, fr);

          for (int i = 0; i < 64; i++) {
            TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-11) );
          }
        }
      }
     
      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric_l[i]);
        gkyl_free(spatial_metric_r[i]);
        gkyl_free(inv_spatial_metric_l[i]);
        gkyl_free(inv_spatial_metric_r[i]);
        gkyl_free(extrinsic_curvature_l[i]);
        gkyl_free(extrinsic_curvature_r[i]);
        gkyl_free(shift_der_l[i]);
        gkyl_free(shift_der_r[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der_l[i][j]);
          gkyl_free(spatial_metric_der_r[i][j]);
        }
        gkyl_free(spatial_metric_der_l[i]);
        gkyl_free(spatial_metric_der_r[i]);
      }
      gkyl_free(spatial_metric_l);
      gkyl_free(spatial_metric_r);
      gkyl_free(inv_spatial_metric_l);
      gkyl_free(inv_spatial_metric_r);
      gkyl_free(extrinsic_curvature_l);
      gkyl_free(extrinsic_curvature_r);
      gkyl_free(shift_l);
      gkyl_free(shift_r);
      gkyl_free(lapse_der_l);
      gkyl_free(lapse_der_r);
      gkyl_free(shift_der_l);
      gkyl_free(shift_der_r);
      gkyl_free(spatial_metric_der_l);
      gkyl_free(spatial_metric_der_r);
    }
  }
}

void
test_vacuum_einstein_waves_kerr()
{
  double excision_threshold = 0.3;
  enum gkyl_spacetime_slicing spacetime_slicing = GKYL_1PLUSLOG_SLICING;
  enum gkyl_spacetime_evolution spacetime_evolution = GKYL_EINSTEIN_EVOLUTION;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *vacuum_einstein = gkyl_wv_vacuum_einstein_new(excision_threshold, spacetime_slicing, spacetime_evolution, false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double spatial_det_l, spatial_det_r;
      double lapse_l, lapse_r;
      double *shift_l = gkyl_malloc(sizeof(double[3]));
      double *shift_r = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region_l, in_excision_region_r;

      double **spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
      double **spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
      double **inv_spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
      double **inv_spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
        spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
        inv_spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
        inv_spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **extrinsic_curvature_l = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        extrinsic_curvature_l[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der_l = gkyl_malloc(sizeof(double[3]));
      double *lapse_der_r = gkyl_malloc(sizeof(double[3]));
      double **shift_der_l = gkyl_malloc(sizeof(double*[3]));
      double **shift_der_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der_l[i] = gkyl_malloc(sizeof(double[3]));
        shift_der_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der_l = gkyl_malloc(sizeof(double**[3]));
      double ***spatial_metric_der_r = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der_l[i] = gkyl_malloc(sizeof(double*[3]));
        spatial_metric_der_r[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der_l[i][j] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_der_r[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_det_l);
      spacetime->spatial_metric_det_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_det_r);
      spacetime->lapse_function_func(spacetime, 0.0, x - 0.1, y, 0.0, &lapse_l);
      spacetime->lapse_function_func(spacetime, 0.0, x + 0.1, y, 0.0, &lapse_r);
      spacetime->shift_vector_func(spacetime, 0.0, x - 0.1, y, 0.0, &shift_l);
      spacetime->shift_vector_func(spacetime, 0.0, x + 0.1, y, 0.0, &shift_r);
      spacetime->excision_region_func(spacetime, 0.0, x - 0.1, y, 0.0, &in_excision_region_l);
      spacetime->excision_region_func(spacetime, 0.0, x + 0.1, y, 0.0, &in_excision_region_r);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_metric_l);
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_metric_r);
      spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, &inv_spatial_metric_l);
      spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, &inv_spatial_metric_r);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_l);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_r);

      spacetime->lapse_function_der_func(spacetime, 0.0, x - 0.1, y, 0.0, 0.1, 0.1, 0.1, &lapse_der_l);
      spacetime->lapse_function_der_func(spacetime, 0.0, x + 0.1, y, 0.0, 0.1, 0.1, 0.1, &lapse_der_r);
      spacetime->shift_vector_der_func(spacetime, 0.0, x - 0.1, y, 0.0, 0.1, 0.1, 0.1, &shift_der_l);
      spacetime->shift_vector_der_func(spacetime, 0.0, x + 0.1, y, 0.0, 0.1, 0.1, 0.1, &shift_der_r);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x - 0.1, y, 0.0, 0.1, 0.1, 0.1, &spatial_metric_der_l);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x + 0.1, y, 0.0, 0.1, 0.1, 0.1, &spatial_metric_der_r);

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            spatial_metric_der_l[i][j][k] = 0.5 * spatial_metric_der_l[i][j][k];
            spatial_metric_der_r[i][j][k] = 0.5 * spatial_metric_der_r[i][j][k];
          }

          shift_der_l[i][j] = 0.5 * shift_der_l[i][j];
          shift_der_r[i][j] = 0.5 * shift_der_r[i][j];
        }
      }

      for (int i = 0; i < 3; i++) {
        lapse_der_l[i] = lapse_der_l[i] / lapse_l;
        lapse_der_r[i] = lapse_der_r[i] / lapse_r;
      }

      double extrinsic_curvature_trace_l = 0.0;
      double extrinsic_curvature_trace_r = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          extrinsic_curvature_trace_l += inv_spatial_metric_l[i][j] * extrinsic_curvature_l[i][j];
          extrinsic_curvature_trace_r += inv_spatial_metric_r[i][j] * extrinsic_curvature_r[i][j];
        }
      }

      double spatial_metric_der_raised1_l[3][3][3];
      double spatial_metric_der_raised1_r[3][3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            spatial_metric_der_raised1_l[k][i][j] = 0.0;
            spatial_metric_der_raised1_r[k][i][j] = 0.0;
            
            for (int l = 0; l < 3; l++) {
              spatial_metric_der_raised1_l[k][i][j] += inv_spatial_metric_l[k][l] * spatial_metric_der_l[l][i][j];
              spatial_metric_der_raised1_r[k][i][j] += inv_spatial_metric_r[k][l] * spatial_metric_der_r[l][i][j];
            }
          }
        }
      }
    
      double spatial_metric_der_raised3_l[3][3][3];
      double spatial_metric_der_raised3_r[3][3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            spatial_metric_der_raised3_l[i][j][k] = 0.0;
            spatial_metric_der_raised3_r[i][j][k] = 0.0;

            for (int l = 0; l < 3; l++) {
              spatial_metric_der_raised3_l[i][j][k] += inv_spatial_metric_l[l][k] * spatial_metric_der_l[i][j][l];
              spatial_metric_der_raised3_r[i][j][k] += inv_spatial_metric_r[l][k] * spatial_metric_der_r[i][j][l];
            }
          }
        }
      }

      double aux_vect_l[3];
      double aux_vect_r[3];
      for (int i = 0; i < 3; i++) {
        aux_vect_l[i] = 0.0;
        aux_vect_r[i] = 0.0;

        for (int s = 0; s < 3; s++) {
          aux_vect_l[i] += spatial_metric_der_raised3_l[i][s][s];
          aux_vect_l[i] -= spatial_metric_der_raised1_l[s][s][i];

          aux_vect_r[i] += spatial_metric_der_raised3_r[i][s][s];
          aux_vect_r[i] -= spatial_metric_der_raised1_r[s][s][i];
        }
      }

      double aux_vect_raised_l[3];
      double aux_vect_raised_r[3];
      for (int k = 0; k < 3; k++) {
        aux_vect_raised_l[k] = 0.0;
        aux_vect_raised_r[k] = 0.0;

        for (int l = 0; l < 3; l++) {
          aux_vect_raised_l[k] += inv_spatial_metric_l[k][l] * aux_vect_l[l];
          aux_vect_raised_r[k] += inv_spatial_metric_r[k][l] * aux_vect_r[l];
        }
      }

      double shift_der_lowered_l[3][3];
      double shift_der_lowered_r[3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          shift_der_lowered_l[i][j] = 0.0;
          shift_der_lowered_r[i][j] = 0.0;

          for (int k = 0; k < 3; k++) {
            shift_der_lowered_l[i][j] += spatial_metric_l[k][j] * shift_der_l[i][k];
            shift_der_lowered_r[i][j] += spatial_metric_r[k][j] * shift_der_r[i][k];
          }
        }
      }

      double shift_der_switched_l[3][3];
      double shift_der_switched_r[3][3];
      for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
          shift_der_switched_l[i][k] = 0.0;
          shift_der_switched_r[i][k] = 0.0;

          for (int l = 0; l < 3; l++) {
            for (int m = 0; m < 3; m++) {
              shift_der_switched_l[i][k] += inv_spatial_metric_l[i][l] * spatial_metric_l[m][k] * shift_der_l[l][m];
              shift_der_switched_r[i][k] += inv_spatial_metric_r[i][l] * spatial_metric_r[m][k] * shift_der_r[l][m];
            }
          }
        }
      }

      double symmetrized_shift_l[3][3];
      double symmetrized_shift_r[3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          symmetrized_shift_l[i][j] = (1.0 / lapse_l) * (shift_der_lowered_l[i][j] + shift_der_lowered_l[j][i]);
          symmetrized_shift_r[i][j] = (1.0 / lapse_r) * (shift_der_lowered_r[i][j] + shift_der_lowered_r[j][i]);
        }
      }

      if (!in_excision_region_l && !in_excision_region_r) {
        double ql[64], qr[64];
        ql[0] = spatial_metric_l[0][0]; ql[1] = spatial_metric_l[0][1]; ql[2] = spatial_metric_l[0][2];
        ql[3] = spatial_metric_l[1][0]; ql[4] = spatial_metric_l[1][1]; ql[5] = spatial_metric_l[1][2];
        ql[6] = spatial_metric_l[2][0]; ql[7] = spatial_metric_l[2][1]; ql[8] = spatial_metric_l[2][2];

        ql[9] = lapse_l;

        ql[10] = extrinsic_curvature_l[0][0]; ql[11] = extrinsic_curvature_l[0][1]; ql[12] = extrinsic_curvature_l[0][2];
        ql[13] = extrinsic_curvature_l[1][0]; ql[14] = extrinsic_curvature_l[1][1]; ql[15] = extrinsic_curvature_l[1][2];
        ql[16] = extrinsic_curvature_l[2][0]; ql[17] = extrinsic_curvature_l[2][1]; ql[18] = extrinsic_curvature_l[2][2];

        ql[19] = spatial_metric_der_l[0][0][0]; ql[20] = spatial_metric_der_l[0][0][1]; ql[21] = spatial_metric_der_l[0][0][2];
        ql[22] = spatial_metric_der_l[0][1][0]; ql[23] = spatial_metric_der_l[0][1][1]; ql[24] = spatial_metric_der_l[0][1][2];
        ql[25] = spatial_metric_der_l[0][2][0]; ql[26] = spatial_metric_der_l[0][2][1]; ql[27] = spatial_metric_der_l[0][2][2];

        ql[28] = spatial_metric_der_l[1][0][0]; ql[29] = spatial_metric_der_l[1][0][1]; ql[30] = spatial_metric_der_l[1][0][2];
        ql[31] = spatial_metric_der_l[1][1][0]; ql[32] = spatial_metric_der_l[1][1][1]; ql[33] = spatial_metric_der_l[1][1][2];
        ql[34] = spatial_metric_der_l[1][2][0]; ql[35] = spatial_metric_der_l[1][2][1]; ql[36] = spatial_metric_der_l[1][2][2];

        ql[37] = spatial_metric_der_l[2][0][0]; ql[38] = spatial_metric_der_l[2][0][1]; ql[39] = spatial_metric_der_l[2][0][2];
        ql[40] = spatial_metric_der_l[2][1][0]; ql[41] = spatial_metric_der_l[2][1][1]; ql[42] = spatial_metric_der_l[2][1][2];
        ql[43] = spatial_metric_der_l[2][2][0]; ql[44] = spatial_metric_der_l[2][2][1]; ql[45] = spatial_metric_der_l[2][2][2];

        ql[46] = lapse_der_l[0]; ql[47] = lapse_der_l[1]; ql[48] = lapse_der_l[2];

        ql[49] = aux_vect_l[0]; ql[50] = aux_vect_l[1]; ql[51] = aux_vect_l[2];

        ql[52] = shift_l[0]; ql[53] = shift_l[1]; ql[54] = shift_l[2];
        
        ql[55] = shift_der_l[0][0]; ql[56] = shift_der_l[0][1]; ql[57] = shift_der_l[0][2];
        ql[58] = shift_der_l[1][0]; ql[59] = shift_der_l[1][1]; ql[60] = shift_der_l[1][2];
        ql[61] = shift_der_l[2][0]; ql[62] = shift_der_l[2][1]; ql[63] = shift_der_l[2][2];

        qr[0] = spatial_metric_r[0][0]; qr[1] = spatial_metric_r[0][1]; qr[2] = spatial_metric_r[0][2];
        qr[3] = spatial_metric_r[1][0]; qr[4] = spatial_metric_r[1][1]; qr[5] = spatial_metric_r[1][2];
        qr[6] = spatial_metric_r[2][0]; qr[7] = spatial_metric_r[2][1]; qr[8] = spatial_metric_r[2][2];

        qr[9] = lapse_r;

        qr[10] = extrinsic_curvature_r[0][0]; qr[11] = extrinsic_curvature_r[0][1]; qr[12] = extrinsic_curvature_r[0][2];
        qr[13] = extrinsic_curvature_r[1][0]; qr[14] = extrinsic_curvature_r[1][1]; qr[15] = extrinsic_curvature_r[1][2];
        qr[16] = extrinsic_curvature_r[2][0]; qr[17] = extrinsic_curvature_r[2][1]; qr[18] = extrinsic_curvature_r[2][2];

        qr[19] = spatial_metric_der_r[0][0][0]; qr[20] = spatial_metric_der_r[0][0][1]; qr[21] = spatial_metric_der_r[0][0][2];
        qr[22] = spatial_metric_der_r[0][1][0]; qr[23] = spatial_metric_der_r[0][1][1]; qr[24] = spatial_metric_der_r[0][1][2];
        qr[25] = spatial_metric_der_r[0][2][0]; qr[26] = spatial_metric_der_r[0][2][1]; qr[27] = spatial_metric_der_r[0][2][2];

        qr[28] = spatial_metric_der_r[1][0][0]; qr[29] = spatial_metric_der_r[1][0][1]; qr[30] = spatial_metric_der_r[1][0][2];
        qr[31] = spatial_metric_der_r[1][1][0]; qr[32] = spatial_metric_der_r[1][1][1]; qr[33] = spatial_metric_der_r[1][1][2];
        qr[34] = spatial_metric_der_r[1][2][0]; qr[35] = spatial_metric_der_r[1][2][1]; qr[36] = spatial_metric_der_r[1][2][2];

        qr[37] = spatial_metric_der_r[2][0][0]; qr[38] = spatial_metric_der_r[2][0][1]; qr[39] = spatial_metric_der_r[2][0][2];
        qr[40] = spatial_metric_der_r[2][1][0]; qr[41] = spatial_metric_der_r[2][1][1]; qr[42] = spatial_metric_der_r[2][1][2];
        qr[43] = spatial_metric_der_r[2][2][0]; qr[44] = spatial_metric_der_r[2][2][1]; qr[45] = spatial_metric_der_r[2][2][2];

        qr[46] = lapse_der_r[0]; qr[47] = lapse_der_r[1]; qr[48] = lapse_der_r[2];

        qr[49] = aux_vect_r[0]; qr[50] = aux_vect_r[1]; qr[51] = aux_vect_r[2];

        qr[52] = shift_r[0]; qr[53] = shift_r[1]; qr[54] = shift_r[2];
        
        qr[55] = shift_der_r[0][0]; qr[56] = shift_der_r[0][1]; qr[57] = shift_der_r[0][2];
        qr[58] = shift_der_r[1][0]; qr[59] = shift_der_r[1][1]; qr[60] = shift_der_r[1][2];
        qr[61] = shift_der_r[2][0]; qr[62] = shift_der_r[2][1]; qr[63] = shift_der_r[2][2];

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

        for (int d = 0; d < 3; d++) {
          double speeds[2], waves[2 * 64], waves_local[2 * 64];

          double ql_local[64], qr_local[64];
          gkyl_wv_eqn_rotate_to_local(vacuum_einstein, tau1[d], tau2[d], norm[d], ql, ql_local);
          gkyl_wv_eqn_rotate_to_local(vacuum_einstein, tau1[d], tau2[d], norm[d], qr, qr_local);

          double delta[64];
          for (int i = 0; i < 64; i++) {
            delta[i] = qr_local[i] - ql_local[i];
          }

          gkyl_wv_eqn_waves(vacuum_einstein, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

          double apdq_local[64], amdq_local[64];
          gkyl_wv_eqn_qfluct(vacuum_einstein, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

          for (int i = 0; i < 2; i++) {
            gkyl_wv_eqn_rotate_to_global(vacuum_einstein, tau1[d], tau2[d], norm[d], &waves_local[i * 64], &waves[i * 64]);
          }

          double apdq[64], amdq[64];
          gkyl_wv_eqn_rotate_to_global(vacuum_einstein, tau1[d], tau2[d], norm[d], apdq_local, apdq);
          gkyl_wv_eqn_rotate_to_global(vacuum_einstein, tau1[d], tau2[d], norm[d], amdq_local, amdq);

          double fl_local[64], fr_local[64];
          gkyl_vacuum_einstein_flux(excision_threshold, spacetime_slicing, spacetime_evolution, ql_local, fl_local);
          gkyl_vacuum_einstein_flux(excision_threshold, spacetime_slicing, spacetime_evolution, qr_local, fr_local);

          double fl[64], fr[64];
          gkyl_wv_eqn_rotate_to_global(vacuum_einstein, tau1[d], tau2[d], norm[d], fl_local, fl);
          gkyl_wv_eqn_rotate_to_global(vacuum_einstein, tau1[d], tau2[d], norm[d], fr_local, fr);

          for (int i = 0; i < 64; i++) {
            TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-12) );
          }
        }
      }
     
      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric_l[i]);
        gkyl_free(spatial_metric_r[i]);
        gkyl_free(inv_spatial_metric_l[i]);
        gkyl_free(inv_spatial_metric_r[i]);
        gkyl_free(extrinsic_curvature_l[i]);
        gkyl_free(extrinsic_curvature_r[i]);
        gkyl_free(shift_der_l[i]);
        gkyl_free(shift_der_r[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der_l[i][j]);
          gkyl_free(spatial_metric_der_r[i][j]);
        }
        gkyl_free(spatial_metric_der_l[i]);
        gkyl_free(spatial_metric_der_r[i]);
      }
      gkyl_free(spatial_metric_l);
      gkyl_free(spatial_metric_r);
      gkyl_free(inv_spatial_metric_l);
      gkyl_free(inv_spatial_metric_r);
      gkyl_free(extrinsic_curvature_l);
      gkyl_free(extrinsic_curvature_r);
      gkyl_free(shift_l);
      gkyl_free(shift_r);
      gkyl_free(lapse_der_l);
      gkyl_free(lapse_der_r);
      gkyl_free(shift_der_l);
      gkyl_free(shift_der_r);
      gkyl_free(spatial_metric_der_l);
      gkyl_free(spatial_metric_der_r);
    }
  }
}

TEST_LIST = {
  { "vacuum_einstein_basic_minkowski", test_vacuum_einstein_basic_minkowski },
  { "vacuum_einstein_basic_schwarzschild", test_vacuum_einstein_basic_schwarzschild },
  { "vacuum_einstein_waves_schwarzschild", test_vacuum_einstein_waves_schwarzschild },
  { "vacuum_einstein_waves_kerr", test_vacuum_einstein_waves_kerr },
  { NULL, NULL },
};