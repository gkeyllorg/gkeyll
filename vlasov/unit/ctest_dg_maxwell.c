#include <acutest.h>

#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_maxwell_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_position_map.h>

void
test_dg_max_ho()
{
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, 1);

  // Identity position map providing the (required) conf-space Jacobian.
  double lower[] = {0.0}, upper[] = {1.0};
  int cells[] = {8}, ghost[] = {1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);
  struct gkyl_range crange, crange_ext;
  gkyl_create_grid_ranges(&grid, ghost, &crange_ext, &crange);
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = {0};
  struct gkyl_vlasov_position_map *pos_map =
    gkyl_vlasov_position_map_new(&grid, &crange, &crange_ext, &basis, inp_pmap, false);

  // Input structure for building the dg eqn object
  struct gkyl_dg_maxwell_inp inp_dg_maxwell = {
    .cbasis = &basis,
    .crange = &crange,
    .jacob_pos = pos_map->jacob_pos,
    .conf_flux_surf = 0,
    .lightSpeed = 1.0,
    .field_id = GKYL_FIELD_E_B,
    .elcErrorSpeedFactor = 0.5,
    .mgnErrorSpeedFactor = 0.25,
    .use_gpu = false,
  };
  struct gkyl_dg_eqn *eqn = gkyl_dg_maxwell_inew(&inp_dg_maxwell);

  TEST_CHECK(eqn->num_equations == 8);

  // this is not possible from user code and should NOT be done. This
  // is for testing only
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);

  TEST_CHECK(maxwell->maxwell_data.c == 1.0);
  TEST_CHECK(maxwell->maxwell_data.chi == 0.5);
  TEST_CHECK(maxwell->maxwell_data.gamma == 0.25);

  struct dg_maxwell *m_on_dev = container_of(eqn->on_dev, struct dg_maxwell, eqn);

  TEST_CHECK(m_on_dev->maxwell_data.c == 1.0);
  TEST_CHECK(m_on_dev->maxwell_data.chi == 0.5);
  TEST_CHECK(m_on_dev->maxwell_data.gamma == 0.25);

  gkyl_dg_eqn_release(eqn);
  gkyl_vlasov_position_map_release(pos_map);
}

// Absorbing boundaries: with the central flux on the interior faces and the boundary surface
// kernels on every edge of the domain, the rate of change of the energy
//   (1/2) int [ |E|^2 + c^2 phi^2 + c^2 |B|^2 + psi^2 ] dx
// of any state is minus one half of the boundary integral of the wave speeds times the energy
// densities of the pairs they carry: c for the tangential fields, c chi for (E_n, phi) and
// c gamma for (B_n, psi). The ghost cells hold large values that must not enter.
static void
test_dg_maxwell_absorb(int cdim, int poly_order, bool use_tensor)
{
  struct gkyl_basis basis;
  if (use_tensor) {
    gkyl_cart_modal_tensor(&basis, cdim, poly_order);
  } else {
    gkyl_cart_modal_serendip(&basis, cdim, poly_order);
  }
  int num_basis = basis.num_basis;

  double lower[] = {0.0, -1.0, 0.5}, upper[] = {1.0, 2.0, 1.75};
  int cells[] = {4, 3, 2}, ghost[] = {1, 1, 1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
  struct gkyl_range crange, crange_ext;
  gkyl_create_grid_ranges(&grid, ghost, &crange_ext, &crange);
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = {0};
  struct gkyl_vlasov_position_map *pos_map =
    gkyl_vlasov_position_map_new(&grid, &crange, &crange_ext, &basis, inp_pmap, false);

  double light_speed = 1.5, elc_error_speed_fact = 0.75, mgn_error_speed_fact = 1.25;
  struct gkyl_dg_maxwell_inp inp_dg_maxwell = {
    .cbasis = &basis,
    .crange = &crange,
    .jacob_pos = pos_map->jacob_pos,
    .conf_flux_surf = 0,
    .lightSpeed = light_speed,
    .field_id = GKYL_FIELD_E_B,
    .elcErrorSpeedFactor = elc_error_speed_fact,
    .mgnErrorSpeedFactor = mgn_error_speed_fact,
    .use_gpu = false,
  };
  struct gkyl_dg_eqn *eqn = gkyl_dg_maxwell_inew(&inp_dg_maxwell);
  double chi = light_speed * elc_error_speed_fact, gamma = light_speed * mgn_error_speed_fact;
  double light_speed_sq = light_speed * light_speed;

  // Boundary surface kernels on every edge.
  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux_flags[2 * GKYL_MAX_DIM] = {0};
  for (int d = 0; d < cdim; ++d) {
    zero_flux_flags[d] = zero_flux_flags[d + cdim] = 1;
  }
  struct gkyl_hyper_dg *slvr =
    gkyl_hyper_dg_new(&grid, &basis, eqn, cdim, up_dirs, zero_flux_flags, 1, false);

  struct gkyl_array *em = gkyl_array_new(GKYL_DOUBLE, 8 * num_basis, crange_ext.volume);
  struct gkyl_array *rhs = gkyl_array_new(GKYL_DOUBLE, 8 * num_basis, crange_ext.volume);
  struct gkyl_array *cflrate = gkyl_array_new(GKYL_DOUBLE, 1, crange_ext.volume);

  // Large values everywhere (the ghost cells keep them), a different value in every coefficient
  // of every interior cell.
  gkyl_array_clear(em, 1.0e3);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &crange);
  while (gkyl_range_iter_next(&iter)) {
    long linc = gkyl_range_idx(&crange, iter.idx);
    double *em_d = gkyl_array_fetch(em, linc);
    for (int k = 0; k < 8 * num_basis; ++k) {
      em_d[k] = sin(1.0 + 0.7 * linc + 1.3 * k);
    }
  }

  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_hyper_dg_advance(slvr, &crange, em, cflrate, rhs);

  // Energy weights of the eight components.
  double weight[8] = {1.0, 1.0, 1.0, light_speed_sq, light_speed_sq, light_speed_sq, light_speed_sq,
                      1.0};

  double cell_fac = grid.cellVolume / pow(2.0, cdim);
  double energy_rate = 0.0;
  gkyl_range_iter_init(&iter, &crange);
  while (gkyl_range_iter_next(&iter)) {
    long linc = gkyl_range_idx(&crange, iter.idx);
    const double *em_d = gkyl_array_cfetch(em, linc);
    const double *rhs_d = gkyl_array_cfetch(rhs, linc);
    for (int c = 0; c < 8; ++c) {
      for (int k = 0; k < num_basis; ++k) {
        energy_rate += cell_fac * weight[c] * em_d[c * num_basis + k] * rhs_d[c * num_basis + k];
      }
    }
  }

  // Boundary integral with a four-point Gauss-Legendre rule in each direction of the face.
  double gl_node[4] = {
    -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526
  };
  double gl_weight[4] = {
    0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538
  };
  int num_face_quad = cdim == 1 ? 1 : (cdim == 2 ? 4 : 16);

  double boundary_rate = 0.0;
  for (int dir = 0; dir < cdim; ++dir) {
    double face_fac = grid.cellVolume / grid.dx[dir] / pow(2.0, cdim - 1);
    for (int edge = 0; edge < 2; ++edge) {
      gkyl_range_iter_init(&iter, &crange);
      while (gkyl_range_iter_next(&iter)) {
        if (iter.idx[dir] != (edge == 0 ? crange.lower[dir] : crange.upper[dir])) {
          continue;
        }
        const double *em_d = gkyl_array_cfetch(em, gkyl_range_idx(&crange, iter.idx));

        for (int n = 0; n < num_face_quad; ++n) {
          // Reference coordinates of the quadrature node on the face.
          double xi[3], quad_weight = 1.0;
          int node = n;
          for (int d = 0; d < cdim; ++d) {
            if (d == dir) {
              xi[d] = edge == 0 ? -1.0 : 1.0;
            } else {
              xi[d] = gl_node[node % 4];
              quad_weight *= gl_weight[node % 4];
              node /= 4;
            }
          }
          double val[8];
          for (int c = 0; c < 8; ++c) {
            val[c] = basis.eval_expand(xi, &em_d[c * num_basis]);
          }

          double density = 0.0;
          for (int c = 0; c < 3; ++c) {
            if (c == dir) {
              density += light_speed * chi * val[c] * val[c];
              density += light_speed * gamma * light_speed_sq * val[3 + c] * val[3 + c];
            } else {
              density += light_speed * val[c] * val[c];
              density += light_speed * light_speed_sq * val[3 + c] * val[3 + c];
            }
          }
          density += light_speed * chi * light_speed_sq * val[6] * val[6];
          density += light_speed * gamma * val[7] * val[7];

          boundary_rate -= 0.5 * face_fac * quad_weight * density;
        }
      }
    }
  }

  TEST_CHECK(boundary_rate < 0.0);
  TEST_CHECK(gkyl_compare(energy_rate, boundary_rate, 1e-12));
  TEST_MSG("energy rate %.15e, boundary integral %.15e", energy_rate, boundary_rate);

  gkyl_array_release(em);
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_hyper_dg_release(slvr);
  gkyl_dg_eqn_release(eqn);
  gkyl_vlasov_position_map_release(pos_map);
}

void
test_dg_maxwell_absorb_1x()
{
  test_dg_maxwell_absorb(1, 1, false);
  test_dg_maxwell_absorb(1, 2, false);
  test_dg_maxwell_absorb(1, 3, false);
}

void
test_dg_maxwell_absorb_2x()
{
  test_dg_maxwell_absorb(2, 1, false);
  test_dg_maxwell_absorb(2, 2, false);
  test_dg_maxwell_absorb(2, 3, false);
  test_dg_maxwell_absorb(2, 2, true);
  test_dg_maxwell_absorb(2, 3, true);
}

void
test_dg_maxwell_absorb_3x()
{
  test_dg_maxwell_absorb(3, 1, false);
  test_dg_maxwell_absorb(3, 2, false);
  test_dg_maxwell_absorb(3, 2, true);
}

#ifdef GKYL_HAVE_CUDA

int cu_maxwell_test(const struct gkyl_dg_eqn *eqn);

void
test_dg_max_dev()
{
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, 1);

  // Identity position map providing the (required) conf-space Jacobian.
  double lower[] = {0.0}, upper[] = {1.0};
  int cells[] = {8}, ghost[] = {1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);
  struct gkyl_range crange, crange_ext;
  gkyl_create_grid_ranges(&grid, ghost, &crange_ext, &crange);
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = {0};
  struct gkyl_vlasov_position_map *pos_map =
    gkyl_vlasov_position_map_new(&grid, &crange, &crange_ext, &basis, inp_pmap, true);

  // Input structure for building the dg eqn object
  struct gkyl_dg_maxwell_inp inp_dg_maxwell = {
    .cbasis = &basis,
    .crange = &crange,
    .jacob_pos = pos_map->jacob_pos,
    .conf_flux_surf = 0,
    .lightSpeed = 1.0,
    .field_id = GKYL_FIELD_E_B,
    .elcErrorSpeedFactor = 0.5,
    .mgnErrorSpeedFactor = 0.25,
    .use_gpu = true,
  };
  struct gkyl_dg_eqn *eqn = gkyl_dg_maxwell_inew(&inp_dg_maxwell);

  // this is not possible from user code and should NOT be done. This
  // is for testing only
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);

  TEST_CHECK(maxwell->maxwell_data.c == 1.0);
  TEST_CHECK(maxwell->maxwell_data.chi == 0.5);
  TEST_CHECK(maxwell->maxwell_data.gamma == 0.25);

  // call CUDA test
  /* int nfail = cu_maxwell_test(eqn->on_dev); */

  /* TEST_CHECK( nfail == 0 ); */

  gkyl_dg_eqn_release(eqn);
  gkyl_vlasov_position_map_release(pos_map);
}

#endif

TEST_LIST = {
  {"dg_max_ho", test_dg_max_ho},
  {"dg_maxwell_absorb_1x", test_dg_maxwell_absorb_1x},
  {"dg_maxwell_absorb_2x", test_dg_maxwell_absorb_2x},
  {"dg_maxwell_absorb_3x", test_dg_maxwell_absorb_3x},
#ifdef GKYL_HAVE_CUDA
  {"dg_max_dev", test_dg_max_dev},
#endif
  {NULL, NULL}
};
