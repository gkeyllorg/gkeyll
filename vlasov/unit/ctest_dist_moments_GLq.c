#include <acutest.h>
#include <math.h>
#include <stdio.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_dist_moments.h>
#include <gkyl_gauss_quad_data.h>

struct experiment {
  int cdim, vdim, scenario;
  double vmin, vmax;
};

struct moment_errors {
  double modal[5], nodal[5], gap[5];
};

static void
profiles(const struct experiment *ctx, const double *xpos,
  double *density, double drift[2], double temperature[2])
{
  drift[0] = 0.3;
  drift[1] = -0.2;
  temperature[0] = 0.8;
  temperature[1] = 1.4;

  if (ctx->cdim == 1) {
    *density = ctx->scenario == 0 ? 1.0 : 1.0 + 0.2*sin(xpos[0]);

    if (ctx->scenario == 2) {
      drift[0] += 0.15*sin(xpos[0]);
      drift[1] += 0.1*cos(xpos[0]);
      temperature[0] += 0.15*cos(xpos[0]);
      temperature[1] += 0.2*sin(xpos[0]);
    }
  }
  else {
    *density = ctx->scenario == 0 ? 1.0 : 1.0 + 0.1*sin(xpos[0])*cos(xpos[1]);

    if (ctx->scenario == 2) {
      drift[0] += 0.15*sin(xpos[0])*cos(xpos[1]);
      drift[1] += 0.1*cos(xpos[0])*sin(xpos[1]);
      temperature[0] += 0.15*cos(xpos[0])*sin(xpos[1]);
      temperature[1] += 0.2*sin(xpos[0])*cos(xpos[1]);
    }
  }
}

static void
eval_maxwellian(double time, const double *xn, double *out, void *vctx)
{
  (void) time;
  const struct experiment *ctx = vctx;

  double density, drift[2], temperature[2];
  profiles(ctx, xn, &density, drift, temperature);

  double value = density;
  const double pi = acos(-1.0);

  for (int dim = 0; dim < ctx->vdim; ++dim) {
    double relative_velocity = xn[ctx->cdim+dim] - drift[dim];
    value *= exp(-relative_velocity*relative_velocity/(2.0*temperature[dim]))
      /sqrt(2.0*pi*temperature[dim]);
  }

  out[0] = value;
}

// Exact primitive moments n,u and temperature.
static void
reference_moments(const struct experiment *ctx, const double *xpos, double ref[5])
{
  double density, drift_input[2], temperature[2];
  profiles(ctx, xpos, &density, drift_input, temperature);

  ref[0] = density;
  const double pi = acos(-1.0);

  for (int dim = 0; dim < ctx->vdim; ++dim) {
    double thermal_speed = sqrt(temperature[dim]);
    double zlo = (ctx->vmin-drift_input[dim])/thermal_speed;
    double zhi = (ctx->vmax-drift_input[dim])/thermal_speed;

    double exp_lo = exp(-0.5*zlo*zlo)/sqrt(2.0*pi);
    double exp_hi = exp(-0.5*zhi*zhi)/sqrt(2.0*pi);

    // Integrals of the normalized 1D Gaussian:
    // zeroth   = integral g(v) dv
    // first  = integral v*g(v) dv
    // second = integral v^2*g(v) dv
    double zeroth = 0.5*(erf(zhi/sqrt(2.0))-erf(zlo/sqrt(2.0)));
    double first = drift_input[dim]*zeroth + thermal_speed*(exp_lo-exp_hi);
    double second = (drift_input[dim]*drift_input[dim]+temperature[dim])*zeroth
      + 2.0*drift_input[dim]*thermal_speed*(exp_lo-exp_hi)
      + temperature[dim]*(zlo*exp_lo-zhi*exp_hi);

    double drift = first/zeroth;

    ref[0] *= zeroth;
    ref[1+dim] = drift;
    ref[1+ctx->vdim+dim] = second/zeroth-drift*drift;
  }
}

static struct moment_errors
run_case(struct experiment *ctx, int nx, int nv, int poly_order)
{
  int num_nodes = poly_order + 1;
  const double *nodes = gkyl_gauss_ordinates[num_nodes];
  const double *weights = gkyl_gauss_weights[num_nodes];

  int cdim = ctx->cdim;
  int vdim = ctx->vdim, num_mom = 1+2*vdim;
  int pdim = cdim + vdim;

  double lower[4] = {0.0}, upper[4] = {0.0};
  int cells[4] = {0}, ghosts[4] = {0};
  const double pi = acos(-1.0);
  for (int dim = 0; dim < cdim; ++dim) {
    lower[dim] = 0;
    upper[dim] = 2.0*pi;
    cells[dim] = nx;
  }

  for (int dim = 0; dim < vdim; ++dim) {
    lower[cdim+dim] = ctx->vmin;
    upper[cdim+dim] = ctx->vmax;
    cells[cdim+dim] = nv;
  }

  struct gkyl_rect_grid phase_grid, conf_grid, vel_grid;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_rect_grid_init(&vel_grid, vdim, lower+cdim, upper+cdim, cells+cdim);

  struct gkyl_range phase_ext, phase_local, conf_ext, conf_local;
  struct gkyl_range vel_ext, vel_local;
  gkyl_create_grid_ranges(&phase_grid, ghosts, &phase_ext, &phase_local);
  gkyl_create_grid_ranges(&conf_grid, ghosts, &conf_ext, &conf_local);
  gkyl_create_grid_ranges(&vel_grid, ghosts, &vel_ext, &vel_local);

  struct gkyl_basis phase_basis, conf_basis, vel_basis;
  gkyl_cart_modal_tensor(&phase_basis, pdim, poly_order);
  gkyl_cart_modal_tensor(&conf_basis, cdim, poly_order);
  gkyl_cart_modal_tensor(&vel_basis, vdim, poly_order);

  int ncb = conf_basis.num_basis;

  struct gkyl_array *distf = gkyl_array_new(
    GKYL_DOUBLE, phase_basis.num_basis, phase_ext.volume);
  struct gkyl_array *moms = gkyl_array_new(
    GKYL_DOUBLE, num_mom*ncb, conf_ext.volume);

  // Project the actual Maxwellian using eight quadrature nodes
  // per phase-space dimension.
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(
    &phase_grid, &phase_basis, 8, 1, eval_maxwellian, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &phase_local, distf);

  struct gkyl_vlasov_dist_moments_inp inp = {
    .phase_grid = &phase_grid,
    .vel_grid = &vel_grid,
    .conf_basis = &conf_basis,
    .vel_basis = &vel_basis,
    .phase_basis = &phase_basis,
    .conf_range = &conf_local,
    .conf_range_ext = &conf_ext,
    .vel_range = &vel_local,
    .phase_range = &phase_local,
    .dist_id = vdim == 1
      ? GKYL_DIST_TYPE_LTE : GKYL_DIST_TYPE_MULTI_TEMP_MAX,
    .model_id = GKYL_MODEL_DEFAULT,
    .mass = 1.0,
    .use_gpu = false,
  };

  // Modal primitive moments.
  gkyl_vlasov_dist_moments *up = gkyl_vlasov_dist_moments_inew(&inp);
  gkyl_vlasov_dist_moments_advance(up, &phase_local, &conf_local, distf, moms);

  double modal_err[5] = {0.0};
  double nodal_err[5] = {0.0};
  double gap[5] = {0.0};

  double jacobian = 1.0;
  for (int dim = 0; dim < vdim; ++dim)
    jacobian *= 0.5*vel_grid.dx[dim];

  int num_conf_nodes = 1;
  for (int dim = 0; dim < cdim; ++dim)
    num_conf_nodes *= num_nodes;

  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, &conf_local);

  // Advance through configuration cells.
  while (gkyl_range_iter_next(&citer)) {
    // Pointer for the primitive modal coefficients.
    const double *mc = gkyl_array_cfetch(moms, gkyl_range_idx(&conf_ext, citer.idx));

    double center[2];
    gkyl_rect_grid_cell_center(&conf_grid, citer.idx, center);

    for (int iq = 0; iq < num_conf_nodes; ++iq) {
      int conf_remaining = iq;
      double zconf[2] = {0.0};
      double xpos[2] = {0.0};

      for (int dim = 0; dim < cdim; ++dim) {
        int node_index = conf_remaining % num_nodes;
        conf_remaining /= num_nodes;

        zconf[dim] = nodes[node_index];
        xpos[dim] = center[dim]+0.5*conf_grid.dx[dim]*zconf[dim];
      }

      double ref[5] = {0.0};  // Analytical moments.
      double raw[5] = {0.0};  // Accumulates nodal velocity moments.
      double nodal[5] = {0.0}; // Resulting nodal primitive moemnts.

      reference_moments(ctx, xpos, ref);

      struct gkyl_range_iter viter;
      gkyl_range_iter_init(&viter, &vel_local);

      // Advance through velocity cells at the current configuration GL node.
      while (gkyl_range_iter_next(&viter)) {
        int pidx[4] = {0};

        for (int dim = 0; dim < cdim; ++dim)
          pidx[dim] = citer.idx[dim];

        for (int dim = 0; dim < vdim; ++dim)
          pidx[cdim+dim] = viter.idx[dim];

        // Pointer for modal coefficients of f_h.
        const double *fc = gkyl_array_cfetch(distf, gkyl_range_idx(&phase_ext, pidx));

        double vcenter[2];
        gkyl_rect_grid_cell_center(&vel_grid, viter.idx, vcenter);

        // Tensor-product rule: num_nodes points in 1v, num_nodes^2 in 2v.
        int num_velocity_nodes = 1;
        for (int dim = 0; dim < vdim; ++dim)
          num_velocity_nodes *= num_nodes;

        for (int qnode = 0; qnode < num_velocity_nodes; ++qnode) {
          int vel_remaining = qnode;
          double quadrature_weight = 1.0;
          double zphase[4] = {0.0};
          double velocity[2];

          // Copy the current reference configuration coordinates.
          for (int dim = 0; dim < cdim; ++dim)
            zphase[dim] = zconf[dim];

          for (int dim = 0; dim < vdim; ++dim) {
            int node_index = vel_remaining % num_nodes;
            vel_remaining /= num_nodes;

            double eta = nodes[node_index];

            // Store the reference coordinate and map it to physical velocity.
            zphase[cdim+dim] = eta;
            velocity[dim] = vcenter[dim] + 0.5*vel_grid.dx[dim]*eta;

            quadrature_weight *= weights[node_index];
          }

          // Evaluate f_h from its modal coefficients and apply the velocity
          // mapping Jacobian and tensor-product quadrature weight.
          double weighted_f = jacobian*quadrature_weight
            *phase_basis.eval_expand(zphase, fc);

          // Accumulate M0, each M1_d, and each diagonal M2_dd.
          raw[0] += weighted_f;

          for (int dim = 0; dim < vdim; ++dim) {
            raw[1+dim] += velocity[dim]*weighted_f;
            raw[1+vdim+dim] += velocity[dim]*velocity[dim]*weighted_f;
          }
        }
      }

      // Pointwise primitive conversion at each configuration GL node.
      nodal[0] = raw[0];
      for (int dim = 0; dim < vdim; ++dim) {
        double drift = raw[1+dim]/raw[0];
        nodal[1+dim] = drift;
        nodal[1+vdim+dim] = raw[1+vdim+dim]/raw[0]-drift*drift;
      }

      for (int comp = 0; comp < num_mom; ++comp) {
        double modal = conf_basis.eval_expand(zconf, mc+comp*ncb);

        modal_err[comp] = fmax(modal_err[comp], fabs(modal-ref[comp]));
        nodal_err[comp] = fmax(nodal_err[comp], fabs(nodal[comp]-ref[comp]));
        gap[comp] = fmax(gap[comp], fabs(modal-nodal[comp]));
      }
    }
  }

  struct moment_errors result = {0};
  for (int comp = 0; comp < num_mom; ++comp) {
    result.modal[comp] = modal_err[comp];
    result.nodal[comp] = nodal_err[comp];
    result.gap[comp] = gap[comp];
  }

  gkyl_vlasov_dist_moments_release(up);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_release(distf);
  gkyl_array_release(moms);
  
  return result;
}

static void
test_maxwellian_GLq(void)
{
  enum {
    num_orders = 2,    // Polynomial orders: p=1 and p=2.
    num_cases = 3,     // Dimension cases: 1x1v, 1x2v, and 2x2v.
    num_scenarios = 3, // Profile choices: constant, varying n, varying n/u/T.
    num_levels = 3    // Grid resolutions: (nx,nv) = (2,16), (4,32), (8,64).
  };

  const struct {int cdim, vdim; const char *name;} 
    cases[num_cases] = {
    {1, 1, "LTE 1x1v"},
    {1, 2, "MTM 1x2v"},
    {2, 2, "MTM 2x2v"}
  };

  struct moment_errors
    results[num_orders][num_cases][num_scenarios][num_levels];

  for (int poly_order = 1; poly_order <= num_orders; ++poly_order)
    for (int case_idx = 0; case_idx < num_cases; ++case_idx)
      for (int scenario = 0; scenario < num_scenarios; ++scenario) {
        struct experiment ctx = {
          .cdim = cases[case_idx].cdim,
          .vdim = cases[case_idx].vdim,
          .scenario = scenario,
          .vmin = -10.0,
          .vmax = 10.0
        };

        for (int level = 0; level < num_levels; ++level)
          results[poly_order-1][case_idx][scenario][level] =
            run_case(&ctx, 2 << level, 16 << level, poly_order);
      }

  const char *scenarios[num_scenarios] = {
    "constant", "varying density", "varying n,u,T"
  };

  printf("\nOverall results: maximum absolute errors at configuration GL nodes\n");
  printf("nx and nv are cell counts per configuration and velocity direction.\n");
  printf("Temperature = T/m.\n");

  printf("%1s | %-12s | %-15s | %3s | %3s | %-20s | %18s | %18s | %18s\n",
    "p", "Distribution", "Scenario", "nx", "nv", "Quantity",
    "Modal vs analytic", "Nodal vs analytic", "Modal vs nodal");

  for (int poly_order = 1; poly_order <= num_orders; ++poly_order)
    for (int case_idx = 0; case_idx < num_cases; ++case_idx)
      for (int scenario = 0; scenario < num_scenarios; ++scenario)
        for (int level = 0; level < num_levels; ++level) {
          int vdim = cases[case_idx].vdim;
          const struct moment_errors *result =
            &results[poly_order-1][case_idx][scenario][level];

          for (int comp = 0; comp < 1+2*vdim; ++comp) {
            char label[32];

            if (comp == 0)
              snprintf(label, sizeof label, "Density");
            else if (comp <= vdim)
              snprintf(label, sizeof label, "Drift velocity %d", comp-1);
            else
              snprintf(label, sizeof label, "Temperature %d (T/m)", comp-vdim-1);

            printf("%1d | %-12s | %-15s | %3d | %3d | %-20s | %18.8e | %18.8e | %18.8e\n",
              poly_order, cases[case_idx].name, scenarios[scenario],
              2 << level, 16 << level, label,
              result->modal[comp], result->nodal[comp],
              result->gap[comp]);
          }

          printf("\n");
        }

  for (int summary = 0; summary < 2; ++summary) {
    int first_level = summary == 0 ? 0 : num_levels-1;
    if (summary == 0)
      printf("\nLARGEST MODAL-NODAL GAPS - ALL RESOLUTIONS\n");
    else
      printf("\nLARGEST MODAL-NODAL GAPS - FINEST GRID ONLY (nx=%d, nv=%d)\n",
        2 << first_level, 16 << first_level);

    printf("%1s | %-12s | %-20s | %-15s | %3s | %3s | %18s | %18s | %18s\n",
      "p", "Distribution", "Quantity", "Scenario", "nx", "nv",
      "Modal vs nodal", "Modal vs analytic", "Nodal vs analytic");

    for (int poly_order = 1; poly_order <= num_orders; ++poly_order)
      for (int case_idx = 0; case_idx < num_cases; ++case_idx) {
        int vdim = cases[case_idx].vdim;
        int best_scenario = 0, best_level = first_level, best_comp = 0;
        double largest_gap = -1.0;

        for (int scenario = 0; scenario < num_scenarios; ++scenario)
          for (int level = first_level; level < num_levels; ++level)
            for (int comp = 0; comp < 1+2*vdim; ++comp) {
              double gap = results[poly_order-1][case_idx][scenario][level].gap[comp];
              if (gap > largest_gap) {
                largest_gap = gap;
                best_scenario = scenario;
                best_level = level;
                best_comp = comp;
              }
            }

        const struct moment_errors *result =
          &results[poly_order-1][case_idx][best_scenario][best_level];
        char label[32];
        if (best_comp == 0)
          snprintf(label, sizeof label, "Density");
        else if (best_comp <= vdim)
          snprintf(label, sizeof label, "Drift velocity %d", best_comp-1);
        else
          snprintf(label, sizeof label, "Temperature %d (T/m)", best_comp-vdim-1);

        printf("%1d | %-12s | %-20s | %-15s | %3d | %3d | %18.8e | %18.8e | %18.8e\n",
          poly_order, cases[case_idx].name, label, scenarios[best_scenario],
          2 << best_level, 16 << best_level, result->gap[best_comp],
          result->modal[best_comp], result->nodal[best_comp]);
      }
  }
}

TEST_LIST = {
  {"maxwellian_modal_nodal_GLq", test_maxwellian_GLq},
  {NULL, NULL},
};
