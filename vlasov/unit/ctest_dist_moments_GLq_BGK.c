#include <acutest.h>
#include <math.h>
#include <stdio.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_bgk_collisions.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_dist_moments.h>

struct experiment {
  int cdim, vdim, scenario;
  double vmin, vmax;
  double separation, collision_frequency;
};

struct moment_errors {
  double modal[5], nodal[5], gap[5];
  // Discrete velocity integrals of C, v0*C, v1*C, and |v|^2*C.
  double modal_conservation[4], nodal_conservation[4];
};

struct equilibrium_context {
  const struct gkyl_rect_grid *conf_grid;
  const struct gkyl_basis *conf_basis;
  const struct gkyl_range *conf_ext;
  const struct gkyl_array *moms;
};

static void
profiles(const struct experiment *ctx, const double *xpos,
  double *density, double drift[2], double *temperature)
{
  *density = ctx->scenario == 0 ? 1.0 : 1.0 + 0.1*sin(xpos[0])*cos(xpos[1]);
  drift[0] = 0.3;
  drift[1] = -0.2;
  *temperature = 0.8;

  if (ctx->scenario == 2) {
    drift[0] += 0.15*sin(xpos[0])*cos(xpos[1]);
    drift[1] += 0.1*cos(xpos[0])*sin(xpos[1]);
    *temperature += 0.15*cos(xpos[0])*sin(xpos[1]);
  }
}

static double
maxwellian(double density, const double drift[2], double temperature,
  const double velocity[2])
{
  const double pi = acos(-1.0);
  double relative_velocity0 = velocity[0]-drift[0];
  double relative_velocity1 = velocity[1]-drift[1];

  return density/(2.0*pi*temperature) * exp(-(relative_velocity0*relative_velocity0
      + relative_velocity1*relative_velocity1)/(2.0*temperature));
}

static void
eval_maxwellian(double time, const double *xn, double *out, void *vctx)
{
  (void) time;
  const struct experiment *ctx = vctx;

  double density, drift[2], temperature;
  profiles(ctx, xn, &density, drift, &temperature);

  double drift_plus[] = {drift[0]+ctx->separation, drift[1]};
  double drift_minus[] = {drift[0]-ctx->separation, drift[1]};
  
  out[0] = 0.5*(maxwellian(density, drift_plus, temperature, xn+ctx->cdim)
    + maxwellian(density, drift_minus, temperature, xn+ctx->cdim));
}

static void
eval_collision_frequency(double time, const double *xn, double *out, void *vctx)
{
  (void) time;
  const struct experiment *ctx = vctx;
  out[0] = ctx->collision_frequency;

  if (ctx->scenario == 2)
    out[0] *= 1.0 + 0.1*sin(xn[0])*cos(xn[1]);
}

// Raw moments of a normalized one-dimensional Gaussian on the finite interval.
static void
gaussian_moments(double drift_input, double temperature,
  double vmin, double vmax, double raw[3])
{
  double thermal_speed = sqrt(temperature);
  double zlo = (vmin-drift_input)/thermal_speed;
  double zhi = (vmax-drift_input)/thermal_speed;

  double exp_lo = exp(-0.5*zlo*zlo)/sqrt(2.0*acos(-1.0));
  double exp_hi = exp(-0.5*zhi*zhi)/sqrt(2.0*acos(-1.0));
 
  // Integrals of the normalized 1D Gaussian:
  // zeroth   = integral g(v) dv
  // first  = integral v*g(v) dv
  // second = integral v^2*g(v) dv
  double zeroth = 0.5*(erf(zhi/sqrt(2.0))-erf(zlo/sqrt(2.0)));
  double first = drift_input*zeroth + thermal_speed*(exp_lo-exp_hi);
  double second = (drift_input*drift_input+temperature)*zeroth
    + 2.0*drift_input*thermal_speed*(exp_lo-exp_hi)
    + temperature*(zlo*exp_lo-zhi*exp_hi);

  raw[0] = zeroth;
  raw[1] = first;
  raw[2] = second;
}

static void
reference_moments(const struct experiment *ctx, const double *xpos, double ref[4])
{
  double density, drift_input[2], temperature;
  profiles(ctx, xpos, &density, drift_input, &temperature);
  double raw[4] = {0.0};

  for (int stream = 0; stream < 2; ++stream) {
    double drift = drift_input[0] + (stream == 0 ? ctx->separation : -ctx->separation);
    double velocity0[3], velocity1[3];

    gaussian_moments(drift, temperature, ctx->vmin, ctx->vmax, velocity0);
    gaussian_moments(drift_input[1], temperature, ctx->vmin, ctx->vmax, velocity1);

    raw[0] += 0.5*density*velocity0[0]*velocity1[0];
    raw[1] += 0.5*density*velocity0[1]*velocity1[0];
    raw[2] += 0.5*density*velocity0[0]*velocity1[1];
    raw[3] += 0.5*density*(velocity0[2]*velocity1[0]+velocity0[0]*velocity1[2]);
  }

  ref[0] = raw[0];
  ref[1] = raw[1]/raw[0];
  ref[2] = raw[2]/raw[0];
  ref[3] = 0.5*(raw[3]/raw[0]-ref[1]*ref[1]-ref[2]*ref[2]);
}

// Reconstruct LTE from modal primitives
static void
eval_equilibrium(double time, const double *xn, double *out, void *vctx)
{
  (void) time;
  const struct equilibrium_context *ctx = vctx;

  int cidx[2];
  for (int dim = 0; dim < 2; ++dim)
    cidx[dim] = 1 + (int) floor((xn[dim]-ctx->conf_grid->lower[dim])/ctx->conf_grid->dx[dim]);
  
  double center[2], zconf[2];
  gkyl_rect_grid_cell_center(ctx->conf_grid, cidx, center);
  for (int dim = 0; dim < 2; ++dim)
    zconf[dim] = 2.0*(xn[dim]-center[dim])/ctx->conf_grid->dx[dim];

  const double *mc = gkyl_array_cfetch(ctx->moms, gkyl_range_idx(ctx->conf_ext, cidx));
  double modal[4];
  for (int comp = 0; comp < 4; ++comp)
    modal[comp] = ctx->conf_basis->eval_expand(zconf, mc+comp*ctx->conf_basis->num_basis);

  out[0] = maxwellian(modal[0], modal+1, modal[3], xn+2);
}

// Enumerate tensor GL nodes and map reference coordinates to physical ones.
static double
quadrature_point(int qnode, int ndim, int num_nodes, const double *nodes,
  const double *weights, const double *center, const double *dx,
  double *zcoord, double *physical)
{
  double quadrature_weight = 1.0;
  for (int dim = 0; dim < ndim; ++dim) {
    int node_index = qnode % num_nodes;
    qnode /= num_nodes;
    zcoord[dim] = nodes[node_index];
    physical[dim] = center[dim]+0.5*dx[dim]*zcoord[dim];
    quadrature_weight *= weights[node_index];
  }
  return quadrature_weight;
}

static struct moment_errors
run_case(struct experiment *ctx, int nx, int nv, int poly_order)
{
  int num_nodes = poly_order+1, num_conf_nodes = num_nodes*num_nodes;
  int num_velocity_nodes = num_nodes*num_nodes;
  const double *nodes = gkyl_gauss_ordinates[num_nodes];
  const double *weights = gkyl_gauss_weights[num_nodes];

  int cdim = ctx->cdim; 
  int vdim = ctx->vdim;
  int pdim = cdim+vdim;

  double lower[] = {0.0, 0.0, ctx->vmin, ctx->vmin};
  double upper[] = {2.0*acos(-1.0), 2.0*acos(-1.0), ctx->vmax, ctx->vmax};
  int cells[] = {nx, nx, nv, nv}, ghosts[] = {0, 0, 0, 0};

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
  
  int ncb = conf_basis.num_basis, npb = phase_basis.num_basis;

  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, npb, phase_ext.volume);
  struct gkyl_array *moms = gkyl_array_new(GKYL_DOUBLE, 4*ncb, conf_ext.volume);
  struct gkyl_array *nu = gkyl_array_new(GKYL_DOUBLE, ncb, conf_ext.volume);
  struct gkyl_array *feq = gkyl_array_new(GKYL_DOUBLE, npb, phase_ext.volume);
  struct gkyl_array *nufM = gkyl_array_new(GKYL_DOUBLE, npb, phase_ext.volume);
  struct gkyl_array *rhs = gkyl_array_new(GKYL_DOUBLE, npb, phase_ext.volume);
  struct gkyl_array *cfl = gkyl_array_new(GKYL_DOUBLE, 1, phase_ext.volume);

  // Project the actual Maxwellian and collision frequency using eight quadrature nodes
  // per phase-space dimension.
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(
    &phase_grid, &phase_basis, 8, 1, eval_maxwellian, ctx);
  gkyl_proj_on_basis *proj_nu = gkyl_proj_on_basis_new(
    &conf_grid, &conf_basis, 8, 1, eval_collision_frequency, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &phase_local, distf);
  gkyl_proj_on_basis_advance(proj_nu, 0.0, &conf_local, nu);

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
    .dist_id = GKYL_DIST_TYPE_LTE, 
    .model_id = GKYL_MODEL_DEFAULT,
    .mass = 1.0, 
    .use_gpu = false
  };

  // Modal primitive moments.
  gkyl_vlasov_dist_moments *up = gkyl_vlasov_dist_moments_inew(&inp);
  gkyl_vlasov_dist_moments_advance(up, &phase_local, &conf_local, distf, moms);

  struct equilibrium_context eq_ctx = {&conf_grid, &conf_basis, &conf_ext, moms};
  gkyl_proj_on_basis *proj_eq = gkyl_proj_on_basis_new(
    &phase_grid, &phase_basis, 8, 1, eval_equilibrium, &eq_ctx);
  gkyl_proj_on_basis_advance(proj_eq, 0.0, &phase_local, feq);
  // Calculate the modal product nu*f_eq.
  gkyl_dg_mul_conf_phase_op_range(&conf_basis, &phase_basis,
    nufM, nu, feq, &conf_local, &phase_local);
  gkyl_bgk_collisions *bgk = gkyl_bgk_collisions_new(&conf_basis, &phase_basis, false);

  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(cfl, 0.0);
  gkyl_bgk_collisions_advance(bgk, &conf_local, &phase_local,
    nu, nufM, distf, false, 0.0, rhs, cfl);

  struct moment_errors result = {0};
  double jacobian = 0.25*vel_grid.dx[0]*vel_grid.dx[1];
  
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, &conf_local);
  
  // Advance through configuration cells.
  while (gkyl_range_iter_next(&citer)) {
    // Pointer for the primitive modal coefficients.
    const double *mc = gkyl_array_cfetch(moms, gkyl_range_idx(&conf_ext, citer.idx));
    const double *nu_c = gkyl_array_cfetch(nu, gkyl_range_idx(&conf_ext, citer.idx));

    double center[2];
    gkyl_rect_grid_cell_center(&conf_grid, citer.idx, center);
    
    for (int iq = 0; iq < num_conf_nodes; ++iq) {
      double zconf[2], xpos[2], ref[4], raw[4] = {0.0}, nodal[4];

      quadrature_point(iq, cdim, num_nodes, nodes, weights,
        center, conf_grid.dx, zconf, xpos);
        
      reference_moments(ctx, xpos, ref);
      double frequency = conf_basis.eval_expand(zconf, nu_c);
      
      struct gkyl_range_iter viter;
      gkyl_range_iter_init(&viter, &vel_local);
      
      // Advance through velocity cells at the current configuration GL node.
      while (gkyl_range_iter_next(&viter)) {
        int pidx[] = {citer.idx[0], citer.idx[1], viter.idx[0], viter.idx[1]};

        // Pointer for modal coefficients of f_h.
        const double *fc = gkyl_array_cfetch(distf, gkyl_range_idx(&phase_ext, pidx));
        
        double vcenter[2];
        gkyl_rect_grid_cell_center(&vel_grid, viter.idx, vcenter);
        
        for (int qnode = 0; qnode < num_velocity_nodes; ++qnode) {
          double zphase[4] = {zconf[0], zconf[1], 0.0, 0.0}, velocity[2];
          double weight = jacobian*quadrature_point(qnode, vdim, num_nodes,
            nodes, weights, vcenter, vel_grid.dx, zphase+cdim, velocity);
          double weighted_f = weight*phase_basis.eval_expand(zphase, fc);
        
          raw[0] += weighted_f;
          raw[1] += velocity[0]*weighted_f;
          raw[2] += velocity[1]*weighted_f;
          raw[3] += (velocity[0]*velocity[0]+velocity[1]*velocity[1])*weighted_f;
        }
      }
      nodal[0] = raw[0];
      nodal[1] = raw[1]/raw[0];
      nodal[2] = raw[2]/raw[0];
      nodal[3] = 0.5*(raw[3]/raw[0]-nodal[1]*nodal[1]-nodal[2]*nodal[2]);

      for (int comp = 0; comp < 4; ++comp) {
        double modal = conf_basis.eval_expand(zconf, mc+comp*ncb);

        result.modal[comp] = fmax(result.modal[comp], fabs(modal-ref[comp]));
        result.nodal[comp] = fmax(result.nodal[comp], fabs(nodal[comp]-ref[comp]));
        result.gap[comp] = fmax(result.gap[comp], fabs(modal-nodal[comp]));
      }

      // Evaluate the BGK collision operator in modal and nodal integrating its velocity moments.
      double modal_conservation[4] = {0.0}, nodal_conservation[4] = {0.0};
      gkyl_range_iter_init(&viter, &vel_local);

      while (gkyl_range_iter_next(&viter)) {
        int pidx[] = {citer.idx[0], citer.idx[1], viter.idx[0], viter.idx[1]};
        long ploc = gkyl_range_idx(&phase_ext, pidx);
        const double *fc = gkyl_array_cfetch(distf, ploc);
        const double *rc = gkyl_array_cfetch(rhs, ploc);

        double vcenter[2];
        gkyl_rect_grid_cell_center(&vel_grid, viter.idx, vcenter);

        for (int qnode = 0; qnode < num_velocity_nodes; ++qnode) {
          double zphase[4] = {zconf[0], zconf[1], 0.0, 0.0};
          double xn[4] = {xpos[0], xpos[1], 0.0, 0.0};
          double weight = jacobian*quadrature_point(qnode, vdim, num_nodes,
            nodes, weights, vcenter, vel_grid.dx, zphase+cdim, xn+cdim);

          double fnode = phase_basis.eval_expand(zphase, fc);

          // Evaluate the modal collision term at this GL point.
          double modal = phase_basis.eval_expand(zphase, rc);

          // Calculate the nodal collision term from the nodal primitives.
          double nodal_rhs = frequency*(maxwellian(nodal[0], nodal+1, nodal[3], xn+cdim)-fnode);

          double original, nu_ref;
          eval_maxwellian(0.0, xn, &original, ctx);
          eval_collision_frequency(0.0, xn, &nu_ref, ctx);
          double reference = nu_ref*(maxwellian(ref[0], ref+1, ref[3], xn+cdim)-original);

          result.modal[4] = fmax(result.modal[4], fabs(modal-reference));
          result.nodal[4] = fmax(result.nodal[4], fabs(nodal_rhs-reference));
          result.gap[4] = fmax(result.gap[4], fabs(modal-nodal_rhs));

          double invariant[] = {1.0, xn[2], xn[3], xn[2]*xn[2]+xn[3]*xn[3]};
          for (int comp = 0; comp < 4; ++comp) {
            modal_conservation[comp] += weight*invariant[comp]*modal;
            nodal_conservation[comp] += weight*invariant[comp]*nodal_rhs;
          }
        }
      }

      for (int comp = 0; comp < 4; ++comp) {
        result.modal_conservation[comp] = fmax(result.modal_conservation[comp], fabs(modal_conservation[comp]));
        result.nodal_conservation[comp] = fmax(result.nodal_conservation[comp], fabs(nodal_conservation[comp]));
      }
    }
  }

  gkyl_bgk_collisions_release(bgk);
  gkyl_vlasov_dist_moments_release(up);
  gkyl_proj_on_basis_release(proj);
  gkyl_proj_on_basis_release(proj_nu);
  gkyl_proj_on_basis_release(proj_eq);
  gkyl_array_release(distf);
  gkyl_array_release(moms);
  gkyl_array_release(nu);
  gkyl_array_release(feq);
  gkyl_array_release(nufM);
  gkyl_array_release(rhs);
  gkyl_array_release(cfl);

  return result;
}

static void
test_maxwellian_GLq_BGK(void)
{
  enum {
    num_orders = 2,    // Polynomial orders: p=1 and p=2.
    num_scenarios = 3, // Constant, varying n, and varying n/u/T/nu.
    num_levels = 2    // (nx,nv) = (2,8), (4,16), per direction.
  };
  struct moment_errors results[num_orders][num_scenarios][num_levels];
  for (int poly_order = 1; poly_order <= num_orders; ++poly_order)
    for (int scenario = 0; scenario < num_scenarios; ++scenario) {
      struct experiment ctx = {
        .cdim = 2, .vdim = 2, .scenario = scenario,
        .vmin = -8.0, .vmax = 8.0,
        .separation = 0.6, .collision_frequency = 0.7
      };
      for (int level = 0; level < num_levels; ++level)
        results[poly_order-1][scenario][level] =
          run_case(&ctx, 2 << level, 8 << level, poly_order);
    }

  const char *scenarios[] = {"constant", "varying density", "varying n,u,T,nu"};
  const char *labels[] = {"Density", "Drift velocity 0", "Drift velocity 1", "Temperature (T/m)", "BGK RHS"};
  const char *invariants[] = {"integral C", "integral v0*C", "integral v1*C", "integral |v|^2*C"};
  printf("\n%1s | %-12s | %-15s | %3s | %3s | %-20s | %18s | %18s | %18s\n",
    "p", "Distribution", "Scenario", "nx", "nv", "Quantity",
    "Modal vs analytic", "Nodal vs analytic", "Modal vs nodal");
  for (int poly_order = 1; poly_order <= num_orders; ++poly_order)
    for (int scenario = 0; scenario < num_scenarios; ++scenario)
      for (int level = 0; level < num_levels; ++level) {
        const struct moment_errors *result = &results[poly_order-1][scenario][level];
        for (int comp = 0; comp < 5; ++comp)
          printf("%1d | %-12s | %-15s | %3d | %3d | %-20s | %18.8e | %18.8e | %18.8e\n",
            poly_order, "LTE BGK 2x2v", scenarios[scenario], 2 << level, 8 << level, labels[comp],
            result->modal[comp], result->nodal[comp], result->gap[comp]);
        printf("\n");
      }

  printf("%1s | %-12s | %-15s | %3s | %3s | %-20s | %18s | %18s\n",
    "p", "Distribution", "Scenario", "nx", "nv", "Quantity",
    "Modal conservation error", "Nodal conservation error");
  for (int poly_order = 1; poly_order <= num_orders; ++poly_order)
    for (int scenario = 0; scenario < num_scenarios; ++scenario)
      for (int level = 0; level < num_levels; ++level) {
        const struct moment_errors *result = &results[poly_order-1][scenario][level];
        for (int comp = 0; comp < 4; ++comp)
          printf("%1d | %-12s | %-15s | %3d | %3d | %-20s | %18.8e | %18.8e\n",
            poly_order, "LTE BGK 2x2v", scenarios[scenario], 2 << level, 8 << level, invariants[comp],
            result->modal_conservation[comp], result->nodal_conservation[comp]);
        printf("\n");
      }
  for (int summary = 0; summary < 2; ++summary) {
    int first_level = summary == 0 ? 0 : num_levels-1;
    if (summary == 0)
      printf("\nLARGEST BGK MODAL-NODAL GAPS - ALL RESOLUTIONS\n");
    else
      printf("\nLARGEST BGK MODAL-NODAL GAPS - FINEST GRID (nx=%d, nv=%d)\n",
        2 << first_level, 8 << first_level);
    printf("%1s | %-12s | %-20s | %-15s | %3s | %3s | %18s | %18s | %18s\n",
      "p", "Distribution", "Quantity", "Scenario", "nx", "nv",
      "Modal vs nodal", "Modal vs analytic", "Nodal vs analytic");
    for (int poly_order = 1; poly_order <= num_orders; ++poly_order) {
      int best_scenario = 0, best_level = first_level;
      double largest_gap = -1.0;
      for (int scenario = 0; scenario < num_scenarios; ++scenario)
        for (int level = first_level; level < num_levels; ++level)
          if (results[poly_order-1][scenario][level].gap[4] > largest_gap) {
            largest_gap = results[poly_order-1][scenario][level].gap[4];
            best_scenario = scenario;
            best_level = level;
          }
      const struct moment_errors *result = &results[poly_order-1][best_scenario][best_level];
      printf("%1d | %-12s | %-20s | %-15s | %3d | %3d | %18.8e | %18.8e | %18.8e\n",
        poly_order, "LTE BGK 2x2v", "BGK RHS", scenarios[best_scenario],
        2 << best_level, 8 << best_level,
        result->gap[4], result->modal[4], result->nodal[4]);
    }
  }

  // Compare scientific-notation exponents, ignoring the coefficients.
  printf("\nConservation magnitude means the power of ten: 4.037e-05 -> e-05.\n");
  printf("Exponent gap = abs(modal exponent - nodal exponent); ties keep the first row.\n");
  printf("Exact zeros have no exponent and are excluded from these summaries.\n");
  for (int summary = 0; summary < 2; ++summary) {
    int first_level = summary == 0 ? 0 : num_levels-1;
    if (summary == 0)
      printf("\nLARGEST CONSERVATION ERROR MAGNITUDE GAPS - ALL RESOLUTIONS\n");
    else
      printf("\nLARGEST CONSERVATION ERROR MAGNITUDE GAPS - FINEST GRID (nx=%d, nv=%d)\n",
        2 << first_level, 8 << first_level);

    printf("%1s | %-12s | %-20s | %-15s | %3s | %3s | %18s | %18s | %18s\n",
      "p", "Distribution", "Quantity", "Scenario", "nx", "nv",
      "Exponent gap", "Modal magnitude", "Nodal magnitude");
    for (int poly_order = 1; poly_order <= num_orders; ++poly_order) {
      int best_scenario = 0, best_level = first_level, best_comp = 0;
      int best_modal_exponent = 0, best_nodal_exponent = 0;
      int largest_gap = -1;
      for (int scenario = 0; scenario < num_scenarios; ++scenario)
        for (int level = first_level; level < num_levels; ++level) {
          const struct moment_errors *result = &results[poly_order-1][scenario][level];
          for (int comp = 0; comp < 4; ++comp) {
            double modal_error = result->modal_conservation[comp];
            double nodal_error = result->nodal_conservation[comp];
            if (modal_error == 0.0 || nodal_error == 0.0)
              continue;

            int modal_exponent = (int) floor(log10(modal_error));
            int nodal_exponent = (int) floor(log10(nodal_error));
            int difference = modal_exponent-nodal_exponent;
            int gap = difference < 0 ? -difference : difference;
            if (gap > largest_gap) {
              largest_gap = gap;
              best_scenario = scenario;
              best_level = level;
              best_comp = comp;
              best_modal_exponent = modal_exponent;
              best_nodal_exponent = nodal_exponent;
            }
          }
        }

      if (largest_gap < 0) {
        printf("%d | LTE BGK 2x2v | No nonzero conservation-error pairs.\n", poly_order);
        continue;
      }
      char modal_magnitude[16], nodal_magnitude[16];
      snprintf(modal_magnitude, sizeof modal_magnitude, "e%+03d", best_modal_exponent);
      snprintf(nodal_magnitude, sizeof nodal_magnitude, "e%+03d", best_nodal_exponent);
      printf("%1d | %-12s | %-20s | %-15s | %3d | %3d | %18d | %18s | %18s\n",
        poly_order, "LTE BGK 2x2v", invariants[best_comp], scenarios[best_scenario],
        2 << best_level, 8 << best_level, largest_gap, modal_magnitude, nodal_magnitude);
    }
  }
}

TEST_LIST = {
  {"maxwellian_modal_nodal_GLq_BGK", test_maxwellian_GLq_BGK},
  {NULL, NULL},
};
