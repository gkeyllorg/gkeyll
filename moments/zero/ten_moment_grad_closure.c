#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_null_comm.h>
#include <gkyl_ten_moment_grad_closure.h>
#include <gkyl_ten_moment_grad_closure_priv.h>

gkyl_ten_moment_grad_closure*
gkyl_ten_moment_grad_closure_new(const struct gkyl_ten_moment_grad_closure_inp *inp)
{
  gkyl_ten_moment_grad_closure *up = gkyl_malloc(sizeof(gkyl_ten_moment_grad_closure));

  up->grid = *(inp->grid);
  up->ndim = up->grid.ndim;
  up->k0 = inp->k0;
  up->cfl = inp->cfl;

  int ndim = inp->update_range->ndim;
  
  create_offsets_vertices(inp->update_range, up->offsets_vertices);
  create_offsets_centers(inp->heat_flux_range, up->offsets_centers);
  
  if (inp->comm)
    up->comm = gkyl_comm_acquire(inp->comm);
  else
    up->comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) { } );

  grad_closure_calc_q_choose(up);
  grad_closure_update_q_choose(up);

  up->cfla = gkyl_malloc(sizeof(double));

  return up;
}

struct gkyl_ten_moment_grad_closure_status
gkyl_ten_moment_grad_closure_advance(const gkyl_ten_moment_grad_closure *gces,
  const struct gkyl_range *heat_flux_range, const struct gkyl_range *update_range,
  const struct gkyl_array *fluid, const struct gkyl_array *em_tot,
  struct gkyl_array *cflrate, double dt, struct gkyl_array *heat_flux,
  struct gkyl_array *rhs)
{
  int ndim = update_range->ndim;
  long sz[] = { 2, 4, 8 };

  double *cfla = gces->cfla;
  double cfl = gces->cfl, cflm = 1.1*cfl;
  double is_cfl_violated = 0.0; // deliberately a double
  
  const double* fluid_d[sz[ndim-1]];
  const double* em_tot_d[sz[ndim-1]];
  double *heat_flux_d;
  const double* heat_flux_up[sz[ndim-1]];
  double *rhs_d;

  struct gkyl_range_iter iter_vertex;
  gkyl_range_iter_init(&iter_vertex, heat_flux_range);
  while (gkyl_range_iter_next(&iter_vertex)) {

    long linc_vertex = gkyl_range_idx(heat_flux_range, iter_vertex.idx);
    long linc_center = gkyl_range_idx(update_range, iter_vertex.idx);

    for (int i=0; i<sz[ndim-1]; ++i) {
      em_tot_d[i] =  gkyl_array_cfetch(em_tot, linc_center + gces->offsets_vertices[i]);
      fluid_d[i] = gkyl_array_cfetch(fluid, linc_center + gces->offsets_vertices[i]);
    }

    heat_flux_d = gkyl_array_fetch(heat_flux, linc_vertex);

    gces->calc_q(gces, fluid_d, gkyl_array_fetch(cflrate, linc_center), dt, heat_flux_d);
  }

  struct gkyl_range_iter iter_center;
  gkyl_range_iter_init(&iter_center, update_range);
  while (gkyl_range_iter_next(&iter_center)) {

    long linc_vertex = gkyl_range_idx(heat_flux_range, iter_center.idx);
    long linc_center = gkyl_range_idx(update_range, iter_center.idx);

    for (int i=0; i<sz[ndim-1]; ++i)
      heat_flux_up[i] = gkyl_array_fetch(heat_flux,
        linc_vertex + gces->offsets_centers[i]);

    rhs_d = gkyl_array_fetch(rhs, linc_center);

    gces->update_q(gces, heat_flux_up, rhs_d);
  }

  gkyl_array_reduce(cfla, cflrate, GKYL_MAX);

  if (cfla[0] > cflm)
    is_cfl_violated = 1.0;

  // compute actual CFL, status & max-speed across all domains
  double red_vars[2] = { cfla[0], is_cfl_violated };
  double red_vars_global[2] = { 0.0, 0.0 };
  gkyl_comm_allreduce(gces->comm, GKYL_DOUBLE, GKYL_MAX, 2, red_vars,
    red_vars_global);

  cfla[0] = red_vars_global[0];
  is_cfl_violated = red_vars_global[1];

  double dt_suggested = dt*cfl/fmax(cfla[0], DBL_MIN);

  if (is_cfl_violated > 0.0) {
    // indicate failure, and return smaller stable time-step
    return (struct gkyl_ten_moment_grad_closure_status) {
      .success = 0,
      .dt_suggested = dt_suggested,
    };
  }
  // on success, suggest only bigger time-step; (Only way dt can
  // reduce is if the update fails. If the code comes here the update
  // succeeded and so we should not allow dt to reduce).
  return (struct gkyl_ten_moment_grad_closure_status) {
    .success = is_cfl_violated > 0.0 ? 0 : 1,
    .dt_suggested = dt_suggested > dt ? dt_suggested : dt,
  };
}

void
gkyl_ten_moment_grad_closure_release(gkyl_ten_moment_grad_closure* up)
{
  gkyl_comm_release(up->comm);
  gkyl_free(up->cfla);
  gkyl_free(up);
}
