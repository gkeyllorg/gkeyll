#include <gkyl_alloc.h>
#include <gkyl_eval_offset_fd.h>

#include <string.h>

struct gkyl_eval_offset_fd {
  
  struct gkyl_rect_grid grid;
  int num_ret_vals; // number of values returned by eval function
  struct gkyl_offset_descr *offsets; // size num_ret_vals
  evalf_t eval; // function to project
  void *ctx; // evaluation context
};

gkyl_eval_offset_fd*
gkyl_eval_offset_fd_new(const struct gkyl_eval_offset_fd_inp *inp)
{
  struct gkyl_eval_offset_fd *up = gkyl_malloc(sizeof(*up));

  up->grid = *inp->grid;
  int num_ret_vals = up->num_ret_vals = inp->num_ret_vals;
  up->eval = inp->eval;
  up->ctx = inp->ctx;

  up->offsets = gkyl_malloc(sizeof(struct gkyl_offset_descr[num_ret_vals]));
  memcpy(up->offsets, inp->offsets, sizeof(struct gkyl_offset_descr[num_ret_vals]));

  return up;
}

static inline void
comp_to_phys(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  for (int d=0; d<ndim; ++d) xout[d] = dx[d]*eta[d]+xc[d];
}

void
gkyl_eval_offset_fd_advance(const gkyl_eval_offset_fd *up,
  double tm, const struct gkyl_range *update_rng, struct gkyl_array *out)
{
  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];

  int num_ret_vals = up->num_ret_vals;
  double fvals[num_ret_vals];
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_rng);
  
  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(&up->grid, iter.idx, xc);

    long lidx = gkyl_range_idx(update_rng, iter.idx);
    double *out_p = gkyl_array_fetch(out, lidx);

    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&up->grid, iter.idx, xc);

    for (int c=0; c<num_ret_vals; ++c) {
      // We need to evaluate the function once for each ret value as
      // each can be on a different location in the cell. This is not
      // too efficient, but likely does not matter.
      double xout[GKYL_MAX_DIM];
      comp_to_phys(up->grid.ndim, up->offsets[c].od_off, up->grid.dx, xc, xout);

      up->eval(tm, xout, fvals, up->ctx);
      out_p[c] = fvals[c];
    }
  }
}

void
gkyl_eval_offset_fd_release(gkyl_eval_offset_fd *up)
{
  gkyl_free(up->offsets);
  gkyl_free(up);
}
