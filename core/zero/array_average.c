#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_average_priv.h>
#include <gkyl_array_average.h>

#include <assert.h>

struct gkyl_array_average* gkyl_array_average_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  const struct gkyl_basis *basis_avg, const struct gkyl_range *local, const struct gkyl_range *local_avg,
  const struct gkyl_range *local_avg_ext, const struct gkyl_array *weight, const int *avg_dim, bool use_gpu)
{
  return gkyl_array_average_inew(
    &(struct gkyl_array_average_inp) {
      .grid = grid,
      .basis = *basis,
      .basis_avg = *basis_avg,
      .local = local,
      .local_avg = local_avg,
      .local_avg_ext = local_avg_ext,
      .weight = weight,
      .avg_dim = avg_dim,
      .use_gpu = use_gpu
    }
  );
}

struct gkyl_array_average*
gkyl_array_average_inew(const struct gkyl_array_average_inp *inp)
{
  // works for p <=2 only due to the gkyl_dg_div_op_range call in advance
  assert(inp->basis.poly_order <= 2); 

  // allocate space for new updater.
  struct gkyl_array_average *up = gkyl_malloc(sizeof(struct gkyl_array_average));

  up->use_gpu = inp->use_gpu;
  up->ndim = inp->basis.ndim;
  up->basis = inp->basis;
  up->basis_avg = inp->basis_avg;
  up->local = *inp->local;
  up->local_avg = *inp->local_avg;

  // set up the array of all dimensions that are conserved after the average (=0 for removed)
  // according to the operation input variable
  up->num_avg_dim = 0;
  for (int d=0; d < up->ndim; ++d) {
    up->avg_dim[d] = inp->avg_dim[d];
    up->num_avg_dim += inp->avg_dim[d]; 
  }
  assert(up->num_avg_dim <= up->ndim);

  up->num_dim_remain = up->ndim - up->num_avg_dim;

  for (int d=0; d < up->ndim; ++d) 
    up->dim_remains[d] = 1 - inp->avg_dim[d];

  int k = 0;
  for (int d=0; d < up->ndim; ++d) {
    if(up->dim_remains[d]){
      up->sub_dir[k] = d;
      k++;
    }
  }

  // compute the inverse of the volume of the averaging space
  up->vol_avg_inv = 1.;
  for (int d = 0; d < up->ndim; d++)
    up->vol_avg_inv *= up->avg_dim[d]? inp->grid->upper[d] - inp->grid->lower[d] : 1.0;
  up->vol_avg_inv = 1./up->vol_avg_inv;

  // compute the cell sub-dimensional volume
  up->subvol = 1.0;
  for (int d=0; d < up->ndim; ++d)
    if (up->avg_dim[d]) up->subvol *= 0.5*inp->grid->dx[d];

  // handle a possible weighted average
  up->isweighted = false;
  up->weight = NULL;
  up->weight_avg = NULL;
  up->div_mem = NULL;
  if (inp->weight) {
    up->isweighted = true;
    // compute the subdim integral of the weight (for volume division after integration)
    up->weight_avg = up->use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis_avg.num_basis, inp->local_avg_ext->volume)
      : gkyl_array_new(GKYL_DOUBLE, up->basis_avg.num_basis, inp->local_avg_ext->volume);
    // create new average routine to integrate the weight
    struct gkyl_array_average_inp inp_integral = {
      .grid = inp->grid,
      .basis = inp->basis,
      .basis_avg = inp->basis_avg,
      .local = inp->local,
      .local_avg = inp->local_avg,
      .weight = NULL, // Recursive call without weights
      .avg_dim = inp->avg_dim,
      .use_gpu = inp->use_gpu
    };
    struct gkyl_array_average *int_w = gkyl_array_average_inew(&inp_integral);
    // run the updater to integrate the weight
    gkyl_array_average_advance(int_w, inp->weight, up->weight_avg);
    gkyl_array_average_release(int_w);
    // multiply by the volume to get the integral and not the average
    gkyl_array_scale(up->weight_avg,1./up->vol_avg_inv);
    // allocate memory to prepare the weak division at the end of the advance routine
    up->div_mem = up->use_gpu? gkyl_dg_bin_op_mem_cu_dev_new(up->local_avg.volume, up->basis_avg.num_basis)
      : gkyl_dg_bin_op_mem_new(up->local_avg.volume, up->basis_avg.num_basis);
    // assign the weight pointer to the input weight array
    up->weight = gkyl_array_acquire(inp->weight);
  } 
  else {
    // assign the weight pointer to the identity weight
    up->weight = up->use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis.num_basis, 1)
      : gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, 1);
    gkyl_array_shiftc(up->weight, pow(sqrt(2.),up->ndim), 0);
    // divide by the total volume through the subvol
    up->subvol *= up->vol_avg_inv;
  }

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_array_average_cu_dev_new(up);
#endif
  
  // choose the kernel that performs the desired operation within the integral.
  gkyl_array_average_choose_kernel(up);

  return up;
}

void gkyl_array_average_advance(const struct gkyl_array_average *up, 
  const struct gkyl_array * fin, struct gkyl_array *avgout)
{

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_array_average_advance_cu(up, fin, avgout);
    return;
  }
#endif

  // clear the array that will contain the result
  gkyl_array_clear_range(avgout, 0.0, &up->local_avg);

  struct gkyl_range_iter iter_cmp, iter_avg;
  struct gkyl_range rng_cmp; // this is the complementary range, sub + cmp = full
  
  // we now loop on the range of the averaged array
  gkyl_range_iter_init(&iter_avg, &up->local_avg);
  while (gkyl_range_iter_next(&iter_avg)) {
    long lidx_avg = gkyl_range_idx(&up->local_avg, iter_avg.idx);
    double *avg_i = gkyl_array_fetch(avgout, lidx_avg);

    // we need to pass the moving index to the deflate operation as a sub dimensional iterator
    int parent_idx[GKYL_MAX_CDIM] = {0};
    int cnter = 0;
    for (int i = 0; i < up->basis.ndim; i++) {
      if (up->dim_remains[i]) {
        parent_idx[i] = iter_avg.idx[cnter];
        cnter ++;
      }
    }

    // set up the complementary range, which is the range through all averaged dimensions
    gkyl_range_deflate(&rng_cmp, &up->local, up->dim_remains, parent_idx);
    gkyl_range_iter_no_split_init(&iter_cmp, &rng_cmp);
    while (gkyl_range_iter_next(&iter_cmp)) {
      long lidx_cmp = gkyl_range_idx(&rng_cmp, iter_cmp.idx);

      const double *fin_i = gkyl_array_cfetch(fin, lidx_cmp);
      const double *win_i = up->isweighted? gkyl_array_cfetch(up->weight, lidx_cmp) : 
        gkyl_array_cfetch(up->weight, 0);

      up->kernel(up->subvol, win_i, fin_i, avg_i);
    }
  }

  // if we provided some weight, we now divide by the integrated weight
  if (up->isweighted)
    gkyl_dg_div_op_range(up->div_mem, &up->basis_avg, 0, avgout, 0, avgout, 0, up->weight_avg, &up->local_avg);

}

void gkyl_array_average_release(struct gkyl_array_average *up)
{
  // release memory associated with this updater.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_cu_free(up->on_dev);
#endif

  gkyl_array_release(up->weight);
  if (up->isweighted) {
    gkyl_array_release(up->weight_avg);
    gkyl_dg_bin_op_mem_release(up->div_mem);
  }

  gkyl_free(up);
}
