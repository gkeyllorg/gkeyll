/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_tensor_field.h>
#include <gkyl_tensor_field_ops.h>
#include <gkyl_util.h>
}


static void
gkyl_get_tensor_field_range_kernel_launch_dims(dim3* dimGrid, dim3* dimBlock, gkyl_range trange, int size)
{
  // Create a 2D thread grid so we launch size*trange.volume number of threads 
  // so we can parallelize over tensor components too
  dimBlock->y = trange.volume; // ncomp *must* be less than 256
  dimGrid->y = 1;
  dimBlock->x = GKYL_DEFAULT_NUM_THREADS/trange.volume;
  dimGrid->x = gkyl_int_div_up(size, dimBlock->x);
}

__global__ static void
tensor_field_raise_or_lower_idx_set_cu_kernel(const struct gkyl_tensor_field *met, int raised_idx, 
  const struct gkyl_tensor_field *ten, struct gkyl_tensor_field *tensor_out)
{

  // iterate over the components of the tensor
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  int iter_tf_out_idx[GKYL_MAX_DIM];
  gkyl_range_inv_idx(&tensor_out->trange, linc2, iter_tf_out_idx); 

  // iterate over the indices of the tensor
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < tensor_out->size; tid += blockDim.x*gridDim.x) {
    const double *metric = gkyl_tensor_field_cfetch(met, tid);

    // loop over the output tensor indices
    double ten_out_elem = 0; 
    int index_raised = iter_tf_out_idx[raised_idx];
    
    // summed over index, j
    for (int j=0; j<ten->ndim; ++j) {
      int idx_met[GKYL_MAX_DIM] = { index_raised, j };
    
      // Get the tensor element we are indexing from
      int idx_tf[GKYL_MAX_DIM];
      for (int k=0; k<GKYL_MAX_DIM; ++k) 
        idx_tf[k] = (raised_idx != k) ?  iter_tf_out_idx[k] : j;

      double met_elem = (double) gkyl_tensor_field_elem_fetch(met, tid, idx_met);
      double ten_elem = (double) gkyl_tensor_field_elem_fetch(ten, tid, idx_tf);
      ten_out_elem += met_elem*ten_elem;
    }

    // set the resulting sum in the temporary tensor
    gkyl_tensor_field_elem_set(tensor_out, tid, iter_tf_out_idx, ten_out_elem);
  }
}

void
tensor_field_raise_or_lower_idx_set_cu(const struct gkyl_tensor_field *met, int raised_idx, 
  const struct gkyl_tensor_field *ten, struct gkyl_tensor_field *tensor_out)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_tensor_field_range_kernel_launch_dims(&dimGrid, &dimBlock, ten->trange, ten->size);

  // ?? There is no met/tensor_out->on_dev at present
  tensor_field_raise_or_lower_idx_set_cu_kernel<<<dimGrid, dimBlock>>>(met->on_dev, raised_idx, ten->on_dev, tensor_out->on_dev);
}