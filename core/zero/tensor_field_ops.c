#include <gkyl_tensor_field_ops.h>
#include <gkyl_tensor_field.h>

#include <assert.h>


static void 
tensor_field_raise_or_lower_idx_in_place(struct gkyl_tensor_field *met, int raised_idx, 
  struct gkyl_tensor_field *ten, struct gkyl_tensor_field *mem)
{

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(ten->tdata)) { 
    // make a temporary amount of memory on device
    enum gkyl_tensor_index_loc iloc[GKYL_MAX_DIM];
    for (int i=0; i<GKYL_MAX_DIM; ++i) 
      mem->iloc[i] = ten->iloc[i];

    // set the input to the memory, mem
    tensor_field_raise_or_lower_idx_set_cu(met, raised_idx, ten, mem); 

    // copy mem to ten
    gkyl_tensor_field_copy(ten, mem);
    return; 
  }
#endif

  // iterate over the field of tensors
  for (long i=0; i<ten->size; ++i){
    const double *metric = gkyl_tensor_field_cfetch(met, i);

    // loop over the output tensor indices
    struct gkyl_range_iter iter_tf_out;
    gkyl_range_iter_init(&iter_tf_out, &ten->trange);
    while (gkyl_range_iter_next(&iter_tf_out)) {
      double ten_out_elem = 0; 
      int index_raised = iter_tf_out.idx[raised_idx];
      
      // summed over index, j
      for (int j=0; j<ten->ndim; ++j) {
        int idx_met[GKYL_MAX_DIM] = { index_raised, j };
        
        // Get the tensor element we are indexing from
        int idx_tf[GKYL_MAX_DIM];
        for (int k=0; k<GKYL_MAX_DIM; ++k) 
          idx_tf[k] = (raised_idx != k) ?  iter_tf_out.idx[k] : j;

        double met_elem = gkyl_tensor_field_elem_fetch(met, i, idx_met);
        double ten_elem = gkyl_tensor_field_elem_fetch(ten, i, idx_tf);
        ten_out_elem += met_elem*ten_elem;
      }

      // set the resulting sum in the temporary tensor
      gkyl_tensor_field_elem_set(mem, i, iter_tf_out.idx, ten_out_elem);
    }

    // Set the old tensor to the new tensor values
    struct gkyl_range_iter iter_tf;
    gkyl_range_iter_init(&iter_tf, &ten->trange);
    while (gkyl_range_iter_next(&iter_tf)) {
      const double ten_out_elem = gkyl_tensor_field_elem_fetch(mem, i, iter_tf.idx);
      gkyl_tensor_field_elem_set(ten, i, iter_tf.idx, ten_out_elem);
    }
  }
}

static void 
tensor_field_raise_or_lower_idx_set(const struct gkyl_tensor_field *met, int raised_idx, 
  const struct gkyl_tensor_field *ten, struct gkyl_tensor_field *tensor_out)
{

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(ten->tdata)) { tensor_field_raise_or_lower_idx_set_cu(met, raised_idx, ten, tensor_out); return; }
#endif

  // iterate over the field of tensors
  for (long i=0; i<ten->size; ++i){
    const double *metric = gkyl_tensor_field_cfetch(met, i);

    // loop over the output tensor indices
    struct gkyl_range_iter iter_tf_out;
    gkyl_range_iter_init(&iter_tf_out, &ten->trange);
    while (gkyl_range_iter_next(&iter_tf_out)) {
      double ten_out_elem = 0; 
      int index_raised = iter_tf_out.idx[raised_idx];
      
      // summed over index, j
      for (int j=0; j<ten->ndim; ++j) {
        int idx_met[GKYL_MAX_DIM] = { index_raised, j };
        
        // Get the tensor element we are indexing from
        int idx_tf[GKYL_MAX_DIM];
        for (int k=0; k<GKYL_MAX_DIM; ++k) 
          idx_tf[k] = (raised_idx != k) ?  iter_tf_out.idx[k] : j;

        double met_elem = gkyl_tensor_field_elem_fetch(met, i, idx_met);
        double ten_elem = gkyl_tensor_field_elem_fetch(ten, i, idx_tf);
        ten_out_elem += met_elem*ten_elem;
      }

      // set the resulting sum in the temporary tensor
      gkyl_tensor_field_elem_set(tensor_out, i, iter_tf_out.idx, ten_out_elem);
    }
  }
}

void 
gkyl_tensor_field_lower_idx_in_place(struct gkyl_tensor_field *metric, int lowered_idx, 
  struct gkyl_tensor_field *ten, struct gkyl_tensor_field *mem)
{
  // Check that we are lowering a contravariant index
  assert( ten->iloc[lowered_idx] == GKYL_TENSOR_INDEX_UPPER );
  assert( metric->iloc[0]  == GKYL_TENSOR_INDEX_LOWER );
  assert( metric->iloc[1]  == GKYL_TENSOR_INDEX_LOWER );
  tensor_field_raise_or_lower_idx_in_place(metric, lowered_idx, ten, mem);
  ten->iloc[lowered_idx] = GKYL_TENSOR_INDEX_LOWER;
}

void 
gkyl_tensor_field_raise_idx_in_place(struct gkyl_tensor_field *metric, int raised_idx,
  struct gkyl_tensor_field *ten, struct gkyl_tensor_field *mem)
{
  // Check that we are lowering a contravariant index
  assert( ten->iloc[raised_idx] == GKYL_TENSOR_INDEX_LOWER );
  assert( metric->iloc[0]  == GKYL_TENSOR_INDEX_UPPER );
  assert( metric->iloc[1]  == GKYL_TENSOR_INDEX_UPPER );
  tensor_field_raise_or_lower_idx_in_place(metric, raised_idx, ten, mem);
  ten->iloc[raised_idx] = GKYL_TENSOR_INDEX_UPPER;
}

void 
gkyl_tensor_field_lower_idx_set(const struct gkyl_tensor_field *metric, int lowered_idx,
  const struct gkyl_tensor_field *ten, struct gkyl_tensor_field *ten_out)
{
  // Check that we are lowering a contravariant index
  assert( ten->iloc[lowered_idx] == GKYL_TENSOR_INDEX_UPPER );
  assert( metric->iloc[0]  == GKYL_TENSOR_INDEX_LOWER );
  assert( metric->iloc[1]  == GKYL_TENSOR_INDEX_LOWER );
  tensor_field_raise_or_lower_idx_set(metric, lowered_idx, ten, ten_out);
  ten_out->iloc[lowered_idx] = GKYL_TENSOR_INDEX_LOWER;
}

void 
gkyl_tensor_field_raise_idx_set(const struct gkyl_tensor_field *metric, int raised_idx, 
  const struct gkyl_tensor_field *ten, struct gkyl_tensor_field *ten_out)
{
  // Check that we are lowering a contravariant index
  assert( ten->iloc[raised_idx] == GKYL_TENSOR_INDEX_LOWER );
  assert( metric->iloc[0]  == GKYL_TENSOR_INDEX_UPPER );
  assert( metric->iloc[1]  == GKYL_TENSOR_INDEX_UPPER );
  tensor_field_raise_or_lower_idx_set(metric, raised_idx, ten, ten_out);
  ten_out->iloc[raised_idx] = GKYL_TENSOR_INDEX_UPPER;
}