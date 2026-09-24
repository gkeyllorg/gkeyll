#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_tensor_field.h>
  
static void
tensor_field_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_tensor_field *tfld = container_of(ref, struct gkyl_tensor_field, ref_count);

  // release the tensorfield data
  gkyl_array_release(tfld->tdata);

  gkyl_free(tfld);
}

struct gkyl_tensor_field *
gkyl_tensor_field_new(size_t rank, size_t ndim, size_t size, const enum gkyl_tensor_index_loc *iloc)
{
  struct gkyl_tensor_field *tfld = gkyl_malloc(sizeof *tfld);

  tfld->ndim = ndim;
  tfld->size = size;
  tfld->rank = rank;

  size_t ncomp = 1;
  int shape[rank];
  for (int i=0; i<rank; ++i) {
    ncomp *= ndim;
    shape[i] = ndim;
  }
  
  tfld->flags = 0;

  tfld->tdata = gkyl_array_new(GKYL_DOUBLE, ncomp, size);
  gkyl_range_init_from_shape(&tfld->trange, rank, shape);

  for (int i=0; i<GKYL_MAX_DIM; ++i) {
    tfld->iloc[i] = iloc[i]; // either upper or lower indices
  }

  tfld->ref_count = gkyl_ref_count_init(tensor_field_free);
  
  return tfld;
}

struct gkyl_tensor_field*
gkyl_tensor_field_acquire(const struct gkyl_tensor_field* tfld)
{
  gkyl_ref_count_inc(&tfld->ref_count);
  return (struct gkyl_tensor_field*) tfld;
}

struct gkyl_tensor_field*
gkyl_tensor_field_copy(struct gkyl_tensor_field* dest, const struct gkyl_tensor_field* src)
{
  for (int i=0; i<GKYL_MAX_DIM; ++i) 
    dest->iloc[i] = src->iloc[i];
  gkyl_array_copy(dest->tdata, src->tdata);
  return dest;
}

void
gkyl_tensor_field_release(const struct gkyl_tensor_field* ten)
{
  if (ten) 
    gkyl_ref_count_dec(&ten->ref_count);
}

bool
gkyl_tensor_field_is_cu_dev(const struct gkyl_tensor_field *tfld)
{
  return GKYL_IS_CU_ALLOC(tfld->flags);  
}


// CUDA specific code

#ifdef GKYL_HAVE_CUDA

struct gkyl_tensor_field*
gkyl_tensor_field_cu_dev_new(size_t rank, size_t ndim, size_t size, const enum gkyl_tensor_index_loc *iloc)
{
  struct gkyl_tensor_field* tfld = gkyl_malloc(sizeof(struct gkyl_tensor_field));

  tfld->rank = rank;
  tfld->ndim = ndim;
  tfld->size = size;

  size_t ncomp = 1;
  int shape[rank];
  for (int i=0; i<rank; ++i) {
    ncomp *= ndim;
    shape[i] = ndim;
  }

  tfld->tdata = gkyl_array_cu_dev_new(GKYL_DOUBLE, ncomp, size); 
  gkyl_range_init_from_shape(&tfld->trange, rank, shape);
  tfld->ref_count = gkyl_ref_count_init(tensor_field_free);

  GKYL_SET_CU_ALLOC(tfld->flags);
  
  for (int i=0; i<GKYL_MAX_DIM; ++i) {
    tfld->iloc[i] = iloc[i]; // either upper or lower indices
  }

  // create a clone of the struct tfld->on_dev that lives on the device,
  // so that the whole tfld->on_dev struct can be passed to a device kernel
  tfld->on_dev = gkyl_cu_malloc(sizeof(struct gkyl_tensor_field));
  gkyl_cu_memcpy(tfld->on_dev, tfld, sizeof(struct gkyl_tensor_field), GKYL_CU_MEMCPY_H2D);
  // set device-side data pointer in tfld->on_dev to tfld->data->on-dev
  // (which is the host-side pointer to the device data)
  gkyl_cu_memcpy(&((tfld->on_dev)->tdata), &tfld->tdata->on_dev, sizeof(void*), GKYL_CU_MEMCPY_H2D);

  return tfld;
}

struct gkyl_tensor_field*
gkyl_tensor_field_cu_host_new(size_t rank, size_t ndim, size_t size, const enum gkyl_tensor_index_loc *iloc)
{
  struct gkyl_tensor_field* tfld = gkyl_cu_malloc_host(sizeof(struct gkyl_tensor_field));

  tfld->rank = rank;
  tfld->ndim = ndim;
  tfld->size = size;

  size_t ncomp = 1;
  int shape[rank];
  for (int i=0; i<rank; ++i) {
    ncomp *= ndim;
    shape[i] = ndim;
  }

  tfld->tdata = gkyl_array_cu_host_new(GKYL_DOUBLE, ncomp, size); // gkyl_cu_malloc_host(tfld->size*tfld->esznc);
  gkyl_range_init_from_shape(&tfld->trange, rank, shape);
  tfld->ref_count = gkyl_ref_count_init(tensor_field_free);

  tfld->flags = 0;

  for (int i=0; i<GKYL_MAX_DIM; ++i) {
    tfld->iloc[i] = iloc[i]; // either upper or lower indices
  }

  tfld->on_dev = tfld; // on_dev reference
  
  return tfld;
}


#else

struct gkyl_tensor_field*
gkyl_tensor_field_cu_dev_new(size_t rank, size_t ndim, size_t size, const enum gkyl_tensor_index_loc *iloc)
{
  assert(false);
  return 0;
}

struct gkyl_tensor_field*
gkyl_tensor_field_cu_host_new(size_t rank, size_t ndim, size_t size, const enum gkyl_tensor_index_loc *iloc)
{
  assert(false);
  return 0;
}

#endif // CUDA specific code