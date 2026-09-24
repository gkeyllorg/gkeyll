#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_util.h>

// undefine this to use non-aligned memory allocations
#define USE_ALIGNED_ALLOC
// alignment boundary is 32 bytes to be compatible with AVX
static const size_t ARRAY_ALIGN_BND = 32;

#define set_arr_dat_zero_ho(arr, data) \
  for (size_t i=0; i<arr->size*arr->ncomp; ++i) data[i] = 0

#define set_arr_dat_zero_dev(arr, data_ho) \
  for (size_t i=0; i<arr->size*arr->ncomp; ++i) data_ho[i] = 0; \
  gkyl_cu_memcpy(arr->data, data_ho, arr->size*arr->esznc, GKYL_CU_MEMCPY_H2D);

static void*
g_array_alloc(size_t num, size_t sz)
{
#ifdef USE_ALIGNED_ALLOC
  return gkyl_aligned_alloc(ARRAY_ALIGN_BND, num*sz);
#else
  return gkyl_calloc(num, sz);
#endif
}
static void
g_array_free(void* ptr)
{
#ifdef USE_ALIGNED_ALLOC  
  gkyl_aligned_free(ptr);
#else
  gkyl_free(ptr);
#endif  
}

// size in bytes for various data-types
static const size_t array_elem_size[] = {
  [GKYL_INT] = sizeof(int),
  [GKYL_LONG] = sizeof(long),
  [GKYL_FLOAT] = sizeof(float),
  [GKYL_DOUBLE] = sizeof(double),
  [GKYL_USER] = 1,
};

static void
array_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_array *arr = container_of(ref, struct gkyl_array, ref_count);

  if (false == GKYL_IS_ALLOC_EXTERN(arr->flags)) {
    // only free if we allocated memory ourselves
  
    if (GKYL_IS_CU_ALLOC(arr->flags)) {
#ifdef GKYL_HAVE_CUDA 
      cudaStreamDestroy(arr->iostream);
#endif
      gkyl_cu_free(arr->data);
      gkyl_cu_free(arr->on_dev);
    }
    else {
      g_array_free(arr->data);
    }
    
  }
  gkyl_free(arr);
}

// internal method to allocate array
static struct gkyl_array*
array_new(enum gkyl_elem_type type, size_t ncomp, size_t size, bool is_alloc_extern, void *buff)
{
  struct gkyl_array* arr = gkyl_malloc(sizeof(struct gkyl_array));

  arr->type = type;
  arr->elemsz = array_elem_size[type];
  arr->ncomp = ncomp;
  arr->size = size;
  arr->flags = 0;

  if (is_alloc_extern)
    GKYL_SET_ALLOC_EXTERN(arr->flags);
  else
    GKYL_CLEAR_ALLOC_EXTERN(arr->flags);
  
  GKYL_CLEAR_CU_ALLOC(arr->flags);
#ifdef USE_ALIGNED_ALLOC  
  GKYL_SET_ALLOC_ALIGNED(arr->flags);
#else
  GKYL_CLEAR_ALLOC_ALIGNED(arr->flags);
#endif
  
  arr->esznc = arr->elemsz*arr->ncomp;
  arr->data = buff;  
  if (!is_alloc_extern)
    arr->data = g_array_alloc(arr->size, arr->esznc);
  
  arr->ref_count = gkyl_ref_count_init(array_free);

  arr->nthreads = 1;
  arr->nblocks = 1;

  arr->on_dev = arr; // on_dev reference

  if (!is_alloc_extern) {
    // Zero out array elements (not for user-defined type).
    if (type == GKYL_INT) {
      int *dat_p = arr->data;
      set_arr_dat_zero_ho(arr, dat_p);
    }
    else if (type == GKYL_LONG) {
      long *dat_p = arr->data;
      set_arr_dat_zero_ho(arr, dat_p);
    }
    else if (type == GKYL_FLOAT) {
      float *dat_p = arr->data;
      set_arr_dat_zero_ho(arr, dat_p);
    }
    else if (type == GKYL_DOUBLE) {
      double *dat_p = arr->data;
      set_arr_dat_zero_ho(arr, dat_p);
    }
  }

  return arr;
}

struct gkyl_array*
gkyl_array_new(enum gkyl_elem_type type, size_t ncomp, size_t size)
{
  return array_new(type, ncomp, size, false, 0);
}

struct gkyl_array*
gkyl_array_new_from_buff(enum gkyl_elem_type type, size_t ncomp, size_t size, void *buff)
{
  return array_new(type, ncomp, size, true, buff);
}

bool
gkyl_array_is_cu_dev(const struct gkyl_array *arr)
{
  return GKYL_IS_CU_ALLOC(arr->flags);  
}

bool
gkyl_array_is_using_buffer(const struct gkyl_array *arr)
{
  return GKYL_IS_ALLOC_EXTERN(arr->flags);
}

struct gkyl_array*
gkyl_array_copy(struct gkyl_array* dest, const struct gkyl_array* src)
{
  assert(dest->esznc == src->esznc);
  
  long ncopy = src->size < dest->size ? src->size : dest->size;

  bool dest_is_cu_dev = gkyl_array_is_cu_dev(dest);
  bool src_is_cu_dev = gkyl_array_is_cu_dev(src);

  if (src_is_cu_dev) {
    // source is on device
    if (dest_is_cu_dev)
      gkyl_cu_memcpy(dest->data, src->data, ncopy*src->esznc, GKYL_CU_MEMCPY_D2D);
    else
      gkyl_cu_memcpy(dest->data, src->data, ncopy*src->esznc, GKYL_CU_MEMCPY_D2H);
  }
  else {
    // source is on host
    if (dest_is_cu_dev)
      gkyl_cu_memcpy(dest->data, src->data, ncopy*src->esznc, GKYL_CU_MEMCPY_H2D);
    else
      memcpy(dest->data, src->data, ncopy*src->esznc);
  }
  
  return dest;
}

struct gkyl_array*
gkyl_array_copy_async(struct gkyl_array* dest, const struct gkyl_array* src)
{
  assert(dest->esznc == src->esznc);
  
  long ncopy = src->size < dest->size ? src->size : dest->size;

  bool dest_is_cu_dev = gkyl_array_is_cu_dev(dest);
  bool src_is_cu_dev = gkyl_array_is_cu_dev(src);

  if (src_is_cu_dev) {
    // source is on device
    if (dest_is_cu_dev)
      gkyl_cu_memcpy_async(dest->data, src->data, ncopy*src->esznc, GKYL_CU_MEMCPY_D2D, src->iostream);
    else
      gkyl_cu_memcpy_async(dest->data, src->data, ncopy*src->esznc, GKYL_CU_MEMCPY_D2H, src->iostream);
  }
  else {
    // source is on host
    if (dest_is_cu_dev)
      gkyl_cu_memcpy_async(dest->data, src->data, ncopy*src->esznc, GKYL_CU_MEMCPY_H2D, dest->iostream);
    else
      memcpy(dest->data, src->data, ncopy*src->esznc);
  }
  
  return dest;
}

struct gkyl_array*
gkyl_array_clone(const struct gkyl_array* src)
{
  struct gkyl_array* arr = gkyl_malloc(sizeof(struct gkyl_array));

  arr->type = src->type;
  arr->elemsz = src->elemsz;
  arr->ncomp = src->ncomp;
  arr->esznc = src->esznc;
  arr->size = src->size;
  arr->flags = src->flags;

  GKYL_CLEAR_ALLOC_EXTERN(arr->flags);

  if (!GKYL_IS_CU_ALLOC(src->flags)) {
    arr->data = g_array_alloc(arr->size, arr->esznc);
    memcpy(arr->data, src->data, arr->size*arr->esznc);
  }
#ifdef GKYL_HAVE_CUDA
  if (GKYL_IS_CU_ALLOC(src->flags)) {
    arr->nthreads = src->nthreads;
    arr->nblocks = src->nblocks;

    cudaStreamCreate(&arr->iostream);

    arr->data = gkyl_cu_malloc(arr->size*arr->esznc);
    arr->on_dev = gkyl_cu_malloc(sizeof(struct gkyl_array));
    gkyl_cu_memcpy(arr->data, src->data, arr->size*arr->esznc, GKYL_CU_MEMCPY_D2D);
    gkyl_cu_memcpy(arr->on_dev, src->on_dev, sizeof(struct gkyl_array), GKYL_CU_MEMCPY_D2D);
    gkyl_cu_memcpy(&((arr->on_dev)->data), &arr->data, sizeof(void*), GKYL_CU_MEMCPY_H2D);
  }
#endif
  
  arr->ref_count = gkyl_ref_count_init(array_free);
  
  return arr;
}

struct gkyl_array*
gkyl_array_acquire(const struct gkyl_array* arr)
{
  gkyl_ref_count_inc(&arr->ref_count);
  return (struct gkyl_array*) arr;
}

void
gkyl_array_release(const struct gkyl_array* arr)
{
  if (arr)
    gkyl_ref_count_dec(&arr->ref_count);
}

// CUDA specific code

#ifdef GKYL_HAVE_CUDA

struct gkyl_array*
gkyl_array_cu_dev_new(enum gkyl_elem_type type, size_t ncomp, size_t size)
{
  struct gkyl_array* arr = gkyl_malloc(sizeof(struct gkyl_array));

  arr->type = type;
  arr->elemsz = array_elem_size[type];
  arr->ncomp = ncomp;
  arr->size = size;
  arr->flags = 0;

  GKYL_CLEAR_ALLOC_EXTERN(arr->flags);
  GKYL_SET_CU_ALLOC(arr->flags);
  GKYL_CLEAR_ALLOC_ALIGNED(arr->flags);
  
  arr->esznc = arr->elemsz*arr->ncomp;
  arr->ref_count = gkyl_ref_count_init(array_free);
  arr->data = gkyl_cu_malloc(arr->size*arr->esznc);
  arr->nthreads = GKYL_DEFAULT_NUM_THREADS;
  arr->nblocks = gkyl_int_div_up(arr->size*arr->ncomp, arr->nthreads);

  cudaStreamCreate(&arr->iostream);

  // create a clone of the struct arr->on_dev that lives on the device,
  // so that the whole arr->on_dev struct can be passed to a device kernel
  arr->on_dev = gkyl_cu_malloc(sizeof(struct gkyl_array));
  gkyl_cu_memcpy(arr->on_dev, arr, sizeof(struct gkyl_array), GKYL_CU_MEMCPY_H2D);
  // set device-side data pointer in arr->on_dev to arr->data 
  // (which is the host-side pointer to the device data)
  gkyl_cu_memcpy(&((arr->on_dev)->data), &arr->data, sizeof(void*), GKYL_CU_MEMCPY_H2D);

  // Zero out array elements (not for user-defined type).
  if (type == GKYL_INT) {
    int *data_ho = gkyl_malloc(arr->size*arr->esznc);
    set_arr_dat_zero_dev(arr, data_ho);
    gkyl_free(data_ho);
  }
  else if (type == GKYL_LONG) {
    long *data_ho = gkyl_malloc(arr->size*arr->esznc);
    set_arr_dat_zero_dev(arr, data_ho);
    gkyl_free(data_ho);
  }
  else if (type == GKYL_FLOAT) {
    float *data_ho = gkyl_malloc(arr->size*arr->esznc);
    set_arr_dat_zero_dev(arr, data_ho);
    gkyl_free(data_ho);
  }
  else if (type == GKYL_DOUBLE) {
    double *data_ho = gkyl_malloc(arr->size*arr->esznc);
    set_arr_dat_zero_dev(arr, data_ho);
    gkyl_free(data_ho);
  }

  return arr;
}

struct gkyl_array*
gkyl_array_cu_host_new(enum gkyl_elem_type type, size_t ncomp, size_t size)
{
  struct gkyl_array* arr = gkyl_cu_malloc_host(sizeof(struct gkyl_array));

  arr->type = type;
  arr->elemsz = array_elem_size[type];
  arr->ncomp = ncomp;
  arr->size = size;
  arr->flags = 0;

  GKYL_CLEAR_ALLOC_EXTERN(arr->flags);
  GKYL_CLEAR_CU_ALLOC(arr->flags);
#ifdef USE_ALIGNED_ALLOC  
  GKYL_SET_ALLOC_ALIGNED(arr->flags);
#else
  GKYL_CLEAR_ALLOC_ALIGNED(arr->flags);
#endif
  
  arr->esznc = arr->elemsz*arr->ncomp;
  arr->data = gkyl_cu_malloc_host(arr->size*arr->esznc);
  arr->ref_count = gkyl_ref_count_init(array_free);

  arr->nthreads = 1;
  arr->nblocks = 1;

  arr->on_dev = arr; // on_dev reference
  
  // Zero out array elements (not for user-defined type).
  if (type == GKYL_INT) {
    int *dat_p = arr->data;
    set_arr_dat_zero_ho(arr, dat_p);
  }
  else if (type == GKYL_LONG) {
    long *dat_p = arr->data;
    set_arr_dat_zero_ho(arr, dat_p);
  }
  else if (type == GKYL_FLOAT) {
    float *dat_p = arr->data;
    set_arr_dat_zero_ho(arr, dat_p);
  }
  else if (type == GKYL_DOUBLE) {
    double *dat_p = arr->data;
    set_arr_dat_zero_ho(arr, dat_p);
  }

  return arr;
}


#else

struct gkyl_array*
gkyl_array_cu_dev_new(enum gkyl_elem_type type, size_t ncomp, size_t size)
{
  assert(false);
  return 0;
}

struct gkyl_array*
gkyl_array_cu_host_new(enum gkyl_elem_type type, size_t ncomp, size_t size)
{
  assert(false);
  return 0;
}

#endif // CUDA specific code
