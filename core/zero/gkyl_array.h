#pragma once

#include <gkyl_ref_count.h>
#include <gkyl_util.h>
#include <gkyl_elem_type.h>

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <math.h>

/**
 * Array object. This is an untype, undimensioned, reference counted
 * array object. All additional structure is provided else where,
 * mainly by the range object.
 */
struct gkyl_array {
  enum gkyl_elem_type type; // type of data stored in array
  size_t elemsz, ncomp; // size of elements, number of 'components'
  size_t size; // number of indices

  size_t esznc; // elemsz*ncomp
  void *data; // pointer to data
  uint32_t flags;  
  struct gkyl_ref_count ref_count;

  int nthreads, nblocks; // threads per block, number of blocks
  struct gkyl_array *on_dev; // pointer to itself or device data
#ifdef GKYL_HAVE_CUDA
  cudaStream_t iostream;
#else
  int iostream;
#endif
};

/**
 * Create new array. Delete using gkyl_array_release method.
 * 
 * @param type Type of data in array
 * @param ncomp Number of components at each index
 * @param size Number of indices 
 * @return Pointer to newly allocated array.
 */
struct gkyl_array* gkyl_array_new(enum gkyl_elem_type type, size_t ncomp, size_t size);

/**
 * Create new array from per-allocated memory. You must ensure that
 * the memory is sufficiently large for use in the array. Delete using
 * gkyl_array_release method.
 *
 * @param type Type of data in array
 * @param ncomp Number of components at each index
 * @param size Number of indices
 * @param buff Buffer to use for array data
 * @return Pointer to newly allocated array.
 */
struct gkyl_array *gkyl_array_new_from_buff(
  enum gkyl_elem_type type, size_t ncomp, size_t size, void *buff);

/**
 * Create new array with data on NV-GPU. Delete using
 * gkyl_array_release method.
 *
 * NOTE: the data member lives on GPU, but the struct lives on the
 * host.  However, the on_dev member for this cal is set to a device
 * clone of the host struct, and is what should be used to pass to
 * CUDA kernels which require the entire array struct on device.
 * 
 * @param type Type of data in array
 * @param ncomp Number of components at each index
 * @param size Number of indices 
 * @return Pointer to newly allocated array.
 */
struct gkyl_array* gkyl_array_cu_dev_new(enum gkyl_elem_type type, size_t ncomp, size_t size);

/**
 * Create new array with host-pinned data for use with NV-GPU. Delete using
 * gkyl_array_release method.
 *
 * @param type Type of data in array
 * @param ncomp Number of components at each index
 * @param size Number of indices 
 * @return Pointer to newly allocated array.
 */
struct gkyl_array* gkyl_array_cu_host_new(enum gkyl_elem_type type, size_t ncomp, size_t size);

/**
 * Returns true if array lives on NV-GPU.
 *
 * @param arr Array to check
 * @return true of array lives on NV-GPU, false otherwise
 */
bool gkyl_array_is_cu_dev(const struct gkyl_array *arr);

/**
 * Returns true if array uses external buffer for storage.
 *
 * @param arr Array to check
 * @return true if array uses external buffer for storage
 */
bool gkyl_array_is_using_buffer(const struct gkyl_array *arr);

/**
 * Copy into array: pointer to dest array is returned. 'dest' and
 * 'src' must not point to same data.
 *
 * @param dest Destination for copy.
 * @param src Source to copy from.
 * @return dest is returned
 */
struct gkyl_array* gkyl_array_copy(struct gkyl_array* dest,
  const struct gkyl_array* src);

/**
 * Copy into array using async methods for cuda arrays: pointer to dest array is returned. 'dest' and
 * 'src' must not point to same data.
 *
 * @param dest Destination for copy.
 * @param src Source to copy from.
 * @return dest is returned
 */
struct gkyl_array* gkyl_array_copy_async(struct gkyl_array* dest,
  const struct gkyl_array* src);

/**
 * Clone array: pointer to newly created array is returned.
 * 
 * @param arr Array to clone
 * @return Pointer to clone
 */
struct gkyl_array* gkyl_array_clone(const struct gkyl_array* arr);


/**
 * Fetches a pointer to the element stored at the index 'loc'.
 *
 * @param arr Array to fetch from
 * @param loc Element to fetch
 * @return Element at location 'loc'
 */

GKYL_CU_DH
static inline void*
gkyl_array_fetch(struct gkyl_array* arr, long loc)
{
  return ((char*) arr->data) + loc*arr->esznc;
}

/** Same as above, except fetches a constant pointer */
GKYL_CU_DH
static inline const void*
gkyl_array_cfetch(const struct gkyl_array* arr, long loc)
{
  return ((const char*) arr->data) + loc*arr->esznc;
}

/**
 * Acquire pointer to array. The pointer must be released using
 * gkyl_array_release method.
 *
 * @param arr Array to which a pointer is needed
 * @return Pointer to acquired array
 */
struct gkyl_array* gkyl_array_acquire(const struct gkyl_array* arr);

/**
 * Release pointer to array
 *
 * @param arr Array to release.
 */
void gkyl_array_release(const struct gkyl_array* arr);
