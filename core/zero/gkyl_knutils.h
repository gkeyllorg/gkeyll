#pragma once

#include <gkyl_ref_count.h>

#include <stdbool.h>
#include <stdint.h>

// Struct to hold data passed to KANN functions: nvec float vectors
// with N elements each. For GPU usage, create a device kn_vec with
// gkyl_kn_vec_cu_dev_new and use gkyl_kn_vec_copy to transfer data.
struct gkyl_kn_vec {
  int nvec; // number of vectors
  int N; // length of each vector
  float **vals; // host pointers for indexing (NULL on device variant)
  float *data; // contiguous store (device mem when CU_ALLOC set)

  uint32_t flags;
  struct gkyl_kn_vec *on_dev; // device-resident struct clone

  struct gkyl_ref_count ref_count;
};

/**
 * Create a new host kn_vec to hold @a nvec vectors, each of size @a N.
 *
 * @param nvec Number of vectors to create
 * @param N Length of each vector
 * @return New kn_vec
 */
struct gkyl_kn_vec* gkyl_kn_vec_new(int nvec, int N);

/**
 * Create a new device kn_vec. The data member lives on GPU.
 * vals is NULL on the device variant.
 *
 * @param nvec Number of vectors
 * @param N Length of each vector
 * @return New device kn_vec
 */
struct gkyl_kn_vec* gkyl_kn_vec_cu_dev_new(int nvec, int N);

/**
 * Copy data between kn_vecs. Direction (H2D, D2H, D2D, H2H) is
 * inferred from the flags on src and dest.
 *
 * @param dest Destination kn_vec
 * @param src Source kn_vec
 * @return dest
 */
struct gkyl_kn_vec* gkyl_kn_vec_copy(struct gkyl_kn_vec *dest,
  const struct gkyl_kn_vec *src);

/**
 * Check if a kn_vec has device data.
 *
 * @param vec Vector to check
 * @return true if device data
 */
bool gkyl_kn_vec_is_cu_dev(const struct gkyl_kn_vec *vec);

/**
 * Acquire pointer to kn_vec.
 *
 * @param vec Vector to which reference is needed
 * @return Pointer to acquired array
 */
struct gkyl_kn_vec* gkyl_kn_vec_acquire(const struct gkyl_kn_vec* vec);

/**
 * Release kn_vec.
 *
 * @param vec Vector to release.
 */
void gkyl_kn_vec_release(struct gkyl_kn_vec *vec);
