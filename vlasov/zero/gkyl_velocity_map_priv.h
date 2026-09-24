#include <gkyl_array.h>

// Allocate double array (filled with zeros).
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

/**
 * Function that actually frees memory associated with this
 * object when the number of references has decreased to zero.
 *
 * @param ref Reference counter for this object.
 */
void
gkyl_velocity_map_free(const struct gkyl_ref_count *ref);

#ifdef GKYL_HAVE_CUDA

/**
 * Allocate a velocity map object with pointers to device memory
 * based on a host-side velocity map.
 *
 * @param gvm_ho Host side velocity map object.
 * @return New velocity map object with device pointers.
 */
struct gkyl_velocity_map* gkyl_velocity_map_new_cu_dev(struct gkyl_velocity_map *gvm_ho);

#endif
