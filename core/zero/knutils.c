#include <gkyl_alloc.h>
#include <gkyl_knutils.h>
#include <gkyl_alloc_flags_priv.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void
kn_vec_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_kn_vec *vec = container_of(ref, struct gkyl_kn_vec, ref_count);

  if (GKYL_IS_CU_ALLOC(vec->flags)) {
    gkyl_cu_free(vec->data);
    gkyl_cu_free(vec->on_dev);
  } else {
    gkyl_free(vec->vals);
    gkyl_free(vec->data);
  }
  gkyl_free(vec);
}

struct gkyl_kn_vec*
gkyl_kn_vec_new(int nvec, int N)
{
  struct gkyl_kn_vec *vec = gkyl_malloc(sizeof(*vec));
  vec->nvec = nvec;
  vec->N = N;
  vec->flags = 0;
  vec->on_dev = 0;

  vec->vals = gkyl_malloc(nvec * sizeof(float*));
  vec->data = gkyl_calloc(nvec * N, sizeof(float));

  vec->vals[0] = &vec->data[0];
  for (int i = 1; i < nvec; ++i)
    vec->vals[i] = vec->vals[i-1] + N;

  vec->ref_count = gkyl_ref_count_init(kn_vec_free);

  return vec;
}

#ifdef GKYL_HAVE_CUDA

struct gkyl_kn_vec*
gkyl_kn_vec_cu_dev_new(int nvec, int N)
{
  struct gkyl_kn_vec *vec = gkyl_malloc(sizeof(*vec));
  vec->nvec = nvec;
  vec->N = N;
  vec->flags = 0;
  vec->vals = 0;

  GKYL_SET_CU_ALLOC(vec->flags);
  vec->data = gkyl_cu_malloc(nvec * N * sizeof(float));

  vec->ref_count = gkyl_ref_count_init(kn_vec_free);

  // device-resident struct clone
  vec->on_dev = gkyl_cu_malloc(sizeof(struct gkyl_kn_vec));
  gkyl_cu_memcpy(vec->on_dev, vec, sizeof(struct gkyl_kn_vec),
    GKYL_CU_MEMCPY_H2D);

  return vec;
}

#else

struct gkyl_kn_vec*
gkyl_kn_vec_cu_dev_new(int nvec, int N)
{
  assert(false);
  return 0;
}

#endif

struct gkyl_kn_vec*
gkyl_kn_vec_copy(struct gkyl_kn_vec *dest, const struct gkyl_kn_vec *src)
{
  assert(dest->nvec == src->nvec);
  assert(dest->N == src->N);

  size_t nbytes = src->nvec * src->N * sizeof(float);
  bool dest_is_cu = GKYL_IS_CU_ALLOC(dest->flags);
  bool src_is_cu = GKYL_IS_CU_ALLOC(src->flags);

  if (src_is_cu) {
    if (dest_is_cu)
      gkyl_cu_memcpy(dest->data, src->data, nbytes, GKYL_CU_MEMCPY_D2D);
    else
      gkyl_cu_memcpy(dest->data, src->data, nbytes, GKYL_CU_MEMCPY_D2H);
  } else {
    if (dest_is_cu)
      gkyl_cu_memcpy(dest->data, src->data, nbytes, GKYL_CU_MEMCPY_H2D);
    else
      memcpy(dest->data, src->data, nbytes);
  }

  return dest;
}

bool
gkyl_kn_vec_is_cu_dev(const struct gkyl_kn_vec *vec)
{
  return GKYL_IS_CU_ALLOC(vec->flags);
}

struct gkyl_kn_vec*
gkyl_kn_vec_acquire(const struct gkyl_kn_vec *vec)
{
  gkyl_ref_count_inc(&vec->ref_count);
  return (struct gkyl_kn_vec*) vec;
}

void
gkyl_kn_vec_release(struct gkyl_kn_vec *vec)
{
  gkyl_ref_count_dec(&vec->ref_count);
}
