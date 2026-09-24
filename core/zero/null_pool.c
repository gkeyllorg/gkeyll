#include <gkyl_null_pool.h>
#include <gkyl_alloc.h>

#include <thpool.h>

static void
null_pool_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_job_pool *np = container_of(ref, struct gkyl_job_pool, ref_count);
  gkyl_free(np);
}

static bool
null_pool_add_work(const struct gkyl_job_pool *jp, jp_work_func func, void *ctx)
{
  func(ctx); // run job immediately
  return true;
}

static void
null_pool_wait(const struct gkyl_job_pool *jp)
{
}

struct gkyl_job_pool*
gkyl_null_pool_new(int nthreads)
{
  struct gkyl_job_pool *np = gkyl_malloc(sizeof(struct gkyl_job_pool));

  np->pool_size = nthreads;
  np->add_work = null_pool_add_work;
  np->wait = null_pool_wait;
  
  // set reference counter
  np->ref_count = gkyl_ref_count_init(null_pool_free);
    
  return np;
}
