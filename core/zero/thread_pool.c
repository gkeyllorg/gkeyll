#include <gkyl_thread_pool.h>
#include <gkyl_alloc.h>

#include <thpool.h>

struct jp_thread_pool {
  struct gkyl_job_pool jp; // base job-pool object
  threadpool thpool; // thread-pool object
};

static void
thread_pool_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_job_pool *base = container_of(ref, struct gkyl_job_pool, ref_count);
  struct jp_thread_pool *th = container_of(base, struct jp_thread_pool, jp);
  thpool_destroy(th->thpool);
  gkyl_free(th);
}

static bool
thread_pool_add_work(const struct gkyl_job_pool *jp, jp_work_func func, void *ctx)
{
  struct jp_thread_pool *th = container_of(jp, struct jp_thread_pool, jp);
  int status = thpool_add_work(th->thpool, func, ctx);
  return status == 0 ? true : false;
}

static void
thread_pool_wait(const struct gkyl_job_pool *jp)
{
  struct jp_thread_pool *th = container_of(jp, struct jp_thread_pool, jp);
  thpool_wait(th->thpool);
}

struct gkyl_job_pool*
gkyl_thread_pool_new(int nthreads)
{
  struct jp_thread_pool *th = gkyl_malloc(sizeof(struct jp_thread_pool));
  // initialize the actual pool object
  th->thpool = thpool_init(nthreads);  

  th->jp.pool_size = nthreads;
  th->jp.add_work = thread_pool_add_work;
  th->jp.wait = thread_pool_wait;

  // set reference counter
  th->jp.ref_count = gkyl_ref_count_init(thread_pool_free);
    
  return &th->jp;
}
