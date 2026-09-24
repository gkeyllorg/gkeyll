#include <gkyl_job_pool.h>

bool
gkyl_job_pool_add_work(const struct gkyl_job_pool *jp, jp_work_func func, void *ctx)
{
  return jp->add_work(jp, func, ctx);
}

void
gkyl_job_pool_wait(const struct gkyl_job_pool *jp)
{
  jp->wait(jp);
}

struct gkyl_job_pool*
gkyl_job_pool_acquire(const struct gkyl_job_pool *jp)
{
  gkyl_ref_count_inc(&jp->ref_count);
  return (struct gkyl_job_pool*) jp;
}

void
gkyl_job_pool_release(const struct gkyl_job_pool* jp)
{
  gkyl_ref_count_dec(&jp->ref_count);
}


