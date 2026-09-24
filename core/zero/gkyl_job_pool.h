#pragma once

#include <gkyl_ref_count.h>

// forward declare for use in function pointers
struct gkyl_job_pool;

// Function pointer sig for function that does the actual work 
typedef void (*jp_work_func)(void *ctx);

// Function sig that adds work to the pool
typedef bool (*jp_add_work)(const struct gkyl_job_pool *jp, jp_work_func func, void *ctx);

// Function sig for "wait" function that halts till all jobs have finished
typedef void (*jp_wait)(const struct gkyl_job_pool *jp);

struct gkyl_job_pool {
  int pool_size; // number of worked in pool
  jp_add_work add_work; // function to add work to pool
  jp_wait wait; // function to wait for jobs to finish

  struct gkyl_ref_count ref_count; // reference count  
};

/**
 * Add work to job pool. Work may start immediately on calling this method.
 *
 * @param jp Job-pool object.
 * @param func Pointer to function that does the work
 * @param ctx Context object to add to func
 * @return True if work was added, false otherwise
 */
bool gkyl_job_pool_add_work(const struct gkyl_job_pool *jp, jp_work_func func, void *ctx);

/**
 * Wait till all jobs are completed
 *
 * @param jp Job-pool object.
 */
void gkyl_job_pool_wait(const struct gkyl_job_pool *jp);

/**
 * Acquire pointer to job-pool. Delete using the release()
 * method
 *
 * @param jp Job-pool object.
 */
struct gkyl_job_pool* gkyl_job_pool_acquire(const struct gkyl_job_pool *jp);

/**
 * Delete job-pool object
 *
 * @param jp Object to delete.
 */
void gkyl_job_pool_release(const struct gkyl_job_pool* jp);
