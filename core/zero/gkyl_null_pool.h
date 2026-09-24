#pragma once

#include <gkyl_job_pool.h>

/**
 * The pool returned by this object merely runs jobs immediately as
 * they are added. 
 *
 * @param nthreads Number of "threads" to create (no real threading)
 * @return Pointer to new job-pool object
 */
struct gkyl_job_pool* gkyl_null_pool_new(int nthreads);
