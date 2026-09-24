#pragma once

#include <gkyl_job_pool.h>

/**
 * Create a new thread-pool object
 *
 * @param nthreads Number of threads to create
 * @return Pointer to new job-pool object
 */
struct gkyl_job_pool* gkyl_thread_pool_new(int nthreads);

