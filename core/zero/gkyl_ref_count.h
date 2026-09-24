#pragma once

// Implementation essentially as presented by Chris Wellons:
// https://nullprogram.com/blog/2015/02/17/

#include <gkyl_util.h>

/**
 * Object holding use count and pointer to destructor function.
 */
struct gkyl_ref_count {
  void (*free)(const struct gkyl_ref_count* );
  int count;
};

/**
 * Create a new ref object with specified function pointer that
 * deletes the container object.
 *
 * @param free Function pointer to the delete function
 * @return Ref object
 */
static inline struct gkyl_ref_count
gkyl_ref_count_init(void (*free)(const struct gkyl_ref_count* ))
{
  return (struct gkyl_ref_count) {
    .free = free,
    .count = 1,
  };
}

/**
 * Increment use count.
 *
 * @param ref Object to increment.
 */
static inline void
gkyl_ref_count_inc(const struct gkyl_ref_count *ref)
{
  ((struct gkyl_ref_count *)ref)->count++;
}

/**
 * Decrement use count. If use count is zero the free (destructor)
 * function is called.
 *
 * @param ref Object to decrement.
 */
static inline void
gkyl_ref_count_dec(const struct gkyl_ref_count *ref)
{
  if (--((struct gkyl_ref_count *)ref)->count == 0)
    ref->free(ref);
}
