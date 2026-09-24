#pragma once

// Private header, not for direct use in user code

// Compute number of elements stored in array 'arr'
#define NELM(arr) (arr->size*arr->ncomp)
// Compute size of 'arr' 
#define NSIZE(arr) (arr->size)
// Compute number of components stored in array 'arr'
#define NCOM(arr) (arr->ncomp)

GKYL_CU_DH
static inline void
array_clear1(long n, double *out, double val)
{
  for (int c=0; c<n; ++c)
    out[c] = val;
}

GKYL_CU_DH
static inline void
array_acc1(long n, double * GKYL_RESTRICT out, double a, const double * GKYL_RESTRICT inp)
{
  for (int c=0; c<n; ++c)
    out[c] += a*inp[c];
}

GKYL_CU_DH
static inline void
array_set1(long n,
  double * GKYL_RESTRICT out, double a, const double * GKYL_RESTRICT inp)
{
  for (int c=0; c<n; ++c)
    out[c] = a*inp[c];
}

GKYL_CU_DH
static inline void
array_set2(long n, long m,
  double * GKYL_RESTRICT out, double a, const double * GKYL_RESTRICT inp)
{
  for (int c=0; c<n; ++c)
    out[c] = a*inp[m+c];
}
