#pragma once

#if defined(__GNUC__) || defined(__GNUG__)
#if defined(__arm__) || defined(__arm64__) ||  defined(__powerpc64__)
// nothing for arm chips / power9
#else
#include <xmmintrin.h>
#endif
#endif

#if defined(__clang__)
#if defined(__APPLE__)
#if defined(__arm__) || defined(__arm64__)
// nothing for Apple m1 chip
#else
#include <fenv.h>
#endif
#endif
#endif  

/** Disable denormalized floats from occuring */
static void
disable_denorm_float(void)
{
#if defined(__GNUC__) || defined(__GNUG__)
#if defined(__arm__) || defined(__arm64__)  ||  defined(__powerpc64__)
// nothing for arm chips
#else
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif
#endif

#if defined(__clang__)
#if defined(__APPLE__)
#if defined(__arm__) || defined(__arm64__)
// nothing for Apple m1 chip
#else
  fesetenv(FE_DFL_DISABLE_SSE_DENORMS_ENV);
#endif
#endif
#endif  
}
