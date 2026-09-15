#pragma once

#include <stddef.h>

// Build gates for the optional heavy kernel sets (tensor p=1 hybrid basis in
// 2x3v and 3x3v; see the top-level Makefile and ./configure --help). Kernel
// tables wrap the affected entries in these macros so the rows are NULL when
// the set is not built, and the wrapper functions that call those kernels are
// compiled out. The equation/moment constructors assert on NULL rows, and the
// apps check the gates up front with an actionable message.
#ifdef GKYL_BUILD_VLASOV_HYB_2X3V
#define GKYL_HYB_2X3V(kern) kern
#else
#define GKYL_HYB_2X3V(kern) NULL
#endif

#ifdef GKYL_BUILD_VLASOV_HYB_3X3V
#define GKYL_HYB_3X3V(kern) kern
#else
#define GKYL_HYB_3X3V(kern) NULL
#endif

// Defined by the build only when both 3x3v knobs are set.
#ifdef GKYL_BUILD_VLASOV_HYB_3X3V_PHASE
#define GKYL_HYB_3X3V_PHASE(kern) kern
#else
#define GKYL_HYB_3X3V_PHASE(kern) NULL
#endif
