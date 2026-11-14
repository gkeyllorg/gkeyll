#pragma once

#include <stddef.h>

// Build gate for the optional heavy kernel set (tensor p=2 basis in 2x3v; see
// the top-level Makefile and ./configure --help). Kernel tables wrap the
// affected entries in this macro so the rows are NULL when the set is not
// built, and the wrapper functions that call those kernels are compiled out.
// The equation/moment constructors assert on NULL rows, and the apps check the
// gate up front with an actionable message.

// Tensor p=2 kernels in 2x3v (velocity-space Hamiltonian Vlasov, moments,
// velocity-map division, LBO/BGK and SR/GR conf-space helpers). The
// phase-space Hamiltonian and triad pieces are p=1 only at 2x3v.
#ifdef GKYL_BUILD_VLASOV_TEN_2X3V_P2
#define GKYL_TEN_2X3V_P2(kern) kern
#else
#define GKYL_TEN_2X3V_P2(kern) NULL
#endif
