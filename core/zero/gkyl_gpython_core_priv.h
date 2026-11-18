#pragma once

#include <gkyl_gpython_core.h>
#include <gkyl_array.h>
#include <gkyl_basis.h>

/* Opaque shim handles point directly to their native objects. */
#define ARR(a) ((struct gkyl_array *)(a))
#define CARR(a) ((const struct gkyl_array *)(a))
#define BAS(b) ((struct gkyl_basis *)(b))
#define CBAS(b) ((const struct gkyl_basis *)(b))
