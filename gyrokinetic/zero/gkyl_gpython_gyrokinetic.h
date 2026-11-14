#pragma once
#include <gkyl_gpython_core.h>

#define GPYTHON_GYROKINETIC_API_VERSION 2
int gpython_gyrokinetic_api_version(void);

/* Uniform physical velocity coordinates, no ghost cells. Returns an owned
 * configuration-space modal array, or NULL for an unsupported basis, moment,
 * grid, or shape. Uses Gkeyll's stored-distribution normalization unchanged. */
gpython_array *gpython_mom_gyrokinetic(
  const gpython_basis *cbasis, const gpython_basis *pbasis, const double *lower,
  const double *upper, const int *cells, double mass, double charge,
  const gpython_array *distribution, const gpython_array *bmag, const gpython_array *phi,
  const char *moment
);
