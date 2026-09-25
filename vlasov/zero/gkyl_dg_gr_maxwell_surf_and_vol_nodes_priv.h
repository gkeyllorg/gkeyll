#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// allocate double array (filled with zeros)
static struct gkyl_array *
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array *a;
  if (on_gpu) {
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  } else {
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  }
  return a;
}

/**
 * Free the surface and volume node arrays.
 * @param ref Reference to the surface and volume node arrays to free.
 */
void gkyl_surf_and_vol_node_arrays_free(const struct gkyl_ref_count *ref);

static void
vol_node_list(struct gkyl_array *nodes, int ndim, int polyorder)
{
  // Compute the tensor product of nodes at Gauss-Legendre Points
  // for the volume
  const double *gl_nodes_1d = gkyl_gauss_ordinates[polyorder + 1];
  double *nodes_vals = gkyl_array_fetch(nodes, 0);
  int nn = (polyorder + 1);

  switch (ndim) {
    case 1:
      for (int i = 0; i < nn; ++i) {
        nodes_vals[i] = gl_nodes_1d[i];
      }
      break;
    case 2:
      for (int i = 0; i < nn; ++i) {
        for (int j = 0; j < nn; ++j) {
          int flat_index = 2 * (nn * i + j);
          nodes_vals[flat_index] = gl_nodes_1d[i];
          nodes_vals[flat_index + 1] = gl_nodes_1d[j];
        }
      }
      break;
    case 3:
      for (int i = 0; i < nn; ++i) {
        for (int j = 0; j < nn; ++j) {
          for (int k = 0; k < nn; ++k) {
            int flat_index = 3 * (nn * nn * i + nn * j + k);
            nodes_vals[flat_index] = gl_nodes_1d[i];
            nodes_vals[flat_index + 1] = gl_nodes_1d[j];
            nodes_vals[flat_index + 2] = gl_nodes_1d[k];
          }
        }
      }
      break;
  }
}

static void
surf_node_list(struct gkyl_array *nodes, int ndim, int dir, int polyorder)
{
  const double *gl_nodes_1d = gkyl_gauss_ordinates[polyorder + 1];
  double *nodes_vals = gkyl_array_fetch(nodes, 0);
  int nn = polyorder + 1;

  // Build nodes in the direction, dir (which is 0 to ndim-1),
  // Only copy one cell per that direction (ensuring no duplication)
  // and place the it at the surface -1 by replacing the value.
  switch (ndim) {
    case 1:
      // Only one lower-face node in 1D at -1
      nodes_vals[0] = -1.0;
      break;

    case 2:
      switch (dir) {
        case 0:
          for (int j = 0; j < nn; ++j) {
            int flat = 2 * j;
            nodes_vals[flat] = -1.0;
            nodes_vals[flat + 1] = gl_nodes_1d[j];
          }
          break;
        case 1:
          for (int i = 0; i < nn; ++i) {
            int flat = 2 * i;
            nodes_vals[flat] = gl_nodes_1d[i];
            nodes_vals[flat + 1] = -1.0;
          }
          break;
      }
      break;

    case 3:
      switch (dir) {
        case 0:
          for (int j = 0; j < nn; ++j) {
            for (int k = 0; k < nn; ++k) {
              int idx = nn * j + k;
              int flat = 3 * idx;
              nodes_vals[flat] = -1.0;
              nodes_vals[flat + 1] = gl_nodes_1d[j];
              nodes_vals[flat + 2] = gl_nodes_1d[k];
            }
          }
          break;
        case 1:
          for (int i = 0; i < nn; ++i) {
            for (int k = 0; k < nn; ++k) {
              int idx = nn * i + k;
              int flat = 3 * idx;
              nodes_vals[flat] = gl_nodes_1d[i];
              nodes_vals[flat + 1] = -1.0;
              nodes_vals[flat + 2] = gl_nodes_1d[k];
            }
          }
          break;
        case 2:
          for (int i = 0; i < nn; ++i) {
            for (int j = 0; j < nn; ++j) {
              int idx = nn * i + j;
              int flat = 3 * idx;
              nodes_vals[flat] = gl_nodes_1d[i];
              nodes_vals[flat + 1] = gl_nodes_1d[j];
              nodes_vals[flat + 2] = -1.0;
            }
          }
          break;
      }
      break;
  }
}
