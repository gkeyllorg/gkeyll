#include <assert.h>
#include <stdint.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_grid_priv.h>
#include <gkyl_util.h>
void
gkyl_rect_grid_init(struct gkyl_rect_grid *grid, int ndim,
  const double *lower, const double *upper, const int *cells)
{
//  // MF 2023/07/07: commenting this out because it causes seg faults in g2.
//  *grid = (struct gkyl_rect_grid) { };
  
  grid->ndim = ndim;  
  grid->cellVolume = 1.0;
  for (int i=0; i<ndim; ++i) {
    grid->lower[i] = lower[i];
    grid->upper[i] = upper[i];
    grid->cells[i] = cells[i];
    grid->dx[i] = (upper[i]-lower[i])/cells[i];
    grid->cellVolume *= grid->dx[i];
  }
}

struct gkyl_rect_grid*
gkyl_rect_grid_new(int ndim, const double *lower, const double *upper, const int *cells)
{
  struct gkyl_rect_grid* out = gkyl_calloc(1, sizeof(*out));
  gkyl_rect_grid_init(out, ndim, lower, upper, cells);
  return out;
}

bool
gkyl_rect_grid_cmp(const struct gkyl_rect_grid *grid1, struct gkyl_rect_grid *grid2)
{
  if (grid1->ndim != grid2->ndim)
    return false;

  for (int i=0; i<grid1->ndim; ++i) {
    if (grid1->cells[i] != grid2->cells[i])
      return false;
    if (!gkyl_compare_double(grid1->lower[i], grid2->lower[i], 1e-14))
      return false;
    if (!gkyl_compare_double(grid1->upper[i], grid2->upper[i], 1e-14))
      return false;    
  }
  return true;
}

GKYL_CU_DH
void
gkyl_rect_grid_find_cell(const struct gkyl_rect_grid *grid, const double *point,
  const bool *pick_lower, const int *known_index, int *cell_index){

  int nDim = grid->ndim;
  int search_num = 0;
  int search_dim[GKYL_MAX_DIM];
  int dim_trans[GKYL_MAX_DIM];
  double low, high;
  
  for (int d=0; d<nDim; d++) {
    if (known_index[d] < 0) {
      search_dim[search_num] = d;
      dim_trans[d] = search_num;
      search_num = search_num + 1;
    } else {
      dim_trans[d] = -1;
      cell_index[d] = known_index[d];
      low = grid->lower[d]+(known_index[d]-1)*grid->dx[d];
      high = grid->lower[d]+(known_index[d])*grid->dx[d];
      assert(low<point[d]&&high>point[d]);
    }
  }

  int start_index[GKYL_MAX_DIM], end_index[GKYL_MAX_DIM], mid_index[GKYL_MAX_DIM], new_index[GKYL_MAX_DIM];
  const int *cells = grid -> cells;
  for (int d=0; d<search_num; d++) {
    start_index[d] = 1;
    end_index[d] = cells[search_dim[d]];
    mid_index[d] = 0;
    new_index[d] = 0;
  }

  int plusminus[2] = {-1,1}, low_high_index[2*GKYL_MAX_DIM];
  bool all_less_eq = true;
  double lower_dir[nDim], upper_dir[nDim];
  /* Below we use a binary search. That is, if the i-th coordinate of the point in
   * question is below(above) the ith-coordinate of the current mid-point, we search
   * the lower(upper) half along that direction in the next iteration.
   */
  while (all_less_eq) {
    for (int d=0; d<search_num; d++)
      mid_index[d] = (start_index[d]+end_index[d])/2;  // Integer division intentional

    if (is_in_cell(grid, point, mid_index, dim_trans, known_index)) {
      // Check if neighboring cells also contain this point.
      for (int k=0; k<search_num; k++){
	new_index[k] = mid_index[k];
	low_high_index[search_dim[k]] = mid_index[k];
	low_high_index[nDim+search_dim[k]] = mid_index[k];
      }
      for (int i=0; i<search_num; i++){	
	for (int j=0; j<2; j++){
	  new_index[i] = GKYL_MAX2(GKYL_MIN2(mid_index[i] + plusminus[j],cells[search_dim[i]]),1);
	  if (is_in_cell(grid, point, new_index, dim_trans, known_index))
	    low_high_index[j*nDim+search_dim[i]] = new_index[i];	    	  
	}
	new_index[i] = mid_index[i];
      }
      for (int d=0; d<search_num; d++) {
	if (pick_lower[search_dim[d]])
	  cell_index[search_dim[d]] = low_high_index[search_dim[d]];
	else
	  cell_index[search_dim[d]] = low_high_index[nDim+search_dim[d]];
      }
      break;
    } else {
      in_dir(grid, mid_index, dim_trans, known_index, lower_dir, upper_dir);
      for (int d=0; d<search_num; d++) {
	if (point[search_dim[d]] < lower_dir[search_dim[d]]) {
	  end_index[d] = mid_index[d]-1;
	} else if(point[search_dim[d]] > upper_dir[search_dim[d]]) {
	  start_index[d] = mid_index[d]+1;
	}
      }
    }

    for (int d=0; d<search_num; d++) {
      if (start_index[d]>end_index[d]) {
	all_less_eq = false;
	break;
      }
    }
  }
}

void
gkyl_rect_grid_write(const struct gkyl_rect_grid *grid, const char *nm, FILE *fp)
{
  if (fp == stdout || fp == stderr) {
    fprintf(fp, "%s = { ndim = %d, ", nm, grid->ndim);

    fprintf(fp, " lower = { ");
    for (int d=0; d<grid->ndim; ++d)
      fprintf(fp, "%.9e%c ", grid->lower[d], d==grid->ndim-1 ? ' ' : ',');
    fprintf(fp, "}, ");

    fprintf(fp, "upper = { ");
    for (int d=0; d<grid->ndim; ++d)
      fprintf(fp, "%.9e%c ", grid->upper[d] , d==grid->ndim-1 ? ' ' : ',');
    fprintf(fp, "}, ");

    fprintf(fp, "cells = { ");
    for (int d=0; d<grid->ndim; ++d)
      fprintf(fp, "%d%c ", grid->cells[d] , d==grid->ndim-1 ? ' ' : ',');
    fprintf(fp, "}, ");

    fprintf(fp, "dx = { ");
    for (int d=0; d<grid->ndim; ++d)
      fprintf(fp, "%.9e%c ", grid->upper[d] , d==grid->ndim-1 ? ' ' : ',');
    fprintf(fp, "}, ");

    fprintf(fp, " cellVolume = %.9e, ", grid->cellVolume );
    
    fprintf(fp, " }\n");
    fflush(fp);
  }
  else {
    // Dimension and shape are written as 64 bit integers.
    uint64_t ndim = grid->ndim;
    uint64_t cells[GKYL_MAX_DIM];
    for (int d=0; d<grid->ndim; ++d)
      cells[d] = grid->cells[d];
  
    fwrite(&ndim, sizeof(uint64_t), 1, fp);
    fwrite(cells, sizeof(uint64_t), grid->ndim, fp);
    fwrite(grid->lower, sizeof(double), grid->ndim, fp);
    fwrite(grid->upper, sizeof(double), grid->ndim, fp);
  }
}

bool
gkyl_rect_grid_read(struct gkyl_rect_grid *grid, FILE *fp)
{
  uint64_t ndim = grid->ndim;
  uint64_t cells64[GKYL_MAX_DIM];

  if (1 != fread(&ndim, sizeof(uint64_t), 1, fp))
    return false;
  if (ndim != fread(cells64, sizeof(uint64_t), ndim, fp))
    return false;

  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  if (ndim != fread(lower, sizeof(double), ndim, fp))
    return false;
  if (ndim != fread(upper, sizeof(double), ndim, fp))
    return false;

  // copy into regular int array
  int cells[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d) cells[d] = cells64[d];

  gkyl_rect_grid_init(grid, ndim, lower, upper, cells);

  return true;
}

void
gkyl_rect_grid_release(struct gkyl_rect_grid *grid)
{
  gkyl_free(grid);
}
