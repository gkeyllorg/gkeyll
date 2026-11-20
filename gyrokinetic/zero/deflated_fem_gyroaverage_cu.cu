/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_deflated_fem_gyroaverage.h>
#include <gkyl_deflated_fem_gyroaverage_priv.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
}

static int
defl_gyroavg_num_blocks(long num_threads)
{
  // Number of blocks of GKYL_DEFAULT_NUM_THREADS threads needed for num_threads threads.
  return (int)((num_threads + GKYL_DEFAULT_NUM_THREADS - 1) / GKYL_DEFAULT_NUM_THREADS);
}

__global__ static void
defl_gyroavg_set_cu_dev_ptrs(struct gkyl_deflated_fem_gyroaverage *up)
{
  // Set the pointers to the kernels on the device.
  for (int edge = 0; edge < 2; edge++) {
    up->deflate[edge] = defl_gyroavg_choose_deflate_kernel(up->cdim, edge);
  }
}

void
defl_gyroavg_dev_new(struct gkyl_deflated_fem_gyroaverage *up)
{
  struct gkyl_deflated_fem_gyroaverage *up_dev =
    (struct gkyl_deflated_fem_gyroaverage *)gkyl_cu_malloc(sizeof(*up_dev));
  gkyl_cu_memcpy(up_dev, up, sizeof(*up_dev), GKYL_CU_MEMCPY_H2D);
  defl_gyroavg_set_cu_dev_ptrs<<<1, 1>>>(up_dev);
  up->on_dev = up_dev;
}

__global__ static void
defl_gyroavg_deflate_planes_cu_kernel(
  const struct gkyl_deflated_fem_gyroaverage *up, const struct gkyl_array *fin,
  struct gkyl_array *planes_fld
)
{
  // Evaluate a field on every plane. Each plane is the lower z surface of a cell, except the
  // last one which is the upper z surface of the last cell.
  int cdim = up->cdim;
  long plane_vol = up->deflated_local.volume, plane_size = up->deflated_local_ext.volume;
  long num_tot = up->num_planes * plane_vol;
  int plane_idx[GKYL_MAX_DIM], conf_idx[GKYL_MAX_DIM];
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < num_tot;
       tid += gridDim.x * blockDim.x) {
    int plane = tid / plane_vol;
    int edge = plane == up->num_planes - 1 ? 1 : 0;
    // The plane range is a sub-range, so use the inverse index of sub-ranges.
    gkyl_sub_range_inv_idx(&up->deflated_local, tid % plane_vol, plane_idx);
    for (int d = 0; d < cdim - 1; d++) {
      conf_idx[d] = plane_idx[d] - up->deflated_local.lower[d] + up->local.lower[d];
    }
    conf_idx[cdim - 1] = up->local.lower[cdim - 1] + plane - edge;

    const double *fin_d =
      (const double *)gkyl_array_cfetch(fin, gkyl_range_idx(&up->local, conf_idx));
    double *plane_d = (double *)gkyl_array_fetch(
      planes_fld, plane * plane_size + gkyl_range_idx(&up->deflated_local, plane_idx)
    );
    up->deflate[edge](fin_d, plane_d);
  }
}

__global__ static void
defl_gyroavg_deflate_cells_cu_kernel(
  const struct gkyl_deflated_fem_gyroaverage *up, const struct gkyl_array *fin,
  struct gkyl_array *planes_dens
)
{
  // For every cell along z, evaluate the field at each magnetic moment at the lower and upper z
  // surfaces of the cell, f_lo and f_up, and store 2*f_lo+f_up and 2*f_up+f_lo (the product of
  // [f_lo; f_up] with the matrix M of gkyl_deflated_fem_gyroaverage_advance_adjoint) as the
  // fields 2*c and 2*c+1, c being the index of the cell along z.
  int cdim = up->cdim, num_basis = up->num_basis, defl_num_basis = up->defl_num_basis;
  long plane_vol = up->deflated_local.volume, plane_size = up->deflated_local_ext.volume;
  long num_tot = up->num_cells_z * plane_vol;
  int plane_idx[GKYL_MAX_DIM], conf_idx[GKYL_MAX_DIM];
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < num_tot;
       tid += gridDim.x * blockDim.x) {
    int cell_z = tid / plane_vol;
    // The plane range is a sub-range, so use the inverse index of sub-ranges.
    gkyl_sub_range_inv_idx(&up->deflated_local, tid % plane_vol, plane_idx);
    for (int d = 0; d < cdim - 1; d++) {
      conf_idx[d] = plane_idx[d] - up->deflated_local.lower[d] + up->local.lower[d];
    }
    conf_idx[cdim - 1] = up->local.lower[cdim - 1] + cell_z;

    const double *fin_d =
      (const double *)gkyl_array_cfetch(fin, gkyl_range_idx(&up->local, conf_idx));
    long plane_linidx = gkyl_range_idx(&up->deflated_local, plane_idx);
    double *dens_lo_d =
      (double *)gkyl_array_fetch(planes_dens, (2 * cell_z) * plane_size + plane_linidx);
    double *dens_up_d =
      (double *)gkyl_array_fetch(planes_dens, (2 * cell_z + 1) * plane_size + plane_linidx);
    for (int n = 0; n < up->num_mu; n++) {
      double fld_lo[DEFL_GYROAVG_MAX_DEFL_BASIS], fld_up[DEFL_GYROAVG_MAX_DEFL_BASIS];
      up->deflate[0](&fin_d[n * num_basis], fld_lo);
      up->deflate[1](&fin_d[n * num_basis], fld_up);
      for (int kd = 0; kd < defl_num_basis; kd++) {
        dens_lo_d[n * defl_num_basis + kd] = 2.0 * fld_lo[kd] + fld_up[kd];
        dens_up_d[n * defl_num_basis + kd] = 2.0 * fld_up[kd] + fld_lo[kd];
      }
    }
  }
}

__global__ static void
defl_gyroavg_deflate_cell_surfs_cu_kernel(
  const struct gkyl_deflated_fem_gyroaverage *up, const struct gkyl_array *fin,
  struct gkyl_array *planes
)
{
  // For every cell along z, evaluate a field (num_basis components) at the lower and upper z
  // surfaces of the cell, and store them as the fields 2*c and 2*c+1, c being the index of the
  // cell along z.
  int cdim = up->cdim;
  long plane_vol = up->deflated_local.volume, plane_size = up->deflated_local_ext.volume;
  long num_tot = up->num_cells_z * plane_vol;
  int plane_idx[GKYL_MAX_DIM], conf_idx[GKYL_MAX_DIM];
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < num_tot;
       tid += gridDim.x * blockDim.x) {
    int cell_z = tid / plane_vol;
    // The plane range is a sub-range, so use the inverse index of sub-ranges.
    gkyl_sub_range_inv_idx(&up->deflated_local, tid % plane_vol, plane_idx);
    for (int d = 0; d < cdim - 1; d++) {
      conf_idx[d] = plane_idx[d] - up->deflated_local.lower[d] + up->local.lower[d];
    }
    conf_idx[cdim - 1] = up->local.lower[cdim - 1] + cell_z;

    const double *fin_d =
      (const double *)gkyl_array_cfetch(fin, gkyl_range_idx(&up->local, conf_idx));
    long plane_linidx = gkyl_range_idx(&up->deflated_local, plane_idx);
    up->deflate[0](
      fin_d, (double *)gkyl_array_fetch(planes, (2 * cell_z) * plane_size + plane_linidx)
    );
    up->deflate[1](
      fin_d, (double *)gkyl_array_fetch(planes, (2 * cell_z + 1) * plane_size + plane_linidx)
    );
  }
}

__global__ static void
defl_gyroavg_inflate_cu_kernel(
  const struct gkyl_deflated_fem_gyroaverage *up, int fld_per_cell, double wgt_same,
  double wgt_other, int num_comp, const struct gkyl_array *planes, struct gkyl_array *fout
)
{
  // In every cell, set each of the num_comp components of fout to the field that is linear in z
  // between combinations of the fields at the lower and upper z surfaces of the cell (see
  // defl_gyroavg_inflate_cell). For the c-th cell along z these are the fields
  // fld_per_cell*c and fld_per_cell*c+1 of planes.
  int cdim = up->cdim;
  long plane_size = up->deflated_local_ext.volume;
  int plane_idx[GKYL_MAX_DIM], conf_idx[GKYL_MAX_DIM];
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x; linc < up->local.volume;
       linc += gridDim.x * blockDim.x) {
    // The local range is a sub-range, so use the inverse index of sub-ranges.
    gkyl_sub_range_inv_idx(&up->local, linc, conf_idx);
    for (int d = 0; d < cdim - 1; d++) {
      plane_idx[d] = conf_idx[d] - up->local.lower[d] + up->deflated_local.lower[d];
    }
    long plane_linidx = gkyl_range_idx(&up->deflated_local, plane_idx);
    long fld_lo = fld_per_cell * (conf_idx[cdim - 1] - up->local.lower[cdim - 1]);
    const double *plane_lo_d =
      (const double *)gkyl_array_cfetch(planes, fld_lo * plane_size + plane_linidx);
    const double *plane_up_d =
      (const double *)gkyl_array_cfetch(planes, (fld_lo + 1) * plane_size + plane_linidx);

    double *fout_d = (double *)gkyl_array_fetch(fout, gkyl_range_idx(&up->local, conf_idx));
    for (int n = 0; n < num_comp; n++) {
      defl_gyroavg_inflate_cell(
        up, wgt_same, wgt_other, &plane_lo_d[n * up->defl_num_basis],
        &plane_up_d[n * up->defl_num_basis], &fout_d[n * up->num_basis]
      );
    }
  }
}

void
gkyl_deflated_fem_gyroaverage_advance_cu(
  struct gkyl_deflated_fem_gyroaverage *up, const struct gkyl_array *fin, struct gkyl_array *fout
)
{
  // Gyroaverage the field on every plane.
  long num_plane_cells = up->num_planes * up->deflated_local.volume;
  defl_gyroavg_deflate_planes_cu_kernel<<<
    defl_gyroavg_num_blocks(num_plane_cells), GKYL_DEFAULT_NUM_THREADS>>>(
    up->on_dev, fin->on_dev, up->planes_fld->on_dev
  );
  gkyl_fem_gyroaverage_batch_advance(up->batch, up->planes_fld, up->planes_fld_gyroavg);

  // Interpolate linearly along z between the planes bounding each cell.
  defl_gyroavg_inflate_cu_kernel<<<
    defl_gyroavg_num_blocks(up->local.volume), GKYL_DEFAULT_NUM_THREADS>>>(
    up->on_dev, 1, 1.0, 0.0, up->num_mu, up->planes_fld_gyroavg->on_dev, fout->on_dev
  );
}

void
gkyl_deflated_fem_gyroaverage_advance_dg_cu(
  struct gkyl_deflated_fem_gyroaverage *up, const struct gkyl_array *fin, struct gkyl_array *fout
)
{
  // Gyroaverage the field at the lower and upper z surfaces of every cell along z, with the
  // batch of the adjoint (its fields 2*c and 2*c+1 are those surfaces of the c-th cell); its
  // buffers have the sizes needed (planes_adj: one field per surface, planes_dens: num_mu).
  long num_plane_cells = up->num_cells_z * up->deflated_local.volume;
  defl_gyroavg_deflate_cell_surfs_cu_kernel<<<
    defl_gyroavg_num_blocks(num_plane_cells), GKYL_DEFAULT_NUM_THREADS>>>(
    up->on_dev, fin->on_dev, up->planes_adj->on_dev
  );
  gkyl_fem_gyroaverage_batch_advance(up->batch_adj, up->planes_adj, up->planes_dens);

  // Interpolate linearly along z between the surfaces of each cell.
  defl_gyroavg_inflate_cu_kernel<<<
    defl_gyroavg_num_blocks(up->local.volume), GKYL_DEFAULT_NUM_THREADS>>>(
    up->on_dev, 2, 1.0, 0.0, up->num_mu, up->planes_dens->on_dev, fout->on_dev
  );
}

void
gkyl_deflated_fem_gyroaverage_advance_adjoint_cu(
  struct gkyl_deflated_fem_gyroaverage *up, const struct gkyl_array *fin, struct gkyl_array *fout
)
{
  // See gkyl_deflated_fem_gyroaverage_advance_adjoint: M*[f_lo; f_up] in every cell along z.
  long num_plane_cells = up->num_cells_z * up->deflated_local.volume;
  defl_gyroavg_deflate_cells_cu_kernel<<<
    defl_gyroavg_num_blocks(num_plane_cells), GKYL_DEFAULT_NUM_THREADS>>>(
    up->on_dev, fin->on_dev, up->planes_dens->on_dev
  );
  gkyl_fem_gyroaverage_batch_advance_adjoint(up->batch_adj, up->planes_dens, up->planes_adj);

  // M^{-1} = [2 -1; -1 2]/3.
  defl_gyroavg_inflate_cu_kernel<<<
    defl_gyroavg_num_blocks(up->local.volume), GKYL_DEFAULT_NUM_THREADS>>>(
    up->on_dev, 2, 2.0 / 3.0, -1.0 / 3.0, 1, up->planes_adj->on_dev, fout->on_dev
  );
}
