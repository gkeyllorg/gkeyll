#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_deflated_fem_gyroaverage.h>
#include <gkyl_deflated_fem_gyroaverage_priv.h>
#include <gkyl_rect_decomp.h>

#include <assert.h>
#include <math.h>

static void
defl_gyroavg_inflate_new(struct gkyl_deflated_fem_gyroaverage *up)
{
  // Compute the configuration space coefficients of psi_k(x)*h_e(z), psi_k
  // being a basis function of a plane and h_e the linear function of z that is
  // 1 at the lower (e=0) or upper (e=1) z surface of a cell and 0 at the other.
  // Uses a 2-point Gauss quadrature in each direction, exact for these integrands.
  int cdim = up->cdim;
  double ordinates[] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
  for (int e = 0; e < 2; e++) {
    for (int k = 0; k < up->num_basis; k++) {
      for (int kd = 0; kd < up->defl_num_basis; kd++) {
        up->inflate[e][k][kd] = 0.0;
      }
    }
  }

  int num_quad = 1 << cdim;
  for (int q = 0; q < num_quad; q++) {
    double xc[GKYL_MAX_CDIM];
    for (int d = 0; d < cdim; d++) {
      xc[d] = ordinates[(q >> d) & 1];
    }
    double basis_val[DEFL_GYROAVG_MAX_BASIS], defl_basis_val[DEFL_GYROAVG_MAX_DEFL_BASIS];
    up->basis.eval(xc, basis_val);
    up->deflated_basis.eval(xc, defl_basis_val);
    double hat_val[] = {0.5 * (1.0 - xc[cdim - 1]), 0.5 * (1.0 + xc[cdim - 1])};
    for (int e = 0; e < 2; e++) {
      for (int k = 0; k < up->num_basis; k++) {
        for (int kd = 0; kd < up->defl_num_basis; kd++) {
          up->inflate[e][k][kd] += basis_val[k] * defl_basis_val[kd] * hat_val[e];
        }
      }
    }
  }
}

static void
defl_gyroavg_inflate_planes(
  const struct gkyl_deflated_fem_gyroaverage *up, int zidx, double wgt_same, double wgt_other,
  int num_fld, const struct gkyl_array *plane_lo, const struct gkyl_array *plane_up,
  struct gkyl_array *fout
)
{
  // In the cells with z index zidx, set each of the num_fld fields in fout to the field that is
  // linear in z between combinations of the fields on the planes at the lower and upper z
  // surfaces of those cells, plane_lo and plane_up (see defl_gyroavg_inflate_cell).
  int cdim = up->cdim;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->deflated_local);
  while (gkyl_range_iter_next(&iter)) {
    long plane_linidx = gkyl_range_idx(&up->deflated_local, iter.idx);
    const double *plane_lo_d = gkyl_array_cfetch(plane_lo, plane_linidx);
    const double *plane_up_d = gkyl_array_cfetch(plane_up, plane_linidx);

    int conf_idx[GKYL_MAX_CDIM];
    for (int d = 0; d < cdim - 1; d++) {
      conf_idx[d] = iter.idx[d] - up->deflated_local.lower[d] + up->local.lower[d];
    }
    conf_idx[cdim - 1] = zidx;
    double *fout_d = gkyl_array_fetch(fout, gkyl_range_idx(&up->local, conf_idx));
    for (int n = 0; n < num_fld; n++) {
      defl_gyroavg_inflate_cell(
        up, wgt_same, wgt_other, &plane_lo_d[n * up->defl_num_basis],
        &plane_up_d[n * up->defl_num_basis], &fout_d[n * up->num_basis]
      );
    }
  }
}

static struct gkyl_array *
defl_gyroavg_host_copy(const struct gkyl_array *arr)
{
  // Copy of an array (which may be on the device) on the host.
  struct gkyl_array *arr_ho = gkyl_array_new(GKYL_DOUBLE, arr->ncomp, arr->size);
  gkyl_array_copy(arr_ho, arr);
  return arr_ho;
}

static void
defl_gyroavg_host_buffers_new(struct gkyl_deflated_fem_gyroaverage *up)
{
  // Buffers to apply the gyroaverage of one plane at a time, on the host.
  long defl_size = up->deflated_local_ext.volume;
  int num_mu_comp = up->num_mu * up->defl_num_basis;
  up->plane_fld = gkyl_array_new(GKYL_DOUBLE, up->defl_num_basis, defl_size);
  up->plane_dens_lo = gkyl_array_new(GKYL_DOUBLE, num_mu_comp, defl_size);
  up->plane_dens_up = gkyl_array_new(GKYL_DOUBLE, num_mu_comp, defl_size);
  up->plane_dens_buff = gkyl_array_new(GKYL_DOUBLE, num_mu_comp, defl_size);
  up->plane_adj_lo = gkyl_array_new(GKYL_DOUBLE, up->defl_num_basis, defl_size);
  up->plane_adj_up = gkyl_array_new(GKYL_DOUBLE, up->defl_num_basis, defl_size);
  up->plane_fld_gyroavg = gkyl_malloc(sizeof(struct gkyl_array * [up->num_planes]));
  for (int p = 0; p < up->num_planes; p++) {
    up->plane_fld_gyroavg[p] = gkyl_array_new(GKYL_DOUBLE, num_mu_comp, defl_size);
  }
  up->on_dev = up;
}

static void
defl_gyroavg_host_buffers_release(struct gkyl_deflated_fem_gyroaverage *up)
{
  for (int p = 0; p < up->num_planes; p++) {
    gkyl_array_release(up->plane_fld_gyroavg[p]);
  }
  gkyl_free(up->plane_fld_gyroavg);
  gkyl_array_release(up->plane_fld);
  gkyl_array_release(up->plane_dens_lo);
  gkyl_array_release(up->plane_dens_up);
  gkyl_array_release(up->plane_dens_buff);
  gkyl_array_release(up->plane_adj_lo);
  gkyl_array_release(up->plane_adj_up);
}

static void
defl_gyroavg_dev_buffers_new(struct gkyl_deflated_fem_gyroaverage *up)
{
  // Batches of the gyroaverages of the planes, and buffers to apply them to the fields on all
  // the planes at once, on the GPU.
  up->batch = gkyl_fem_gyroaverage_batch_new(up->gyroavg, up->num_planes);
  int num_fld_adj = 2 * up->num_cells_z;
  struct gkyl_fem_gyroaverage **ups_adj =
    gkyl_malloc(sizeof(struct gkyl_fem_gyroaverage *[num_fld_adj]));
  for (int cell_z = 0; cell_z < up->num_cells_z; cell_z++) {
    ups_adj[2 * cell_z] = up->gyroavg[cell_z];
    ups_adj[2 * cell_z + 1] = up->gyroavg[cell_z + 1];
  }
  up->batch_adj = gkyl_fem_gyroaverage_batch_new(ups_adj, num_fld_adj);
  gkyl_free(ups_adj);

  long defl_size = up->deflated_local_ext.volume;
  int num_mu_comp = up->num_mu * up->defl_num_basis;
  up->planes_fld =
    gkyl_array_cu_dev_new(GKYL_DOUBLE, up->defl_num_basis, up->num_planes * defl_size);
  up->planes_fld_gyroavg =
    gkyl_array_cu_dev_new(GKYL_DOUBLE, num_mu_comp, up->num_planes * defl_size);
  up->planes_dens = gkyl_array_cu_dev_new(GKYL_DOUBLE, num_mu_comp, num_fld_adj * defl_size);
  up->planes_adj = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->defl_num_basis, num_fld_adj * defl_size);
#ifdef GKYL_HAVE_CUDA
  defl_gyroavg_dev_new(up);
#endif
}

static void
defl_gyroavg_dev_buffers_release(struct gkyl_deflated_fem_gyroaverage *up)
{
  gkyl_fem_gyroaverage_batch_release(up->batch);
  gkyl_fem_gyroaverage_batch_release(up->batch_adj);
  gkyl_array_release(up->planes_fld);
  gkyl_array_release(up->planes_fld_gyroavg);
  gkyl_array_release(up->planes_dens);
  gkyl_array_release(up->planes_adj);
  gkyl_cu_free(up->on_dev);
}

struct gkyl_deflated_fem_gyroaverage *
gkyl_deflated_fem_gyroaverage_inew(const struct gkyl_deflated_fem_gyroaverage_inp *inp)
{
  struct gkyl_deflated_fem_gyroaverage *up = gkyl_malloc(sizeof(*up));

  up->cdim = inp->grid->ndim;
  up->basis = *inp->basis;
  up->num_basis = up->basis.num_basis;
  up->local = *inp->local;
  up->num_mu = inp->num_mu;
  up->use_gpu = inp->use_gpu;

  // Only implemented for piecewise linear bases in 2x and 3x.
  int cdim = up->cdim;
  assert(cdim == 2 || cdim == 3);
  assert(up->basis.poly_order == 1);

  up->num_cells_z = up->local.upper[cdim - 1] - up->local.lower[cdim - 1] + 1;
  up->num_planes = up->num_cells_z + 1;

  // Grid, ranges and basis of a plane.
  double deflated_lower[GKYL_MAX_CDIM] = {0.0}, deflated_upper[GKYL_MAX_CDIM] = {0.0};
  int deflated_cells[GKYL_MAX_CDIM] = {0};
  for (int d = 0; d < cdim - 1; d++) {
    deflated_lower[d] = inp->grid->lower[d];
    deflated_upper[d] = inp->grid->upper[d];
    deflated_cells[d] = inp->grid->cells[d];
  }
  gkyl_rect_grid_init(&up->deflated_grid, cdim - 1, deflated_lower, deflated_upper, deflated_cells);
  int deflated_nghost[GKYL_MAX_CDIM] = {1, 1, 1};
  gkyl_create_grid_ranges(
    &up->deflated_grid, deflated_nghost, &up->deflated_local_ext, &up->deflated_local
  );
  gkyl_cart_modal_serendip(&up->deflated_basis, cdim - 1, up->basis.poly_order);
  up->defl_num_basis = up->deflated_basis.num_basis;

  defl_gyroavg_inflate_new(up);

  // The gyroaverage on each plane is assembled on the host, with the geometry evaluated on the
  // plane there.
  up->deflator_lo = gkyl_deflate_zsurf_new(&up->basis, &up->deflated_basis, 0, false);
  up->deflator_up = gkyl_deflate_zsurf_new(&up->basis, &up->deflated_basis, 1, false);

  // Each plane is the lower z surface of a cell, except the last one which is
  // the upper z surface of the last cell.
  up->plane_zidx = gkyl_malloc(sizeof(int[up->num_planes]));
  up->plane_deflator = gkyl_malloc(sizeof(struct gkyl_deflate_zsurf * [up->num_planes]));
  for (int p = 0; p < up->num_planes; p++) {
    bool is_last = p == up->num_planes - 1;
    up->plane_zidx[p] = up->local.lower[cdim - 1] + (is_last ? p - 1 : p);
    up->plane_deflator[p] = is_last ? up->deflator_up : up->deflator_lo;
  }

  // Gyroaverage updater of each plane, with the geometry evaluated on the plane.
  long defl_size = up->deflated_local_ext.volume;
  struct gkyl_array *gxx_ho = defl_gyroavg_host_copy(inp->gxx);
  struct gkyl_array *bmag_ho = defl_gyroavg_host_copy(inp->bmag);
  struct gkyl_array *gxy_ho = cdim == 3 ? defl_gyroavg_host_copy(inp->gxy) : 0;
  struct gkyl_array *gyy_ho = cdim == 3 ? defl_gyroavg_host_copy(inp->gyy) : 0;
  struct gkyl_array *plane_gxx = gkyl_array_new(GKYL_DOUBLE, up->defl_num_basis, defl_size);
  struct gkyl_array *plane_gxy = gkyl_array_new(GKYL_DOUBLE, up->defl_num_basis, defl_size);
  struct gkyl_array *plane_gyy = gkyl_array_new(GKYL_DOUBLE, up->defl_num_basis, defl_size);
  struct gkyl_array *plane_bmag = gkyl_array_new(GKYL_DOUBLE, up->defl_num_basis, defl_size);

  up->gyroavg = gkyl_malloc(sizeof(struct gkyl_fem_gyroaverage * [up->num_planes]));
  for (int p = 0; p < up->num_planes; p++) {
    const struct gkyl_deflate_zsurf *deflator = up->plane_deflator[p];
    int zidx = up->plane_zidx[p];
    gkyl_deflate_zsurf_advance(
      deflator, zidx, &up->local, &up->deflated_local, gxx_ho, plane_gxx, 1
    );
    gkyl_deflate_zsurf_advance(
      deflator, zidx, &up->local, &up->deflated_local, bmag_ho, plane_bmag, 1
    );
    if (cdim == 3) {
      gkyl_deflate_zsurf_advance(
        deflator, zidx, &up->local, &up->deflated_local, gxy_ho, plane_gxy, 1
      );
      gkyl_deflate_zsurf_advance(
        deflator, zidx, &up->local, &up->deflated_local, gyy_ho, plane_gyy, 1
      );
    }

    struct gkyl_fem_gyroaverage_inp plane_inp = {
      .grid = &up->deflated_grid,
      .basis = &up->deflated_basis,
      .local = &up->deflated_local,
      .is_periodic = {inp->is_periodic[0], inp->is_periodic[1]},
      .reflect_bnd = inp->reflect_bnd,
      .gxx = plane_gxx,
      .gxy = plane_gxy,
      .gyy = plane_gyy,
      .bmag = plane_bmag,
      .mass = inp->mass,
      .charge = inp->charge,
      .num_mu = inp->num_mu,
      .mu = inp->mu,
      .use_gpu = false,
    };
    up->gyroavg[p] = gkyl_fem_gyroaverage_inew(&plane_inp);
  }

  gkyl_array_release(plane_gxx);
  gkyl_array_release(plane_gxy);
  gkyl_array_release(plane_gyy);
  gkyl_array_release(plane_bmag);
  gkyl_array_release(gxx_ho);
  gkyl_array_release(bmag_ho);
  if (cdim == 3) {
    gkyl_array_release(gxy_ho);
    gkyl_array_release(gyy_ho);
  }

  if (up->use_gpu) {
    defl_gyroavg_dev_buffers_new(up);
  } else {
    defl_gyroavg_host_buffers_new(up);
  }

  return up;
}

void
gkyl_deflated_fem_gyroaverage_advance(
  struct gkyl_deflated_fem_gyroaverage *up, const struct gkyl_array *fin, struct gkyl_array *fout
)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return gkyl_deflated_fem_gyroaverage_advance_cu(up, fin, fout);
  }
#endif

  int cdim = up->cdim;

  // Gyroaverage the field on each plane.
  for (int p = 0; p < up->num_planes; p++) {
    gkyl_deflate_zsurf_advance(
      up->plane_deflator[p], up->plane_zidx[p], &up->local, &up->deflated_local, fin, up->plane_fld,
      1
    );
    gkyl_fem_gyroaverage_advance(up->gyroavg[p], up->plane_fld, up->plane_fld_gyroavg[p]);
  }

  // Interpolate linearly along z between the planes bounding each cell.
  for (int cell_z = 0; cell_z < up->num_cells_z; cell_z++) {
    defl_gyroavg_inflate_planes(
      up, up->local.lower[cdim - 1] + cell_z, 1.0, 0.0, up->num_mu, up->plane_fld_gyroavg[cell_z],
      up->plane_fld_gyroavg[cell_z + 1], fout
    );
  }
}

void
gkyl_deflated_fem_gyroaverage_advance_dg(
  struct gkyl_deflated_fem_gyroaverage *up, const struct gkyl_array *fin, struct gkyl_array *fout
)
{
  // In each cell the input is f = f_lo*h_lo(z) + f_up*h_up(z) (see
  // gkyl_deflated_fem_gyroaverage_advance_adjoint), and the output is
  // G_lo(f_lo)*h_lo(z) + G_up(f_up)*h_up(z), G_lo/G_up being the gyroaverages on the lower/upper
  // z surface of the cell.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return gkyl_deflated_fem_gyroaverage_advance_dg_cu(up, fin, fout);
  }
#endif

  int cdim = up->cdim;
  for (int cell_z = 0; cell_z < up->num_cells_z; cell_z++) {
    int zidx = up->local.lower[cdim - 1] + cell_z;
    gkyl_deflate_zsurf_advance(
      up->deflator_lo, zidx, &up->local, &up->deflated_local, fin, up->plane_adj_lo, 1
    );
    gkyl_deflate_zsurf_advance(
      up->deflator_up, zidx, &up->local, &up->deflated_local, fin, up->plane_adj_up, 1
    );
    gkyl_fem_gyroaverage_advance(up->gyroavg[cell_z], up->plane_adj_lo, up->plane_dens_lo);
    gkyl_fem_gyroaverage_advance(up->gyroavg[cell_z + 1], up->plane_adj_up, up->plane_dens_up);
    defl_gyroavg_inflate_planes(
      up, zidx, 1.0, 0.0, up->num_mu, up->plane_dens_lo, up->plane_dens_up, fout
    );
  }
}

void
gkyl_deflated_fem_gyroaverage_advance_adjoint(
  struct gkyl_deflated_fem_gyroaverage *up, const struct gkyl_array *fin, struct gkyl_array *fout
)
{
  // In each cell the input is f = f_lo*h_lo(z) + f_up*h_up(z), with h_lo/h_up
  // the linear functions of z that are 1 at the lower/upper z surface and 0 at
  // the other. The adjoint of interpolating the gyroaverage along z between
  // the planes is, with M=[2 1; 1 2] proportional to the mass matrix of h,
  //   [g_lo; g_up] = M^{-1} [G_lo^dagger 0; 0 G_up^dagger] M [f_lo; f_up].
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return gkyl_deflated_fem_gyroaverage_advance_adjoint_cu(up, fin, fout);
  }
#endif

  int cdim = up->cdim;
  for (int cell_z = 0; cell_z < up->num_cells_z; cell_z++) {
    int zidx = up->local.lower[cdim - 1] + cell_z;
    gkyl_deflate_zsurf_advance(
      up->deflator_lo, zidx, &up->local, &up->deflated_local, fin, up->plane_dens_lo, up->num_mu
    );
    gkyl_deflate_zsurf_advance(
      up->deflator_up, zidx, &up->local, &up->deflated_local, fin, up->plane_dens_up, up->num_mu
    );

    // M*[f_lo; f_up].
    gkyl_array_set(up->plane_dens_buff, 2.0, up->plane_dens_lo);
    gkyl_array_accumulate(up->plane_dens_buff, 1.0, up->plane_dens_up);
    gkyl_array_scale(up->plane_dens_up, 2.0);
    gkyl_array_accumulate(up->plane_dens_up, 1.0, up->plane_dens_lo);

    gkyl_fem_gyroaverage_advance_adjoint(up->gyroavg[cell_z], up->plane_dens_buff, up->plane_adj_lo);
    gkyl_fem_gyroaverage_advance_adjoint(
      up->gyroavg[cell_z + 1], up->plane_dens_up, up->plane_adj_up
    );

    // M^{-1} = [2 -1; -1 2]/3.
    defl_gyroavg_inflate_planes(
      up, zidx, 2.0 / 3.0, -1.0 / 3.0, 1, up->plane_adj_lo, up->plane_adj_up, fout
    );
  }
}

void
gkyl_deflated_fem_gyroaverage_release(struct gkyl_deflated_fem_gyroaverage *up)
{
  if (up->use_gpu) {
    defl_gyroavg_dev_buffers_release(up);
  } else {
    defl_gyroavg_host_buffers_release(up);
  }
  for (int p = 0; p < up->num_planes; p++) {
    gkyl_fem_gyroaverage_release(up->gyroavg[p]);
  }
  gkyl_free(up->gyroavg);
  gkyl_free(up->plane_zidx);
  gkyl_free(up->plane_deflator);
  gkyl_deflate_zsurf_release(up->deflator_lo);
  gkyl_deflate_zsurf_release(up->deflator_up);
  gkyl_free(up);
}
