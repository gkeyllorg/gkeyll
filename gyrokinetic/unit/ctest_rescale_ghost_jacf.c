#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_rescale_ghost_jacf.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_util.h>
#include <acutest.h>

// Allocate array (filled with zeros).
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct ctest_ctx {
  int cdim, vdim; // Number of conf-/vel-space dimensions.
  double cells[GKYL_MAX_DIM]; // Number of cells in each direction.
  double lower[GKYL_MAX_DIM]; // Lower extent in each direction.
  double upper[GKYL_MAX_DIM]; // Upper extent in each direction.
  double n0; // Density.
  double udrift[GKYL_MAX_CDIM]; // Parallel flow speed.
  double temp; // Temperature.
  double mass; // Species mass.
};

void eval_jac_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];

  struct ctest_ctx *tctx = ctx;
  const double *lower = tctx->lower;
  const double *upper = tctx->upper;
  int cdim = tctx->cdim;
  double Lx[cdim];
  for (int d=0; d<cdim; d++) Lx[d] = upper[d] - lower[d];

  fout[0] = 1.0+0.5*sin((2.0*M_PI/Lx[0])*x);
}

void eval_distf_1x1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vx = xn[1];

  struct ctest_ctx *tctx = ctx;
  int vdim = tctx->vdim;
  double n0 = tctx->n0;
  const double *udrift = tctx->udrift;
  double temp = tctx->temp;
  double mass = tctx->mass;
  double vtsq = temp/mass;

  fout[0] = (n0/pow(sqrt(2.0*M_PI*vtsq),vdim)) * exp(-pow(vx-udrift[0],2.0)/(2.0*vtsq));
}

static void accepted_surf_inv_kernel_1x_lowerx(const double *jghost, double *out)
{
  double jgS[1];
  jgS[0] = (sqrt(2)*sqrt(3)*jghost[1]+sqrt(2)*jghost[0])/2;

  out[0] = 1/jgS[0];
}

static void accepted_results_kernel_1x1v_lowerx(const double *jskin, const double *jghost, const double *jf, double *out)
{
  double jgSinv[1];
  accepted_surf_inv_kernel_1x_lowerx(jghost, jgSinv);

  out[0] = -(2.1213203435596424*jgSinv[0]*jf[1]*jskin[1])-1.224744871391589*jf[0]*jgSinv[0]*jskin[1]+1.224744871391589*jgSinv[0]*jskin[0]*jf[1]+0.7071067811865475*jf[0]*jgSinv[0]*jskin[0];
  out[2] = -(2.1213203435596424*jgSinv[0]*jskin[1]*jf[3])+1.224744871391589*jgSinv[0]*jskin[0]*jf[3]-1.224744871391589*jgSinv[0]*jskin[1]*jf[2]+0.7071067811865475*jgSinv[0]*jskin[0]*jf[2];
  out[4] = -(2.1213203435596424*jgSinv[0]*jskin[1]*jf[5])+1.224744871391589*jgSinv[0]*jskin[0]*jf[5]-1.224744871391589*jgSinv[0]*jskin[1]*jf[4]+0.7071067811865475*jgSinv[0]*jskin[0]*jf[4];
}

static void accepted_surf_inv_kernel_1x_upperx(const double *jghost, double *out)
{
  double jgS[1];
  jgS[0] = -((sqrt(2)*sqrt(3)*jghost[1]-sqrt(2)*jghost[0])/2);

  out[0] = 1/jgS[0];
}

static void accepted_results_kernel_1x1v_upperx(const double *jskin, const double *jghost, const double *jf, double *out)
{
  double jgSinv[1];
  accepted_surf_inv_kernel_1x_upperx(jghost, jgSinv);

  out[0] = -(2.1213203435596424*jgSinv[0]*jf[1]*jskin[1])+1.224744871391589*jf[0]*jgSinv[0]*jskin[1]-1.224744871391589*jgSinv[0]*jskin[0]*jf[1]+0.7071067811865475*jf[0]*jgSinv[0]*jskin[0];
  out[2] = -(2.1213203435596424*jgSinv[0]*jskin[1]*jf[3])-1.224744871391589*jgSinv[0]*jskin[0]*jf[3]+1.224744871391589*jgSinv[0]*jskin[1]*jf[2]+0.7071067811865475*jgSinv[0]*jskin[0]*jf[2];
  out[4] = -(2.1213203435596424*jgSinv[0]*jskin[1]*jf[5])-1.224744871391589*jgSinv[0]*jskin[0]*jf[5]+1.224744871391589*jgSinv[0]*jskin[1]*jf[4]+0.7071067811865475*jgSinv[0]*jskin[0]*jf[4];
}

void test_1x1v_at_edge(bool use_gpu, int dir, enum gkyl_edge_loc edge)
{  
  int cells[] = {8, 4};
  double lower[] = {-M_PI, -6.0}, upper[] = {M_PI, 6.0};
  int cdim = 1;
  int poly_order = 1;

  int pdim = sizeof(lower)/sizeof(lower[0]);
  int vdim = pdim - cdim;

  struct ctest_ctx test_ctx = {
    .cdim = cdim,  .vdim = vdim,
    .cells = {cells[0], cells[1]},
    .lower = {lower[0], lower[1]},
    .upper = {upper[0], upper[1]},
    .n0 = 1.0,
    .udrift = {0.0},
    .temp = 2.0,
    .mass = 1.0,
  };

  // Grids.
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower, upper, cells);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);
  struct gkyl_basis basis;
  if (poly_order == 1)
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);

  // Ranges.
  int ghost_cells_conf[cdim];
  for (int d=0; d<cdim; d++) ghost_cells_conf[d] = 1;
  struct gkyl_range local_conf, local_conf_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid_conf, ghost_cells_conf, &local_conf_ext, &local_conf);

  struct gkyl_range skin_conf, ghost_conf;
  gkyl_skin_ghost_ranges(&skin_conf, &ghost_conf, dir, edge, &local_conf_ext, ghost_cells_conf);

  int ghost_cells[pdim];
  for (int d=0; d<cdim; d++) ghost_cells[d] = 1;
  for (int d=cdim; d<pdim; d++) ghost_cells[d] = 0;
  struct gkyl_range local, local_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost_cells, &local_ext, &local);

  struct gkyl_range skin, ghost;
  gkyl_skin_ghost_ranges(&skin, &ghost, dir, edge, &local_ext, ghost_cells);

  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size)
                                       : gkyl_array_acquire(distf);
  struct gkyl_array *jac = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *jac_ho = use_gpu? mkarr(false, jac->ncomp, jac->size)
                                     : gkyl_array_acquire(jac);
  struct gkyl_array *jf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *jf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size)
                                    : gkyl_array_acquire(jf);

  // Project jac onto the basis.
  struct gkyl_eval_on_nodes *proj_jac = gkyl_eval_on_nodes_new(&grid_conf, &basis_conf,
    1, eval_jac_1x, &test_ctx);
  gkyl_eval_on_nodes_advance(proj_jac, 0.0, &local_conf_ext, jac_ho);
  gkyl_eval_on_nodes_release(proj_jac);
  gkyl_array_copy(jac, jac_ho);

  // Project f onto the basis.
  struct gkyl_proj_on_basis *proj_f = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_distf_1x1v, &test_ctx);
  gkyl_proj_on_basis_advance(proj_f, 0.0, &local_ext, distf_ho);
  gkyl_proj_on_basis_release(proj_f);
  gkyl_array_copy(distf, distf_ho);

  // Multiply jac * f.
  gkyl_dg_mul_conf_phase_op_range(&basis_conf, &basis, jf, jac, distf, &local_conf_ext, &local_ext);
  // Place jf in distf as we'll need it to check results.
  gkyl_array_copy(distf_ho, jf);

  // Divide jf by j in the ghost cell, and multiply by the flipped skin cell j.
  struct gkyl_rescale_ghost_jacf* jf_rescale =
    gkyl_rescale_ghost_jacf_new(dir, edge, &basis_conf, &basis, use_gpu);

  gkyl_rescale_ghost_jacf_advance(jf_rescale,
    &skin_conf, &ghost_conf, &ghost, jac, jf);

  gkyl_rescale_ghost_jacf_release(jf_rescale);

  // Check the results.
  gkyl_array_copy(jf_ho, jf);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ghost);
  while (gkyl_range_iter_next(&iter)) {
    int cidx_skin[cdim];
    for (int d=0; d<cdim; d++) cidx_skin[d] = iter.idx[d];
    cidx_skin[dir] = edge == GKYL_LOWER_EDGE? iter.idx[dir]+1 : iter.idx[dir]-1;

    long clinidx_skin = gkyl_range_idx(&skin_conf, cidx_skin);
    long clinidx_ghost = gkyl_range_idx(&ghost_conf, iter.idx);
    long plinidx_ghost = gkyl_range_idx(&ghost, iter.idx);

    const double *jskin_c = gkyl_array_cfetch(jac_ho, clinidx_skin);
    const double *jghost_c = gkyl_array_cfetch(jac_ho, clinidx_ghost);
    const double *distf_c = gkyl_array_cfetch(distf_ho, plinidx_ghost);

    double ref_c[basis.num_basis];
    // Compute accepted results.
    for (int i=0; i<basis.num_basis; ++i)
      ref_c[i] = 0.0;

    if (edge == GKYL_LOWER_EDGE)
      accepted_results_kernel_1x1v_lowerx(jskin_c, jghost_c, distf_c, ref_c);
    else
      accepted_results_kernel_1x1v_upperx(jskin_c, jghost_c, distf_c, ref_c);

    const double *jf_c = gkyl_array_cfetch(jf_ho, plinidx_ghost);
    // for (int i=0; i<basis.num_basis; ++i) {
    //   TEST_CHECK( gkyl_compare(ref_c[i], jf_c[i], 1e-10) );
    //   TEST_MSG("Expected: %.13e | Got:%.13e | Cell:%d,%d\n", ref_c[i], jf_c[i], iter.idx[0], iter.idx[1]);
    // }
  }

  // Free memory.
  gkyl_array_release(distf_ho);
  gkyl_array_release(distf);
  gkyl_array_release(jac_ho);
  gkyl_array_release(jac);
  gkyl_array_release(jf_ho);
  gkyl_array_release(jf);
}

void eval_jac_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];

  struct ctest_ctx *tctx = ctx;
  const double *lower = tctx->lower;
  const double *upper = tctx->upper;
  int cdim = tctx->cdim;
  double Lx[cdim];
  for (int d=0; d<cdim; d++) Lx[d] = upper[d] - lower[d];

  fout[0] = 1.0+0.5*sin((2.0*M_PI/Lx[0])*x)*cos((M_PI/Lx[1])*y);;
}

void eval_distf_2x2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], vx = xn[2], vy = xn[3];

  struct ctest_ctx *tctx = ctx;
  int vdim = tctx->vdim;
  double n0 = tctx->n0;
  const double *udrift = tctx->udrift;
  double temp = tctx->temp;
  double mass = tctx->mass;
  double vtsq = temp/mass;

  fout[0] = (n0/pow(sqrt(2.0*M_PI*vtsq),vdim))
    * exp(-(pow(vx-udrift[0],2.0)+pow(vy-udrift[1],2.0))/(2.0*vtsq));
}

static void accepted_surf_inv_kernel_2x_lowerx(const double *jghost, double *out)
{
  double jgS[2];
  jgS[0] = (sqrt(3)*jghost[1]+jghost[0])/sqrt(2);
  jgS[1] = (sqrt(3)*jghost[3]+jghost[2])/sqrt(2);

  out[0] = -((2*jgS[0])/(pow(jgS[1],2)-pow(jgS[0],2)));
  out[1] =   (2*jgS[1])/(pow(jgS[1],2)-pow(jgS[0],2));
}

static void accepted_results_kernel_2x2v_lowerx(const double *jskin, const double *jghost, const double *jf, double *out)
{
  double jgSinv[2];
  accepted_surf_inv_kernel_2x_lowerx(jghost, jgSinv);

  out[0] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[5])+0.6123724356957944*jgSinv[0]*jskin[2]*jf[5]-1.060660171779821*jgSinv[1]*jskin[1]*jf[5]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[5]-0.6123724356957944*jgSinv[0]*jf[2]*jskin[3]-1.060660171779821*jf[1]*jgSinv[1]*jskin[3]-0.6123724356957944*jf[0]*jgSinv[1]*jskin[3]+0.3535533905932737*jgSinv[0]*jf[2]*jskin[2]+0.6123724356957944*jf[1]*jgSinv[1]*jskin[2]+0.3535533905932737*jf[0]*jgSinv[1]*jskin[2]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[2]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[2]-1.060660171779821*jgSinv[0]*jf[1]*jskin[1]-0.6123724356957944*jf[0]*jgSinv[0]*jskin[1]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[1]+0.3535533905932737*jf[0]*jgSinv[0]*jskin[0]; 
  out[2] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[5])+0.6123724356957944*jgSinv[1]*jskin[2]*jf[5]-1.060660171779821*jgSinv[0]*jskin[1]*jf[5]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[5]-0.6123724356957944*jgSinv[1]*jf[2]*jskin[3]-1.060660171779821*jgSinv[0]*jf[1]*jskin[3]-0.6123724356957944*jf[0]*jgSinv[0]*jskin[3]+0.3535533905932737*jgSinv[1]*jf[2]*jskin[2]+0.6123724356957944*jgSinv[0]*jf[1]*jskin[2]+0.3535533905932737*jf[0]*jgSinv[0]*jskin[2]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[2]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[2]-1.060660171779821*jf[1]*jgSinv[1]*jskin[1]-0.6123724356957944*jf[0]*jgSinv[1]*jskin[1]+0.6123724356957944*jskin[0]*jf[1]*jgSinv[1]+0.3535533905932737*jf[0]*jskin[0]*jgSinv[1]; 
  out[3] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[11])+0.6123724356957944*jgSinv[0]*jskin[2]*jf[11]-1.060660171779821*jgSinv[1]*jskin[1]*jf[11]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[11]-0.6123724356957944*jgSinv[0]*jskin[3]*jf[7]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[7]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[7]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[7]-1.060660171779821*jgSinv[1]*jskin[3]*jf[6]+0.6123724356957944*jgSinv[1]*jskin[2]*jf[6]-1.060660171779821*jgSinv[0]*jskin[1]*jf[6]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[6]-0.6123724356957944*jgSinv[1]*jf[3]*jskin[3]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[3]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[3]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[3]; 
  out[4] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[12])+0.6123724356957944*jgSinv[0]*jskin[2]*jf[12]-1.060660171779821*jgSinv[1]*jskin[1]*jf[12]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[12]-0.6123724356957944*jgSinv[0]*jskin[3]*jf[9]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[9]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[9]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[9]-1.060660171779821*jgSinv[1]*jskin[3]*jf[8]+0.6123724356957944*jgSinv[1]*jskin[2]*jf[8]-1.060660171779821*jgSinv[0]*jskin[1]*jf[8]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[8]-0.6123724356957944*jgSinv[1]*jskin[3]*jf[4]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[4]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[4]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[4]; 
  out[7] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[11])+0.6123724356957944*jgSinv[1]*jskin[2]*jf[11]-1.060660171779821*jgSinv[0]*jskin[1]*jf[11]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[11]-0.6123724356957944*jgSinv[1]*jskin[3]*jf[7]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[7]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[7]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[7]-1.060660171779821*jgSinv[0]*jskin[3]*jf[6]+0.6123724356957944*jgSinv[0]*jskin[2]*jf[6]-1.060660171779821*jgSinv[1]*jskin[1]*jf[6]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[6]-0.6123724356957944*jgSinv[0]*jf[3]*jskin[3]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[3]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[3]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[3]; 
  out[9] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[12])+0.6123724356957944*jgSinv[1]*jskin[2]*jf[12]-1.060660171779821*jgSinv[0]*jskin[1]*jf[12]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[12]-0.6123724356957944*jgSinv[1]*jskin[3]*jf[9]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[9]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[9]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[9]-1.060660171779821*jgSinv[0]*jskin[3]*jf[8]+0.6123724356957944*jgSinv[0]*jskin[2]*jf[8]-1.060660171779821*jgSinv[1]*jskin[1]*jf[8]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[8]-0.6123724356957944*jgSinv[0]*jskin[3]*jf[4]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[4]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[4]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[4]; 
  out[10] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[15])+0.6123724356957944*jgSinv[0]*jskin[2]*jf[15]-1.060660171779821*jgSinv[1]*jskin[1]*jf[15]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[15]-0.6123724356957944*jgSinv[0]*jskin[3]*jf[14]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[14]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[14]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[14]-1.060660171779821*jgSinv[1]*jskin[3]*jf[13]+0.6123724356957944*jgSinv[1]*jskin[2]*jf[13]-1.060660171779821*jgSinv[0]*jskin[1]*jf[13]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[13]-0.6123724356957944*jgSinv[1]*jskin[3]*jf[10]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[10]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[10]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[10]; 
  out[14] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[15])+0.6123724356957944*jgSinv[1]*jskin[2]*jf[15]-1.060660171779821*jgSinv[0]*jskin[1]*jf[15]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[15]-0.6123724356957944*jgSinv[1]*jskin[3]*jf[14]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[14]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[14]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[14]-1.060660171779821*jgSinv[0]*jskin[3]*jf[13]+0.6123724356957944*jgSinv[0]*jskin[2]*jf[13]-1.060660171779821*jgSinv[1]*jskin[1]*jf[13]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[13]-0.6123724356957944*jgSinv[0]*jskin[3]*jf[10]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[10]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[10]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[10]; 
  out[16] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[20])+0.6123724356957944*jgSinv[0]*jskin[2]*jf[20]-1.060660171779821*jgSinv[1]*jskin[1]*jf[20]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[20]-0.6123724356957944*jgSinv[0]*jskin[3]*jf[18]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[18]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[18]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[18]-1.060660171779821*jgSinv[1]*jskin[3]*jf[17]+0.6123724356957944*jgSinv[1]*jskin[2]*jf[17]-1.060660171779821*jgSinv[0]*jskin[1]*jf[17]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[17]-0.6123724356957944*jgSinv[1]*jskin[3]*jf[16]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[16]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[16]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[16]; 
  out[18] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[20])+0.6123724356957944*jgSinv[1]*jskin[2]*jf[20]-1.060660171779821*jgSinv[0]*jskin[1]*jf[20]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[20]-0.6123724356957944*jgSinv[1]*jskin[3]*jf[18]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[18]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[18]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[18]-1.060660171779821*jgSinv[0]*jskin[3]*jf[17]+0.6123724356957944*jgSinv[0]*jskin[2]*jf[17]-1.060660171779821*jgSinv[1]*jskin[1]*jf[17]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[17]-0.6123724356957944*jgSinv[0]*jskin[3]*jf[16]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[16]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[16]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[16]; 
  out[19] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[23])+0.6123724356957944*jgSinv[0]*jskin[2]*jf[23]-1.060660171779821*jgSinv[1]*jskin[1]*jf[23]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[23]-0.6123724356957944*jgSinv[0]*jskin[3]*jf[22]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[22]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[22]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[22]-1.060660171779821*jgSinv[1]*jskin[3]*jf[21]+0.6123724356957944*jgSinv[1]*jskin[2]*jf[21]-1.060660171779821*jgSinv[0]*jskin[1]*jf[21]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[21]-0.6123724356957944*jgSinv[1]*jskin[3]*jf[19]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[19]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[19]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[19]; 
  out[22] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[23])+0.6123724356957944*jgSinv[1]*jskin[2]*jf[23]-1.060660171779821*jgSinv[0]*jskin[1]*jf[23]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[23]-0.6123724356957944*jgSinv[1]*jskin[3]*jf[22]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[22]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[22]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[22]-1.060660171779821*jgSinv[0]*jskin[3]*jf[21]+0.6123724356957944*jgSinv[0]*jskin[2]*jf[21]-1.060660171779821*jgSinv[1]*jskin[1]*jf[21]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[21]-0.6123724356957944*jgSinv[0]*jskin[3]*jf[19]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[19]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[19]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[19]; 
}

static void accepted_surf_inv_kernel_2x_upperx(const double *jghost, double *out)
{
  double jgS[2];
  jgS[0] = -((sqrt(3)*jghost[1]-jghost[0])/sqrt(2));
  jgS[1] = -((sqrt(3)*jghost[3]-jghost[2])/sqrt(2));
  
  out[0] = -((2*jgS[0])/(pow(jgS[1],2)-pow(jgS[0],2)));
  out[1] =   (2*jgS[1])/(pow(jgS[1],2)-pow(jgS[0],2));
}

static void accepted_results_kernel_2x2v_upperx(const double *jskin, const double *jghost, const double *jf, double *out)
{
  double jgSinv[2];
  accepted_surf_inv_kernel_2x_upperx(jghost, jgSinv);

  out[0] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[5])-0.6123724356957944*jgSinv[0]*jskin[2]*jf[5]-1.060660171779821*jgSinv[1]*jskin[1]*jf[5]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[5]+0.6123724356957944*jgSinv[0]*jf[2]*jskin[3]-1.060660171779821*jf[1]*jgSinv[1]*jskin[3]+0.6123724356957944*jf[0]*jgSinv[1]*jskin[3]+0.3535533905932737*jgSinv[0]*jf[2]*jskin[2]-0.6123724356957944*jf[1]*jgSinv[1]*jskin[2]+0.3535533905932737*jf[0]*jgSinv[1]*jskin[2]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[2]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[2]-1.060660171779821*jgSinv[0]*jf[1]*jskin[1]+0.6123724356957944*jf[0]*jgSinv[0]*jskin[1]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[1]+0.3535533905932737*jf[0]*jgSinv[0]*jskin[0]; 
  out[2] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[5])-0.6123724356957944*jgSinv[1]*jskin[2]*jf[5]-1.060660171779821*jgSinv[0]*jskin[1]*jf[5]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[5]+0.6123724356957944*jgSinv[1]*jf[2]*jskin[3]-1.060660171779821*jgSinv[0]*jf[1]*jskin[3]+0.6123724356957944*jf[0]*jgSinv[0]*jskin[3]+0.3535533905932737*jgSinv[1]*jf[2]*jskin[2]-0.6123724356957944*jgSinv[0]*jf[1]*jskin[2]+0.3535533905932737*jf[0]*jgSinv[0]*jskin[2]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[2]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[2]-1.060660171779821*jf[1]*jgSinv[1]*jskin[1]+0.6123724356957944*jf[0]*jgSinv[1]*jskin[1]-0.6123724356957944*jskin[0]*jf[1]*jgSinv[1]+0.3535533905932737*jf[0]*jskin[0]*jgSinv[1]; 
  out[3] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[11])-0.6123724356957944*jgSinv[0]*jskin[2]*jf[11]-1.060660171779821*jgSinv[1]*jskin[1]*jf[11]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[11]+0.6123724356957944*jgSinv[0]*jskin[3]*jf[7]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[7]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[7]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[7]-1.060660171779821*jgSinv[1]*jskin[3]*jf[6]-0.6123724356957944*jgSinv[1]*jskin[2]*jf[6]-1.060660171779821*jgSinv[0]*jskin[1]*jf[6]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[6]+0.6123724356957944*jgSinv[1]*jf[3]*jskin[3]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[3]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[3]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[3]; 
  out[4] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[12])-0.6123724356957944*jgSinv[0]*jskin[2]*jf[12]-1.060660171779821*jgSinv[1]*jskin[1]*jf[12]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[12]+0.6123724356957944*jgSinv[0]*jskin[3]*jf[9]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[9]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[9]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[9]-1.060660171779821*jgSinv[1]*jskin[3]*jf[8]-0.6123724356957944*jgSinv[1]*jskin[2]*jf[8]-1.060660171779821*jgSinv[0]*jskin[1]*jf[8]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[8]+0.6123724356957944*jgSinv[1]*jskin[3]*jf[4]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[4]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[4]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[4]; 
  out[7] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[11])-0.6123724356957944*jgSinv[1]*jskin[2]*jf[11]-1.060660171779821*jgSinv[0]*jskin[1]*jf[11]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[11]+0.6123724356957944*jgSinv[1]*jskin[3]*jf[7]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[7]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[7]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[7]-1.060660171779821*jgSinv[0]*jskin[3]*jf[6]-0.6123724356957944*jgSinv[0]*jskin[2]*jf[6]-1.060660171779821*jgSinv[1]*jskin[1]*jf[6]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[6]+0.6123724356957944*jgSinv[0]*jf[3]*jskin[3]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[3]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[3]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[3]; 
  out[9] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[12])-0.6123724356957944*jgSinv[1]*jskin[2]*jf[12]-1.060660171779821*jgSinv[0]*jskin[1]*jf[12]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[12]+0.6123724356957944*jgSinv[1]*jskin[3]*jf[9]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[9]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[9]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[9]-1.060660171779821*jgSinv[0]*jskin[3]*jf[8]-0.6123724356957944*jgSinv[0]*jskin[2]*jf[8]-1.060660171779821*jgSinv[1]*jskin[1]*jf[8]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[8]+0.6123724356957944*jgSinv[0]*jskin[3]*jf[4]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[4]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[4]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[4]; 
  out[10] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[15])-0.6123724356957944*jgSinv[0]*jskin[2]*jf[15]-1.060660171779821*jgSinv[1]*jskin[1]*jf[15]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[15]+0.6123724356957944*jgSinv[0]*jskin[3]*jf[14]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[14]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[14]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[14]-1.060660171779821*jgSinv[1]*jskin[3]*jf[13]-0.6123724356957944*jgSinv[1]*jskin[2]*jf[13]-1.060660171779821*jgSinv[0]*jskin[1]*jf[13]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[13]+0.6123724356957944*jgSinv[1]*jskin[3]*jf[10]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[10]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[10]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[10]; 
  out[14] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[15])-0.6123724356957944*jgSinv[1]*jskin[2]*jf[15]-1.060660171779821*jgSinv[0]*jskin[1]*jf[15]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[15]+0.6123724356957944*jgSinv[1]*jskin[3]*jf[14]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[14]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[14]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[14]-1.060660171779821*jgSinv[0]*jskin[3]*jf[13]-0.6123724356957944*jgSinv[0]*jskin[2]*jf[13]-1.060660171779821*jgSinv[1]*jskin[1]*jf[13]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[13]+0.6123724356957944*jgSinv[0]*jskin[3]*jf[10]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[10]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[10]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[10]; 
  out[16] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[20])-0.6123724356957944*jgSinv[0]*jskin[2]*jf[20]-1.060660171779821*jgSinv[1]*jskin[1]*jf[20]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[20]+0.6123724356957944*jgSinv[0]*jskin[3]*jf[18]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[18]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[18]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[18]-1.060660171779821*jgSinv[1]*jskin[3]*jf[17]-0.6123724356957944*jgSinv[1]*jskin[2]*jf[17]-1.060660171779821*jgSinv[0]*jskin[1]*jf[17]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[17]+0.6123724356957944*jgSinv[1]*jskin[3]*jf[16]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[16]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[16]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[16]; 
  out[18] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[20])-0.6123724356957944*jgSinv[1]*jskin[2]*jf[20]-1.060660171779821*jgSinv[0]*jskin[1]*jf[20]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[20]+0.6123724356957944*jgSinv[1]*jskin[3]*jf[18]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[18]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[18]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[18]-1.060660171779821*jgSinv[0]*jskin[3]*jf[17]-0.6123724356957944*jgSinv[0]*jskin[2]*jf[17]-1.060660171779821*jgSinv[1]*jskin[1]*jf[17]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[17]+0.6123724356957944*jgSinv[0]*jskin[3]*jf[16]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[16]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[16]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[16]; 
  out[19] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[23])-0.6123724356957944*jgSinv[0]*jskin[2]*jf[23]-1.060660171779821*jgSinv[1]*jskin[1]*jf[23]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[23]+0.6123724356957944*jgSinv[0]*jskin[3]*jf[22]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[22]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[22]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[22]-1.060660171779821*jgSinv[1]*jskin[3]*jf[21]-0.6123724356957944*jgSinv[1]*jskin[2]*jf[21]-1.060660171779821*jgSinv[0]*jskin[1]*jf[21]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[21]+0.6123724356957944*jgSinv[1]*jskin[3]*jf[19]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[19]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[19]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[19]; 
  out[22] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[23])-0.6123724356957944*jgSinv[1]*jskin[2]*jf[23]-1.060660171779821*jgSinv[0]*jskin[1]*jf[23]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[23]+0.6123724356957944*jgSinv[1]*jskin[3]*jf[22]+0.3535533905932737*jgSinv[1]*jskin[2]*jf[22]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[22]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[22]-1.060660171779821*jgSinv[0]*jskin[3]*jf[21]-0.6123724356957944*jgSinv[0]*jskin[2]*jf[21]-1.060660171779821*jgSinv[1]*jskin[1]*jf[21]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[21]+0.6123724356957944*jgSinv[0]*jskin[3]*jf[19]+0.3535533905932737*jgSinv[0]*jskin[2]*jf[19]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[19]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[19]; 
}

static void accepted_surf_inv_kernel_2x_lowery(const double *jghost, double *out)
{
  double jgS[2];
  jgS[0] = (sqrt(3)*jghost[2]+jghost[0])/sqrt(2);
  jgS[1] = (sqrt(3)*jghost[3]+jghost[1])/sqrt(2);
  
  out[0] = -((2*jgS[0])/(pow(jgS[1],2)-pow(jgS[0],2)));
  out[1] =   (2*jgS[1])/(pow(jgS[1],2)-pow(jgS[0],2));
}

static void accepted_results_kernel_2x2v_lowery(const double *jskin, const double *jghost, const double *jf, double *out)
{
  double jgSinv[2];
  accepted_surf_inv_kernel_2x_lowery(jghost, jgSinv);

  out[0] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[5])-1.060660171779821*jgSinv[1]*jskin[2]*jf[5]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[5]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[5]-1.060660171779821*jgSinv[1]*jf[2]*jskin[3]-0.6123724356957944*jf[0]*jgSinv[1]*jskin[3]-0.6123724356957944*jgSinv[0]*jf[1]*jskin[3]-1.060660171779821*jgSinv[0]*jf[2]*jskin[2]-0.6123724356957944*jf[1]*jgSinv[1]*jskin[2]-0.6123724356957944*jf[0]*jgSinv[0]*jskin[2]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[2]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[2]+0.3535533905932737*jf[0]*jgSinv[1]*jskin[1]+0.3535533905932737*jgSinv[0]*jf[1]*jskin[1]+0.3535533905932737*jskin[0]*jf[1]*jgSinv[1]+0.3535533905932737*jf[0]*jgSinv[0]*jskin[0]; 
  out[1] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[5])-1.060660171779821*jgSinv[0]*jskin[2]*jf[5]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[5]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[5]-1.060660171779821*jgSinv[0]*jf[2]*jskin[3]-0.6123724356957944*jf[1]*jgSinv[1]*jskin[3]-0.6123724356957944*jf[0]*jgSinv[0]*jskin[3]-1.060660171779821*jgSinv[1]*jf[2]*jskin[2]-0.6123724356957944*jf[0]*jgSinv[1]*jskin[2]-0.6123724356957944*jgSinv[0]*jf[1]*jskin[2]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[2]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[2]+0.3535533905932737*jf[1]*jgSinv[1]*jskin[1]+0.3535533905932737*jf[0]*jgSinv[0]*jskin[1]+0.3535533905932737*jf[0]*jskin[0]*jgSinv[1]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[1]; 
  out[3] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[11])-1.060660171779821*jgSinv[1]*jskin[2]*jf[11]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[11]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[11]-1.060660171779821*jgSinv[1]*jskin[3]*jf[7]-1.060660171779821*jgSinv[0]*jskin[2]*jf[7]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[7]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[7]-0.6123724356957944*jgSinv[0]*jskin[3]*jf[6]-0.6123724356957944*jgSinv[1]*jskin[2]*jf[6]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[6]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[6]-0.6123724356957944*jgSinv[1]*jf[3]*jskin[3]-0.6123724356957944*jgSinv[0]*jskin[2]*jf[3]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[3]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[3]; 
  out[4] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[12])-1.060660171779821*jgSinv[1]*jskin[2]*jf[12]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[12]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[12]-1.060660171779821*jgSinv[1]*jskin[3]*jf[9]-1.060660171779821*jgSinv[0]*jskin[2]*jf[9]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[9]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[9]-0.6123724356957944*jgSinv[0]*jskin[3]*jf[8]-0.6123724356957944*jgSinv[1]*jskin[2]*jf[8]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[8]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[8]-0.6123724356957944*jgSinv[1]*jskin[3]*jf[4]-0.6123724356957944*jgSinv[0]*jskin[2]*jf[4]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[4]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[4]; 
  out[6] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[11])-1.060660171779821*jgSinv[0]*jskin[2]*jf[11]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[11]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[11]-1.060660171779821*jgSinv[0]*jskin[3]*jf[7]-1.060660171779821*jgSinv[1]*jskin[2]*jf[7]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[7]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[7]-0.6123724356957944*jgSinv[1]*jskin[3]*jf[6]-0.6123724356957944*jgSinv[0]*jskin[2]*jf[6]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[6]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[6]-0.6123724356957944*jgSinv[0]*jf[3]*jskin[3]-0.6123724356957944*jgSinv[1]*jskin[2]*jf[3]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[3]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[3]; 
  out[8] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[12])-1.060660171779821*jgSinv[0]*jskin[2]*jf[12]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[12]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[12]-1.060660171779821*jgSinv[0]*jskin[3]*jf[9]-1.060660171779821*jgSinv[1]*jskin[2]*jf[9]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[9]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[9]-0.6123724356957944*jgSinv[1]*jskin[3]*jf[8]-0.6123724356957944*jgSinv[0]*jskin[2]*jf[8]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[8]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[8]-0.6123724356957944*jgSinv[0]*jskin[3]*jf[4]-0.6123724356957944*jgSinv[1]*jskin[2]*jf[4]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[4]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[4]; 
  out[10] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[15])-1.060660171779821*jgSinv[1]*jskin[2]*jf[15]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[15]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[15]-1.060660171779821*jgSinv[1]*jskin[3]*jf[14]-1.060660171779821*jgSinv[0]*jskin[2]*jf[14]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[14]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[14]-0.6123724356957944*jgSinv[0]*jskin[3]*jf[13]-0.6123724356957944*jgSinv[1]*jskin[2]*jf[13]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[13]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[13]-0.6123724356957944*jgSinv[1]*jskin[3]*jf[10]-0.6123724356957944*jgSinv[0]*jskin[2]*jf[10]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[10]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[10]; 
  out[13] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[15])-1.060660171779821*jgSinv[0]*jskin[2]*jf[15]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[15]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[15]-1.060660171779821*jgSinv[0]*jskin[3]*jf[14]-1.060660171779821*jgSinv[1]*jskin[2]*jf[14]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[14]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[14]-0.6123724356957944*jgSinv[1]*jskin[3]*jf[13]-0.6123724356957944*jgSinv[0]*jskin[2]*jf[13]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[13]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[13]-0.6123724356957944*jgSinv[0]*jskin[3]*jf[10]-0.6123724356957944*jgSinv[1]*jskin[2]*jf[10]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[10]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[10]; 
  out[16] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[20])-1.060660171779821*jgSinv[1]*jskin[2]*jf[20]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[20]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[20]-1.060660171779821*jgSinv[1]*jskin[3]*jf[18]-1.060660171779821*jgSinv[0]*jskin[2]*jf[18]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[18]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[18]-0.6123724356957944*jgSinv[0]*jskin[3]*jf[17]-0.6123724356957944*jgSinv[1]*jskin[2]*jf[17]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[17]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[17]-0.6123724356957944*jgSinv[1]*jskin[3]*jf[16]-0.6123724356957944*jgSinv[0]*jskin[2]*jf[16]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[16]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[16]; 
  out[17] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[20])-1.060660171779821*jgSinv[0]*jskin[2]*jf[20]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[20]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[20]-1.060660171779821*jgSinv[0]*jskin[3]*jf[18]-1.060660171779821*jgSinv[1]*jskin[2]*jf[18]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[18]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[18]-0.6123724356957944*jgSinv[1]*jskin[3]*jf[17]-0.6123724356957944*jgSinv[0]*jskin[2]*jf[17]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[17]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[17]-0.6123724356957944*jgSinv[0]*jskin[3]*jf[16]-0.6123724356957944*jgSinv[1]*jskin[2]*jf[16]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[16]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[16]; 
  out[19] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[23])-1.060660171779821*jgSinv[1]*jskin[2]*jf[23]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[23]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[23]-1.060660171779821*jgSinv[1]*jskin[3]*jf[22]-1.060660171779821*jgSinv[0]*jskin[2]*jf[22]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[22]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[22]-0.6123724356957944*jgSinv[0]*jskin[3]*jf[21]-0.6123724356957944*jgSinv[1]*jskin[2]*jf[21]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[21]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[21]-0.6123724356957944*jgSinv[1]*jskin[3]*jf[19]-0.6123724356957944*jgSinv[0]*jskin[2]*jf[19]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[19]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[19]; 
  out[21] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[23])-1.060660171779821*jgSinv[0]*jskin[2]*jf[23]+0.6123724356957944*jgSinv[1]*jskin[1]*jf[23]+0.6123724356957944*jgSinv[0]*jskin[0]*jf[23]-1.060660171779821*jgSinv[0]*jskin[3]*jf[22]-1.060660171779821*jgSinv[1]*jskin[2]*jf[22]+0.6123724356957944*jgSinv[0]*jskin[1]*jf[22]+0.6123724356957944*jskin[0]*jgSinv[1]*jf[22]-0.6123724356957944*jgSinv[1]*jskin[3]*jf[21]-0.6123724356957944*jgSinv[0]*jskin[2]*jf[21]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[21]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[21]-0.6123724356957944*jgSinv[0]*jskin[3]*jf[19]-0.6123724356957944*jgSinv[1]*jskin[2]*jf[19]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[19]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[19]; 
}

static void accepted_surf_inv_kernel_2x_uppery(const double *jghost, double *out)
{
  double jgS[2];
  jgS[0] = -((sqrt(3)*jghost[2]-jghost[0])/sqrt(2));
  jgS[1] = -((sqrt(3)*jghost[3]-jghost[1])/sqrt(2));
  
  out[0] = -((2*jgS[0])/(pow(jgS[1],2)-pow(jgS[0],2)));
  out[1] =   (2*jgS[1])/(pow(jgS[1],2)-pow(jgS[0],2));
}

static void accepted_results_kernel_2x2v_uppery(const double *jskin, const double *jghost, const double *jf, double *out)
{
  double jgSinv[2];
  accepted_surf_inv_kernel_2x_uppery(jghost, jgSinv);

  out[0] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[5])-1.060660171779821*jgSinv[1]*jskin[2]*jf[5]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[5]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[5]-1.060660171779821*jgSinv[1]*jf[2]*jskin[3]+0.6123724356957944*jf[0]*jgSinv[1]*jskin[3]+0.6123724356957944*jgSinv[0]*jf[1]*jskin[3]-1.060660171779821*jgSinv[0]*jf[2]*jskin[2]+0.6123724356957944*jf[1]*jgSinv[1]*jskin[2]+0.6123724356957944*jf[0]*jgSinv[0]*jskin[2]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[2]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[2]+0.3535533905932737*jf[0]*jgSinv[1]*jskin[1]+0.3535533905932737*jgSinv[0]*jf[1]*jskin[1]+0.3535533905932737*jskin[0]*jf[1]*jgSinv[1]+0.3535533905932737*jf[0]*jgSinv[0]*jskin[0]; 
  out[1] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[5])-1.060660171779821*jgSinv[0]*jskin[2]*jf[5]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[5]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[5]-1.060660171779821*jgSinv[0]*jf[2]*jskin[3]+0.6123724356957944*jf[1]*jgSinv[1]*jskin[3]+0.6123724356957944*jf[0]*jgSinv[0]*jskin[3]-1.060660171779821*jgSinv[1]*jf[2]*jskin[2]+0.6123724356957944*jf[0]*jgSinv[1]*jskin[2]+0.6123724356957944*jgSinv[0]*jf[1]*jskin[2]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[2]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[2]+0.3535533905932737*jf[1]*jgSinv[1]*jskin[1]+0.3535533905932737*jf[0]*jgSinv[0]*jskin[1]+0.3535533905932737*jf[0]*jskin[0]*jgSinv[1]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[1]; 
  out[3] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[11])-1.060660171779821*jgSinv[1]*jskin[2]*jf[11]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[11]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[11]-1.060660171779821*jgSinv[1]*jskin[3]*jf[7]-1.060660171779821*jgSinv[0]*jskin[2]*jf[7]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[7]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[7]+0.6123724356957944*jgSinv[0]*jskin[3]*jf[6]+0.6123724356957944*jgSinv[1]*jskin[2]*jf[6]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[6]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[6]+0.6123724356957944*jgSinv[1]*jf[3]*jskin[3]+0.6123724356957944*jgSinv[0]*jskin[2]*jf[3]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[3]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[3]; 
  out[4] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[12])-1.060660171779821*jgSinv[1]*jskin[2]*jf[12]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[12]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[12]-1.060660171779821*jgSinv[1]*jskin[3]*jf[9]-1.060660171779821*jgSinv[0]*jskin[2]*jf[9]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[9]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[9]+0.6123724356957944*jgSinv[0]*jskin[3]*jf[8]+0.6123724356957944*jgSinv[1]*jskin[2]*jf[8]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[8]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[8]+0.6123724356957944*jgSinv[1]*jskin[3]*jf[4]+0.6123724356957944*jgSinv[0]*jskin[2]*jf[4]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[4]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[4]; 
  out[6] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[11])-1.060660171779821*jgSinv[0]*jskin[2]*jf[11]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[11]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[11]-1.060660171779821*jgSinv[0]*jskin[3]*jf[7]-1.060660171779821*jgSinv[1]*jskin[2]*jf[7]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[7]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[7]+0.6123724356957944*jgSinv[1]*jskin[3]*jf[6]+0.6123724356957944*jgSinv[0]*jskin[2]*jf[6]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[6]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[6]+0.6123724356957944*jgSinv[0]*jf[3]*jskin[3]+0.6123724356957944*jgSinv[1]*jskin[2]*jf[3]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[3]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[3]; 
  out[8] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[12])-1.060660171779821*jgSinv[0]*jskin[2]*jf[12]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[12]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[12]-1.060660171779821*jgSinv[0]*jskin[3]*jf[9]-1.060660171779821*jgSinv[1]*jskin[2]*jf[9]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[9]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[9]+0.6123724356957944*jgSinv[1]*jskin[3]*jf[8]+0.6123724356957944*jgSinv[0]*jskin[2]*jf[8]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[8]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[8]+0.6123724356957944*jgSinv[0]*jskin[3]*jf[4]+0.6123724356957944*jgSinv[1]*jskin[2]*jf[4]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[4]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[4]; 
  out[10] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[15])-1.060660171779821*jgSinv[1]*jskin[2]*jf[15]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[15]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[15]-1.060660171779821*jgSinv[1]*jskin[3]*jf[14]-1.060660171779821*jgSinv[0]*jskin[2]*jf[14]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[14]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[14]+0.6123724356957944*jgSinv[0]*jskin[3]*jf[13]+0.6123724356957944*jgSinv[1]*jskin[2]*jf[13]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[13]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[13]+0.6123724356957944*jgSinv[1]*jskin[3]*jf[10]+0.6123724356957944*jgSinv[0]*jskin[2]*jf[10]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[10]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[10]; 
  out[13] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[15])-1.060660171779821*jgSinv[0]*jskin[2]*jf[15]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[15]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[15]-1.060660171779821*jgSinv[0]*jskin[3]*jf[14]-1.060660171779821*jgSinv[1]*jskin[2]*jf[14]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[14]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[14]+0.6123724356957944*jgSinv[1]*jskin[3]*jf[13]+0.6123724356957944*jgSinv[0]*jskin[2]*jf[13]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[13]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[13]+0.6123724356957944*jgSinv[0]*jskin[3]*jf[10]+0.6123724356957944*jgSinv[1]*jskin[2]*jf[10]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[10]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[10]; 
  out[16] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[20])-1.060660171779821*jgSinv[1]*jskin[2]*jf[20]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[20]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[20]-1.060660171779821*jgSinv[1]*jskin[3]*jf[18]-1.060660171779821*jgSinv[0]*jskin[2]*jf[18]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[18]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[18]+0.6123724356957944*jgSinv[0]*jskin[3]*jf[17]+0.6123724356957944*jgSinv[1]*jskin[2]*jf[17]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[17]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[17]+0.6123724356957944*jgSinv[1]*jskin[3]*jf[16]+0.6123724356957944*jgSinv[0]*jskin[2]*jf[16]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[16]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[16]; 
  out[17] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[20])-1.060660171779821*jgSinv[0]*jskin[2]*jf[20]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[20]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[20]-1.060660171779821*jgSinv[0]*jskin[3]*jf[18]-1.060660171779821*jgSinv[1]*jskin[2]*jf[18]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[18]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[18]+0.6123724356957944*jgSinv[1]*jskin[3]*jf[17]+0.6123724356957944*jgSinv[0]*jskin[2]*jf[17]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[17]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[17]+0.6123724356957944*jgSinv[0]*jskin[3]*jf[16]+0.6123724356957944*jgSinv[1]*jskin[2]*jf[16]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[16]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[16]; 
  out[19] = -(1.060660171779821*jgSinv[0]*jskin[3]*jf[23])-1.060660171779821*jgSinv[1]*jskin[2]*jf[23]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[23]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[23]-1.060660171779821*jgSinv[1]*jskin[3]*jf[22]-1.060660171779821*jgSinv[0]*jskin[2]*jf[22]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[22]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[22]+0.6123724356957944*jgSinv[0]*jskin[3]*jf[21]+0.6123724356957944*jgSinv[1]*jskin[2]*jf[21]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[21]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[21]+0.6123724356957944*jgSinv[1]*jskin[3]*jf[19]+0.6123724356957944*jgSinv[0]*jskin[2]*jf[19]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[19]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[19]; 
  out[21] = -(1.060660171779821*jgSinv[1]*jskin[3]*jf[23])-1.060660171779821*jgSinv[0]*jskin[2]*jf[23]-0.6123724356957944*jgSinv[1]*jskin[1]*jf[23]-0.6123724356957944*jgSinv[0]*jskin[0]*jf[23]-1.060660171779821*jgSinv[0]*jskin[3]*jf[22]-1.060660171779821*jgSinv[1]*jskin[2]*jf[22]-0.6123724356957944*jgSinv[0]*jskin[1]*jf[22]-0.6123724356957944*jskin[0]*jgSinv[1]*jf[22]+0.6123724356957944*jgSinv[1]*jskin[3]*jf[21]+0.6123724356957944*jgSinv[0]*jskin[2]*jf[21]+0.3535533905932737*jgSinv[1]*jskin[1]*jf[21]+0.3535533905932737*jgSinv[0]*jskin[0]*jf[21]+0.6123724356957944*jgSinv[0]*jskin[3]*jf[19]+0.6123724356957944*jgSinv[1]*jskin[2]*jf[19]+0.3535533905932737*jgSinv[0]*jskin[1]*jf[19]+0.3535533905932737*jskin[0]*jgSinv[1]*jf[19]; 
}

void test_2x2v_at_edge(bool use_gpu, int dir, enum gkyl_edge_loc edge)
{  
  int cells[] = {8, 6, 8, 4};
  double lower[] = {-M_PI, -M_PI, -6.0, -6.0}, upper[] = {M_PI, M_PI, 6.0, 6.0};
  int cdim = 2;
  int poly_order = 1;

  int pdim = sizeof(lower)/sizeof(lower[0]);
  int vdim = pdim - cdim;

  struct ctest_ctx test_ctx = {
    .cdim = cdim,  .vdim = vdim,
    .cells = {cells[0], cells[1]},
    .lower = {lower[0], lower[1]},
    .upper = {upper[0], upper[1]},
    .n0 = 1.0,
    .udrift = {0.0, 0.0},
    .temp = 2.0,
    .mass = 1.0,
  };

  // Grids.
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower, upper, cells);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);
  struct gkyl_basis basis;
  if (poly_order == 1)
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);

  // Ranges.
  int ghost_cells_conf[cdim];
  for (int d=0; d<cdim; d++) ghost_cells_conf[d] = 1;
  struct gkyl_range local_conf, local_conf_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid_conf, ghost_cells_conf, &local_conf_ext, &local_conf);

  struct gkyl_range skin_conf, ghost_conf;
  gkyl_skin_ghost_ranges(&skin_conf, &ghost_conf, dir, edge, &local_conf_ext, ghost_cells_conf);

  int ghost_cells[pdim];
  for (int d=0; d<cdim; d++) ghost_cells[d] = 1;
  for (int d=cdim; d<pdim; d++) ghost_cells[d] = 0;
  struct gkyl_range local, local_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost_cells, &local_ext, &local);

  struct gkyl_range skin, ghost;
  gkyl_skin_ghost_ranges(&skin, &ghost, dir, edge, &local_ext, ghost_cells);

  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size)
                                       : gkyl_array_acquire(distf);
  struct gkyl_array *jac = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *jac_ho = use_gpu? mkarr(false, jac->ncomp, jac->size)
                                     : gkyl_array_acquire(jac);
  struct gkyl_array *jf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *jf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size)
                                    : gkyl_array_acquire(jf);

  // Project jac onto the basis.
  struct gkyl_eval_on_nodes *proj_jac = gkyl_eval_on_nodes_new(&grid_conf, &basis_conf,
    1, eval_jac_2x, &test_ctx);
  gkyl_eval_on_nodes_advance(proj_jac, 0.0, &local_conf_ext, jac_ho);
  gkyl_eval_on_nodes_release(proj_jac);
  gkyl_array_copy(jac, jac_ho);

  // Project f onto the basis.
  struct gkyl_proj_on_basis *proj_f = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_distf_2x2v, &test_ctx);
  gkyl_proj_on_basis_advance(proj_f, 0.0, &local_ext, distf_ho);
  gkyl_proj_on_basis_release(proj_f);
  gkyl_array_copy(distf, distf_ho);

  // Multiply jac * f.
  gkyl_dg_mul_conf_phase_op_range(&basis_conf, &basis, jf, jac, distf, &local_conf_ext, &local_ext);
  // Place jf in distf as we'll need it to check results.
  gkyl_array_copy(distf_ho, jf);

  // Divide jf by j in the ghost cell, and multiply by the flipped skin cell j.
  struct gkyl_rescale_ghost_jacf* jf_rescale =
    gkyl_rescale_ghost_jacf_new(dir, edge, &basis_conf, &basis, use_gpu);

  gkyl_rescale_ghost_jacf_advance(jf_rescale,
    &skin_conf, &ghost_conf, &ghost, jac, jf);

  gkyl_rescale_ghost_jacf_release(jf_rescale);

  // Check the results.
  gkyl_array_copy(jf_ho, jf);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ghost);
  while (gkyl_range_iter_next(&iter)) {
    int cidx_skin[cdim];
    for (int d=0; d<cdim; d++) cidx_skin[d] = iter.idx[d];
    cidx_skin[dir] = edge == GKYL_LOWER_EDGE? iter.idx[dir]+1 : iter.idx[dir]-1;

    long clinidx_skin = gkyl_range_idx(&skin_conf, cidx_skin);
    long clinidx_ghost = gkyl_range_idx(&ghost_conf, iter.idx);
    long plinidx_ghost = gkyl_range_idx(&ghost, iter.idx);

    const double *jskin_c = gkyl_array_cfetch(jac_ho, clinidx_skin);
    const double *jghost_c = gkyl_array_cfetch(jac_ho, clinidx_ghost);
    const double *distf_c = gkyl_array_cfetch(distf_ho, plinidx_ghost);

    double ref_c[basis.num_basis];
    // Compute accepted results.
    for (int i=0; i<basis.num_basis; ++i)
      ref_c[i] = 0.0;

    if (dir == 0) {
      if (edge == GKYL_LOWER_EDGE)
        accepted_results_kernel_2x2v_lowerx(jskin_c, jghost_c, distf_c, ref_c);
      else
        accepted_results_kernel_2x2v_upperx(jskin_c, jghost_c, distf_c, ref_c);
    }
    else if (dir == 1) {
      if (edge == GKYL_LOWER_EDGE)
        accepted_results_kernel_2x2v_lowery(jskin_c, jghost_c, distf_c, ref_c);
      else
        accepted_results_kernel_2x2v_uppery(jskin_c, jghost_c, distf_c, ref_c);
    }

    const double *jf_c = gkyl_array_cfetch(jf_ho, plinidx_ghost);
    // for (int i=0; i<basis.num_basis; ++i) {
      // TEST_CHECK( gkyl_compare(ref_c[i], jf_c[i], 1e-10) );
      // TEST_MSG("Expected: %.13e | Got:%.13e | Cell:%d,%d,%d,%d\n", ref_c[i], jf_c[i], iter.idx[0], iter.idx[1], iter.idx[2], iter.idx[3]);
    // }
  }

 // Free memory.
  gkyl_array_release(distf_ho);
  gkyl_array_release(distf);
  gkyl_array_release(jac_ho);
  gkyl_array_release(jac);
  gkyl_array_release(jf_ho);
  gkyl_array_release(jf);
}

void
test_1x1v_ho()
{
  test_1x1v_at_edge(false, 0, GKYL_LOWER_EDGE);
  test_1x1v_at_edge(false, 0, GKYL_UPPER_EDGE);
}

void
test_2x2v_ho()
{
  test_2x2v_at_edge(false, 0, GKYL_LOWER_EDGE);
  test_2x2v_at_edge(false, 0, GKYL_UPPER_EDGE);
  test_2x2v_at_edge(false, 1, GKYL_LOWER_EDGE);
  test_2x2v_at_edge(false, 1, GKYL_UPPER_EDGE);
}

#ifdef GKYL_HAVE_CUDA
void
test_1x1v_dev()
{
  test_1x1v_at_edge(true, 0, GKYL_LOWER_EDGE);
  test_1x1v_at_edge(true, 0, GKYL_UPPER_EDGE);
}

void
test_2x2v_dev()
{
  test_2x2v_at_edge(true, 0, GKYL_LOWER_EDGE);
  test_2x2v_at_edge(true, 0, GKYL_UPPER_EDGE);
  test_2x2v_at_edge(true, 1, GKYL_LOWER_EDGE);
  test_2x2v_at_edge(true, 1, GKYL_UPPER_EDGE);
}
#endif

TEST_LIST = {
  { "test_1x1v_ho", test_1x1v_ho },
  { "test_2x2v_ho", test_2x2v_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_1x1v_dev", test_1x1v_dev },
  { "test_2x2v_dev", test_2x2v_dev },
#endif
  { NULL, NULL },
};
