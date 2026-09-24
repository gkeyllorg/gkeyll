/* -*- c++ -*- */
extern "C" {
#include <gkyl_bc_emission_spectrum.h>
#include <gkyl_bc_emission_spectrum_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
}

// start ID for use in various loops
#define START_ID (threadIdx.x + blockIdx.x*blockDim.x)

__global__ static void
gkyl_bc_emission_spectrum_set_exterm_params_cu_ker(struct gkyl_emission_spectrum_model *spectrum_model,
  struct gkyl_emission_yield_model *yield_model, int cdim, int vdim, double mass_in, double mass_out)
{
  spectrum_model->cdim = cdim;
  spectrum_model->vdim = vdim;
  spectrum_model->mass = mass_out;

  yield_model->cdim = cdim;
  yield_model->vdim = vdim;
  yield_model->mass = mass_in;
}

__global__ static void
gkyl_bc_emission_spectrum_sey_calc_cu_ker(struct gkyl_rect_grid grid,
  const struct gkyl_range ghost_r, struct gkyl_array *yield,
  struct gkyl_emission_yield_model *yield_model)
{

  double xc[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM];

  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < ghost_r.volume; linc += blockDim.x*gridDim.x) {

    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&ghost_r, linc, pidx);
    
    long loc = gkyl_range_idx(&ghost_r, pidx);
    double *out = (double *) gkyl_array_fetch(yield, loc);
    gkyl_rect_grid_cell_center(&grid, pidx, xc);
    yield_model->function(out, yield_model, xc);
  }
}

__global__ static void
gkyl_bc_emission_spectrum_advance_cu_weight_ker(int cdim, int dir, enum gkyl_edge_loc edge,
  const struct gkyl_array *bflux, struct gkyl_array *weight, struct gkyl_rect_grid grid,
  struct gkyl_array *yield, const struct gkyl_range impact_buff_r,
  const struct gkyl_range impact_cbuff_r)
{
  double xc[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < impact_buff_r.volume; tid += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&impact_buff_r, tid, pidx);
    
    gkyl_rect_grid_cell_center(&grid, pidx, xc);

    long lincP = gkyl_range_idx(&impact_buff_r, pidx);
    
    const double* inp = (const double*) gkyl_array_cfetch(bflux, lincP);
    const double* gain = (const double*) gkyl_array_cfetch(yield, lincP);
    double wLocal[2];
    for (unsigned int k=0; k<weight->ncomp; ++k)
      wLocal[k] = 0.0;

    bc_weighted_delta(inp, cdim, dir, edge, xc, gain, &wLocal[0]);
    
    // get conf-space linear index.
    for (unsigned int i = 0; i < impact_cbuff_r.ndim; i++)
      cidx[i] = pidx[i];
    long lincC = gkyl_range_idx(&impact_cbuff_r, cidx);

    double* wptr = (double*) gkyl_array_fetch(weight, lincC);
    for (unsigned int k = 0; k < weight->ncomp; ++k) {
       atomicAdd(&wptr[k], wLocal[k]);
    }
  }
}

__global__ static void
gkyl_bc_emission_spectrum_advance_cu_accumulate_ker(const struct gkyl_array *spectrum,
  struct gkyl_array *f_emit, struct gkyl_array *weight, struct gkyl_array *k,
  const struct gkyl_array *flux, const struct gkyl_range emit_buff_r,
  const struct gkyl_range impact_cbuff_r, struct gkyl_emission_spectrum_model *spectrum_model)
{
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < emit_buff_r.volume; tid += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&emit_buff_r, tid, pidx);
    
    // get conf-space linear index.
    for (unsigned int i = 0; i < impact_cbuff_r.ndim; i++)
      cidx[i] = pidx[i];
    long lincC = gkyl_range_idx(&impact_cbuff_r, cidx);

    long lincP = gkyl_range_idx(&emit_buff_r, pidx);
    
    const double* inp = (const double*) gkyl_array_cfetch(spectrum, lincP);
    const double* w = (const double*) gkyl_array_cfetch(weight, lincC);
    const double* boundary_flux = (const double*) gkyl_array_cfetch(flux, lincC);
    double* fac = (double*) gkyl_array_fetch(k, lincC);
    double* out = (double*) gkyl_array_fetch(f_emit, lincP);
    double effective_delta = w[0]/w[1];
    spectrum_model->normalization(fac, spectrum_model, boundary_flux, effective_delta);

    for (int c=0; c<spectrum->ncomp; ++c)
      out[c] += fac[0]*inp[c];

  }
}

void
gkyl_bc_emission_spectrum_set_extern_params_cu(const struct gkyl_bc_emission_spectrum *up,
  int cdim, int vdim, double mass_in, double mass_out)
{
  gkyl_bc_emission_spectrum_set_exterm_params_cu_ker<<<1, 1>>>(up->spectrum_model->on_dev,
    up->yield_model->on_dev, cdim, vdim, mass_in, mass_out);
}

void
gkyl_bc_emission_spectrum_sey_calc_cu(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_array *yield, struct gkyl_rect_grid *grid, const struct gkyl_range *impact_buff_r)
{
  int nblocks = impact_buff_r->nblocks, nthreads = impact_buff_r->nthreads;

  gkyl_bc_emission_spectrum_sey_calc_cu_ker<<<nblocks, nthreads>>>(*grid,
    *impact_buff_r, yield->on_dev, up->yield_model->on_dev);
}

void
gkyl_bc_emission_spectrum_advance_cu(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_range *impact_buff_r, struct gkyl_range *impact_cbuff_r,
  struct gkyl_range *emit_buff_r, struct gkyl_array *bflux, struct gkyl_array *f_emit,
  struct gkyl_array *yield, struct gkyl_array *spectrum, struct gkyl_array *weight,
  struct gkyl_array *flux, struct gkyl_array *k)
{
  int nblocks = impact_buff_r->nblocks, nthreads = impact_buff_r->nthreads;

  // Calculate weighted mean numerator and denominator
  gkyl_bc_emission_spectrum_advance_cu_weight_ker<<<nblocks, nthreads>>>(up->cdim, up->dir,
    up->edge, bflux->on_dev, weight->on_dev, *up->grid, yield->on_dev, *impact_buff_r,
    *impact_cbuff_r);

  nblocks = emit_buff_r->nblocks;
  nthreads = emit_buff_r->nthreads;

  // Finish weighted mean calculation and accumulate to buffer
  gkyl_bc_emission_spectrum_advance_cu_accumulate_ker<<<nblocks, nthreads>>>(spectrum->on_dev,
    f_emit->on_dev, weight->on_dev, k->on_dev, flux->on_dev, *emit_buff_r, *impact_cbuff_r,
    up->spectrum_model->on_dev);
}

