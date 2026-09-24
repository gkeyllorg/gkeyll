/* -*- c++ -*- */

#include <cstdio>
#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_util.h>
}

#include <cstdio>
#include <cassert>

// start ID for use in various loops
#define START_ID (threadIdx.x + blockIdx.x*blockDim.x)

// NOTE: This is duplicated in dg_bin_ops_cu. Should be cleaned up 01/05/22
static void
gkyl_get_array_range_kernel_launch_dims(dim3* dimGrid, dim3* dimBlock, gkyl_range range, int ncomp)
{
  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension)
  int ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;
  // CUDA Max block size in x is 2^31 - 1, Max block size in y is 2^16-1
  // Thus, x block size should be bigger to avoid max block size limits
  dimBlock->y = GKYL_MIN2(ncomp*ac1, GKYL_DEFAULT_NUM_THREADS);
  dimGrid->y = gkyl_int_div_up(ncomp*ac1, dimBlock->y);
  dimBlock->x = gkyl_int_div_up(GKYL_DEFAULT_NUM_THREADS, ncomp*ac1);
  dimGrid->x = gkyl_int_div_up(range.volume, ac1*dimBlock->x);
}

__global__ void
gkyl_array_clear_cu_kernel(struct gkyl_array* out, double val)
{
  double *out_d = (double*) out->data;
  for (unsigned long linc = START_ID; linc < NELM(out); linc += blockDim.x*gridDim.x)
    out_d[linc] = val;
}

__global__ void
gkyl_array_accumulate_cu_kernel(struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  double *out_d = (double*) out->data;
  const double *inp_d = (const double*) inp->data;
  for (unsigned long linc = START_ID; linc < NELM(out); linc += blockDim.x*gridDim.x)
    out_d[linc] += a*inp_d[linc];
}

__global__ void
gkyl_array_accumulate_offset_cu_kernel(struct gkyl_array* out, double a,
  const struct gkyl_array* inp, int coff)
{
  double *out_d = (double*) out->data;
  const double *inp_d = (const double*) inp->data;
  if (NCOM(out) < NCOM(inp)) {
    for (unsigned long linc = START_ID; linc < NSIZE(out); linc += blockDim.x*gridDim.x)
      for (unsigned k=0; k<NCOM(out); ++k)
        out_d[linc*NCOM(out)+k] += a*inp_d[linc*NCOM(inp)+coff+k];
  } else {
    for (unsigned long linc = START_ID; linc < NSIZE(out); linc += blockDim.x*gridDim.x)
      for (unsigned k=0; k<NCOM(inp); ++k)
        out_d[linc*NCOM(out)+coff+k] += a*inp_d[linc*NCOM(inp)+k];
  }
}

__global__ void
gkyl_array_set_cu_kernel(struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  double *out_d = (double*) out->data;
  const double *inp_d = (const double*) inp->data;
  for (unsigned long linc = START_ID; linc < NELM(out); linc += blockDim.x*gridDim.x)
    out_d[linc] = a*inp_d[linc];
}

__global__ void
gkyl_array_set_offset_cu_kernel(struct gkyl_array* out, double a,
  const struct gkyl_array* inp, int coff)
{
  double *out_d = (double*) out->data;
  const double *inp_d = (const double*) inp->data;
  if (NCOM(out) < NCOM(inp)) {
    for (unsigned long linc = START_ID; linc < NSIZE(out); linc += blockDim.x*gridDim.x)
      for (unsigned k=0; k<NCOM(out); ++k)
        out_d[linc*NCOM(out)+k] = a*inp_d[linc*NCOM(inp)+coff+k];
  } else {
    for (unsigned long linc = START_ID; linc < NSIZE(out); linc += blockDim.x*gridDim.x)
      for (unsigned k=0; k<NCOM(inp); ++k)
        out_d[linc*NCOM(out)+coff+k] = a*inp_d[linc*NCOM(inp)+k];
  }
}

__global__ void
gkyl_array_scale_by_cell_cu_kernel(struct gkyl_array* out, const struct gkyl_array* a)
{
  double *out_d = (double*) out->data;
  const double *a_d = (double*) a->data;
  for (unsigned long linc = START_ID; linc < NELM(out); linc += blockDim.x*gridDim.x)
    out_d[linc] = a_d[linc/out->ncomp]*out_d[linc];
}

__global__ void
gkyl_array_divide_by_cell_cu_kernel(struct gkyl_array* out, const struct gkyl_array* a)
{
  double *out_d = (double*) out->data;
  const double *a_d = (double*) a->data;
  for (unsigned long linc = START_ID; linc < NELM(out); linc += blockDim.x*gridDim.x)
    out_d[linc] = out_d[linc]/a_d[linc/out->ncomp];
}

__global__ void
gkyl_array_invert_by_cell_cu_kernel(struct gkyl_array* out, const struct gkyl_array *inp)
{
  double *out_d = (double*) out->data;
  const double *inp_d = (const double*) inp->data;
  for (unsigned long linc = START_ID; linc < NELM(out); linc += blockDim.x*gridDim.x)
    out_d[linc] = 1.0/inp_d[linc];
}

__global__ void
gkyl_array_shiftc_cu_kernel(struct gkyl_array* out, double a, unsigned k)
{
  double *out_d = (double*) out->data;
  for (unsigned long linc = START_ID; linc < NSIZE(out); linc += blockDim.x*gridDim.x)
    out_d[linc*out->ncomp+k] = a+out_d[linc*out->ncomp+k];
}

__global__ void
gkyl_array_min_by_cell_cu_kernel(struct gkyl_array* out, const struct gkyl_array *inp, double a)
{
  double *out_d = (double*) out->data;
  const double *inp_d = (const double*) inp->data;
  for (unsigned long linc = START_ID; linc < NELM(out); linc += blockDim.x*gridDim.x)
    out_d[linc] = fmin(inp_d[linc], a);
}

// Host-side wrappers for array operations
void
gkyl_array_clear_cu(struct gkyl_array* out, double val)
{
  gkyl_array_clear_cu_kernel<<<out->nblocks, out->nthreads>>>(out->on_dev, val);
}

void
gkyl_array_accumulate_cu(struct gkyl_array* out, double a, const struct gkyl_array* inp)
{
  gkyl_array_accumulate_cu_kernel<<<out->nblocks, out->nthreads>>>(out->on_dev, a, inp->on_dev);
}

void
gkyl_array_accumulate_offset_cu(struct gkyl_array* out, double a, const struct gkyl_array* inp, int coff)
{
  int nblocks = gkyl_int_div_up(out->size, out->nthreads);
  gkyl_array_accumulate_offset_cu_kernel<<<nblocks, out->nthreads>>>(out->on_dev, a, inp->on_dev, coff);
}

void
gkyl_array_set_cu(struct gkyl_array* out, double a, const struct gkyl_array* inp)
{
  gkyl_array_set_cu_kernel<<<out->nblocks, out->nthreads>>>(out->on_dev, a, inp->on_dev);
}

void
gkyl_array_set_offset_cu(struct gkyl_array* out, double a, const struct gkyl_array* inp, int coff)
{
  int nblocks = gkyl_int_div_up(out->size, out->nthreads);
  gkyl_array_set_offset_cu_kernel<<<nblocks, out->nthreads>>>(out->on_dev, a, inp->on_dev, coff);
}

void
gkyl_array_scale_cu(struct gkyl_array* out, double a)
{
  gkyl_array_set_cu_kernel<<<out->nblocks, out->nthreads>>>(out->on_dev, a, out->on_dev);
}

void
gkyl_array_scale_by_cell_cu(struct gkyl_array* out, const struct gkyl_array* a)
{
  gkyl_array_scale_by_cell_cu_kernel<<<out->nblocks, out->nthreads>>>(out->on_dev, a->on_dev);
}

void
gkyl_array_divide_by_cell_cu(struct gkyl_array* out, const struct gkyl_array* a)
{
  gkyl_array_divide_by_cell_cu_kernel<<<out->nblocks, out->nthreads>>>(out->on_dev, a->on_dev);
}

void
gkyl_array_invert_by_cell_cu(struct gkyl_array* out, const struct gkyl_array *inp)
{
  gkyl_array_invert_by_cell_cu_kernel<<<out->nblocks, out->nthreads>>>(out->on_dev, inp->on_dev);
}

void
gkyl_array_shiftc_cu(struct gkyl_array* out, double a, unsigned k)
{
  gkyl_array_shiftc_cu_kernel<<<out->nblocks, out->nthreads>>>(out->on_dev, a, k);
}

void
gkyl_array_min_by_cell_cu(struct gkyl_array* out, const struct gkyl_array *inp, double a)
{
  gkyl_array_min_by_cell_cu_kernel<<<out->nblocks, out->nthreads>>>(out->on_dev, inp->on_dev, a);
}

// Range-based methods
// Range-based methods need to inverse index from linc to idx.
// Must use gkyl_sub_range_inv_idx so that linc=0 maps to idxc={1,1,...}
// since range can be a subrange.
// Then, convert back to a linear index on the super-range.
// This super range can include ghost cells and thus linear index will have
// have jumps over ghost cells.

__global__ void
gkyl_array_clear_range_cu_kernel(struct gkyl_array *out, double val, struct gkyl_range range)
{
  long n = NCOM(out);
  int idx[GKYL_MAX_DIM];
  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension)
  long ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;

  // 2D thread grid
  // linc2 = c + n*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  // linc1 = idx2 + ac2*idx3 + ...
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume/ac1;
      linc1 += gridDim.x*blockDim.x)
  {
    // full linear cell index (not including components) is 
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc1.
    // we want to find the start linear index of each contiguous data block, 
    // which corresponds to idx1 = 0. 
    // so linear index of start of contiguous block is ac1*linc2.
    gkyl_sub_range_inv_idx(&range, ac1*linc1, idx);
    long start = gkyl_range_idx(&range, idx);
    
    double* out_d = (double*) gkyl_array_fetch(out, start);
    // do operation on contiguous data block
    if (linc2 < n*ac1)
      out_d[linc2] = val;
  }
}

__global__ void
gkyl_array_accumulate_range_cu_kernel(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, struct gkyl_range range)
{
  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n = outnc<inpnc ? outnc : inpnc;
  int idx[GKYL_MAX_DIM];

  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension)
  long ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;

  // 2D thread grid
  // linc2 = c + n*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  long c = linc2 % n;
  long idx1 = linc2 / n;
  // get corresponding linc1 index for inp and out 
  // (one of these will not be contiguous if outnc!=inpnc)
  long linc2_in = c + inpnc*idx1; 
  long linc2_out = c + outnc*idx1; 
  // linc1 = idx2 + ac2*idx3 + ...
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume/ac1;
      linc1 += gridDim.x*blockDim.x)
  {
    // full linear cell index (not including components) is 
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc1.
    // we want to find the start linear index of each contiguous data block, 
    // which corresponds to idx1 = 0. 
    // so linear index of start of contiguous block is ac1*linc2.
    gkyl_sub_range_inv_idx(&range, ac1*linc1, idx);
    long start = gkyl_range_idx(&range, idx);
    
    double* out_d = (double*) gkyl_array_fetch(out, start);
    const double* inp_d = (const double*) gkyl_array_cfetch(inp, start);
    // do operation on contiguous data block
    if (linc2 < n*ac1)
      out_d[linc2_out] += a*inp_d[linc2_in];
  }
}

__global__ void
gkyl_array_accumulate_offset_range_cu_kernel(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, int coff, struct gkyl_range range)
{
  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n;
  int outoff, inpoff;
  if (outnc < inpnc) {
     n = outnc;
     outoff = 0;
     inpoff = coff;
  } else {
     n = inpnc;
     outoff = coff;
     inpoff = 0;
  }
  int idx[GKYL_MAX_DIM];

  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension)
  long ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;

  // 2D thread grid
  // linc2 = c + n*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  long c = linc2 % n;
  long idx1 = linc2 / n;
  // get corresponding linc1 index for inp and out 
  // (one of these will not be contiguous if outnc!=inpnc)
  long linc2_in = c + inpnc*idx1; 
  long linc2_out = c + outnc*idx1; 
  // linc1 = idx2 + ac2*idx3 + ...
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume/ac1;
      linc1 += gridDim.x*blockDim.x)
  {
    // full linear cell index (not including components) is 
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc1.
    // we want to find the start linear index of each contiguous data block, 
    // which corresponds to idx1 = 0. 
    // so linear index of start of contiguous block is ac1*linc2.
    gkyl_sub_range_inv_idx(&range, ac1*linc1, idx);
    long start = gkyl_range_idx(&range, idx);
    
    double* out_d = (double*) gkyl_array_fetch(out, start);
    const double* inp_d = (const double*) gkyl_array_cfetch(inp, start);
    // do operation on contiguous data block
    if (linc2 < n*ac1)
      out_d[linc2_out+outoff] += a*inp_d[linc2_in+inpoff];
  }
}

__global__ void
gkyl_array_set_range_cu_kernel(struct gkyl_array *out, double a,
  const struct gkyl_array* inp, struct gkyl_range out_range, struct gkyl_range inp_range)
{
  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n = outnc<inpnc ? outnc : inpnc;
  int idx_out[GKYL_MAX_DIM], idx_inp[GKYL_MAX_DIM];

  int ndim = inp_range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension)
  long ac1 = inp_range.iac[ndim-1] > 0 ? inp_range.iac[ndim-1] : 1;

  // 2D thread grid
  // linc2 = c + n*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  long c = linc2 % n;
  long idx1 = linc2 / n;
  // get corresponding linc1 index for inp and out 
  // (one of these will not be contiguous if outnc!=inpnc)
  long linc2_in = c + inpnc*idx1; 
  long linc2_out = c + outnc*idx1; 
  // linc1 = idx2 + ac2*idx3 + ...
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < inp_range.volume/ac1;
      linc1 += gridDim.x*blockDim.x)
  {
    // full linear cell index (not including components) is 
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc1.
    // we want to find the start linear index of each contiguous data block, 
    // which corresponds to idx1 = 0. 
    // so linear index of start of contiguous block is ac1*linc2.
    gkyl_sub_range_inv_idx(&out_range, ac1*linc1, idx_out);
    gkyl_sub_range_inv_idx(&inp_range, ac1*linc1, idx_inp);
    long start_out = gkyl_range_idx(&out_range, idx_out);
    long start_inp = gkyl_range_idx(&inp_range, idx_inp);
    
    double* out_d = (double*) gkyl_array_fetch(out, start_out);
    const double* inp_d = (const double*) gkyl_array_cfetch(inp, start_inp);
    // do operation on contiguous data block
    if (linc2 < n*ac1)
      out_d[linc2_out] = a*inp_d[linc2_in];
  }
}

__global__ void
gkyl_array_set_offset_range_cu_kernel(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, int coff, struct gkyl_range range)
{
  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n;
  int outoff, inpoff;
  if (outnc < inpnc) {
     n = outnc;
     outoff = 0;
     inpoff = coff;
  } else {
     n = inpnc;
     outoff = coff;
     inpoff = 0;
  }
  int idx[GKYL_MAX_DIM];

  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension)
  long ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;

  // 2D thread grid
  // linc2 = c + n*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  long c = linc2 % n;
  long idx1 = linc2 / n;
  // get corresponding linc1 index for inp and out 
  // (one of these will not be contiguous if outnc!=inpnc)
  long linc2_in = c + inpnc*idx1; 
  long linc2_out = c + outnc*idx1; 
  // linc1 = idx2 + ac2*idx3 + ...
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume/ac1;
      linc1 += gridDim.x*blockDim.x)
  {
    // full linear cell index (not including components) is 
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc1.
    // we want to find the start linear index of each contiguous data block, 
    // which corresponds to idx1 = 0. 
    // so linear index of start of contiguous block is ac1*linc2.
    gkyl_sub_range_inv_idx(&range, ac1*linc1, idx);
    long start = gkyl_range_idx(&range, idx);
    
    double* out_d = (double*) gkyl_array_fetch(out, start);
    const double* inp_d = (const double*) gkyl_array_cfetch(inp, start);
    // do operation on contiguous data block
    if (linc2 < n*ac1)
      out_d[linc2_out+outoff] = a*inp_d[linc2_in+inpoff];
  }
}

__global__ void
gkyl_array_shiftc_range_cu_kernel(struct gkyl_array* out, double a,
  unsigned k, struct gkyl_range range)
{
  long ncomp = NCOM(out);
  int idx[GKYL_MAX_DIM];
  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension).
  long ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;

  // 2D thread grid
  // linc2 = c + ncomp*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  // linc1 = idx2 + ac2*idx3 + ...
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume/ac1;
      linc1 += gridDim.x*blockDim.x)
  {
    // full linear cell index (not including components) is
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc1.
    // we want to find the start linear index of each contiguous data block,
    // which corresponds to idx1 = 0.
    // so linear index of start of contiguous block is ac1*linc2.
    gkyl_sub_range_inv_idx(&range, ac1*linc1, idx);
    long start = gkyl_range_idx(&range, idx);

    double* out_d = (double*) gkyl_array_fetch(out, start);
    // do operation on contiguous data block
    if (linc2*ncomp < ncomp*ac1)
      out_d[linc2*ncomp+k] += a;
  }
}

__global__ void
gkyl_array_min_by_cell_range_cu_kernel(struct gkyl_array* out, const struct gkyl_array* inp,
  double a, struct gkyl_range range)
{
  long ncomp = NCOM(out);
  int idx[GKYL_MAX_DIM];
  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension).
  long ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;

  // 2D thread grid
  // linc2 = c + ncomp*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  // linc1 = idx2 + ac2*idx3 + ...
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume/ac1;
      linc1 += gridDim.x*blockDim.x)
  {
    // full linear cell index (not including components) is
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc1.
    // we want to find the start linear index of each contiguous data block,
    // which corresponds to idx1 = 0.
    // so linear index of start of contiguous block is ac1*linc2.
    gkyl_sub_range_inv_idx(&range, ac1*linc1, idx);
    long start = gkyl_range_idx(&range, idx);

    double* out_d = (double*) gkyl_array_fetch(out, start);
    const double* inp_d = (const double*) gkyl_array_cfetch(inp, start);
    // do operation on contiguous data block
    if (linc2 < ncomp*ac1)
      out_d[linc2] = fmin(inp_d[linc2], a);
  }
}

__global__ void 
gkyl_array_copy_range_cu_kernel(struct gkyl_array *out, const struct gkyl_array* inp,
  struct gkyl_range out_range, struct gkyl_range inp_range)
{
  int idx_out[GKYL_MAX_DIM], idx_inp[GKYL_MAX_DIM];
  long n = NCOM(out); // assume ncomp_in == ncomp_out
  int ndim = inp_range.ndim;
  // ac1 = size of last dimension in input range (fastest moving dimension)
  long ac1 = inp_range.iac[ndim-1] > 0 ? inp_range.iac[ndim-1] : 1;

  // 2D thread grid
  // linc2 = c + n*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  // linc1 = idx2 + ac2*idx3 + ...
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < inp_range.volume/ac1;
      linc1 += gridDim.x*blockDim.x)
  {
    // full linear cell index (not including components) is 
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc1.
    // we want to find the start linear index of each contiguous data block, 
    // which corresponds to idx1 = 0. 
    // so linear index of start of contiguous block is ac1*linc1.
    // NOTE: the above necessarily applies only to inp_range
    gkyl_sub_range_inv_idx(&out_range, ac1*linc1, idx_out);
    gkyl_sub_range_inv_idx(&inp_range, ac1*linc1, idx_inp);
    long start_out = gkyl_range_idx(&out_range, idx_out);
    long start_inp = gkyl_range_idx(&inp_range, idx_inp);
    
    double* out_d = (double*) gkyl_array_fetch(out, start_out);
    const double* inp_d = (const double*) gkyl_array_cfetch(inp, start_inp);
    // do operation on contiguous data block
    if (linc2 < n*ac1)
      out_d[linc2] = inp_d[linc2];
  }
}

__global__ void 
gkyl_array_copy_to_buffer_cu_kernel(void *data, const struct gkyl_array *arr,
  struct gkyl_range range)
{
  double *d_data = (double*) data;
  int idx[GKYL_MAX_DIM];
  long n = NCOM(arr); // assume ncomp_in == ncomp_out
  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension)
  long ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;

  // 2D thread grid
  // linc2 = c + n*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  // linc1 = idx2 + ac2*idx3 + ...
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume/ac1;
      linc1 += gridDim.x*blockDim.x)
  {
    // full linear cell index (not including components) is 
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc1.
    // we want to find the start linear index of each contiguous data block, 
    // which corresponds to idx1 = 0. 
    // so linear index of start of contiguous block is ac1*linc2.
    gkyl_sub_range_inv_idx(&range, ac1*linc1, idx);
    long start = gkyl_range_idx(&range, idx);
    
    const double* arr_d = (const double*) gkyl_array_cfetch(arr, start);
    // read from contiguous data block
    if (linc2 < n*ac1)
      d_data[linc2 + n*ac1*linc1] = arr_d[linc2];
  }
}

__global__ void 
gkyl_array_copy_from_buffer_cu_kernel(struct gkyl_array *arr, const void *data,
  struct gkyl_range range)
{
  int idx[GKYL_MAX_DIM];
  const double *d_data = (const double*) data;
  long n = NCOM(arr);
  
  // since input data is just a linear array, just stream through data linearly
  // linc = c + n*idx1 + n*ac1*idx2 + ...
  for (unsigned long linc = START_ID; linc < range.volume*n; linc += blockDim.x*gridDim.x) {
    int c = linc % n;
    long linc2 = linc / n; // = idx1 + ac1*idx2 + ...
    gkyl_sub_range_inv_idx(&range, linc2, idx);
    long start = gkyl_range_idx(&range, idx);
    
    double *arr_data = (double*) gkyl_array_fetch(arr, start);
    arr_data[c] = d_data[linc];
  }
}

__global__ void 
gkyl_array_copy_to_buffer_fn_cu_kernel(void *data, const struct gkyl_array *arr,
  struct gkyl_range range, struct gkyl_array_copy_func *cf)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume; linc1 += blockDim.x*gridDim.x) {
    // inverse index from linc1 to idxc
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idxc={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&range, idx);

    const double *inp = (const double*) gkyl_array_cfetch(arr, linc);
    double *out = (double*) gkyl_flat_fetch(data, arr->esznc*linc1);
    cf->func(arr->ncomp, out, inp, cf->ctx);
  }
}

__global__ void 
gkyl_array_flip_copy_to_buffer_fn_cu_kernel(void *data, const struct gkyl_array *arr,
  int dir, struct gkyl_range range, struct gkyl_range buff_range,
  struct gkyl_array_copy_func *cf)
{
  int idx[GKYL_MAX_DIM];
  int fidx[GKYL_MAX_DIM]; // flipped index

  int uplo = range.upper[dir]+range.lower[dir];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume; linc1 += blockDim.x*gridDim.x) {
    // inverse index from linc1 to idxc
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idxc={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    gkyl_copy_int_arr(range.ndim, idx, fidx);
    fidx[dir] = uplo - idx[dir];

    // convert back to a linear index on the super-range (with ghost cells)
    // linc and flipped linc (flinc) will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&range, idx);
    long flinc = gkyl_range_idx(&buff_range, fidx);
    
    const double *inp = (const double*) gkyl_array_cfetch(arr, linc);
    double *out = (double*) gkyl_flat_fetch(data, arr->esznc*flinc);
    cf->func(arr->ncomp, out, inp, cf->ctx);
  }
}

__global__ void
gkyl_array_max_by_cell_per_cell_avg_range_cu_kernel(struct gkyl_array *out,
  const struct gkyl_array* inp, struct gkyl_range range)
{
  long outnc = NCOM(out), inpnc = NCOM(inp);
  // For ceil, we assume component counts match, or we limit to outnc
  long n = outnc; 
  int idx[GKYL_MAX_DIM];

  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension)
  long ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;

  // 2D thread grid
  // In this kernel, linc2 corresponds directly to idx1 (the cell index in the contiguous strip)
  // We do not split by component 'c' here because the check on c=0 must gate the copy of c=1..N
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  
  // Note: linc2 acts as idx1 here. 
  // We effectively handle all 'n' components for the cell at 'idx1' in this single thread.

  // linc1 = idx2 + ac2*idx3 + ...
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume/ac1;
      linc1 += gridDim.x*blockDim.x)
  {
    // full linear cell index calculation (same as accumulate)
    // idx1 + ac1*idx2 + ...
    // We find the start linear index of the contiguous block (idx1=0)
    gkyl_sub_range_inv_idx(&range, ac1*linc1, idx);
    long start = gkyl_range_idx(&range, idx);
    
    double* out_d = (double*) gkyl_array_fetch(out, start);
    const double* inp_d = (const double*) gkyl_array_cfetch(inp, start);

    // do operation on contiguous data block
    if (linc2 < ac1) {
      // Calculate specific cell pointers based on idx1 (linc2)
      // Note: out_d points to the start of the 'ac1' strip.
      // We offset by linc2 * n to get to the specific cell.
      long cell_offset_out = linc2 * outnc;
      long cell_offset_inp = linc2 * inpnc;

      // Check the condition on the 0-th component (Cell Average)
      if (out_d[cell_offset_out] < inp_d[cell_offset_inp]) {
        
        // If condition passed, copy ALL components for this cell
        for (int k = 0; k < n; ++k) {
          out_d[cell_offset_out + k] = inp_d[cell_offset_inp + k];
        }
      }
    }
  }
}

// Host-side wrappers for range-based array operations
void
gkyl_array_clear_range_cu(struct gkyl_array *out, double val, const struct gkyl_range *range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, *range, out->ncomp);

  gkyl_array_clear_range_cu_kernel<<<dimGrid, dimBlock>>>(out->on_dev, val, *range);
}

void
gkyl_array_accumulate_range_cu(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, const struct gkyl_range *range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, *range, min(out->ncomp, inp->ncomp));

  gkyl_array_accumulate_range_cu_kernel<<<dimGrid, dimBlock>>>(out->on_dev, a, inp->on_dev, *range);
}

void
gkyl_array_accumulate_offset_range_cu(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, int coff, const struct gkyl_range *range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, *range, min(out->ncomp, inp->ncomp));

  gkyl_array_accumulate_offset_range_cu_kernel<<<dimGrid, dimBlock>>>(out->on_dev, a, inp->on_dev, coff, *range);
}

void
gkyl_array_set_range_cu(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, const struct gkyl_range *range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, *range, min(out->ncomp, inp->ncomp));

  gkyl_array_set_range_cu_kernel<<<dimGrid, dimBlock>>>(out->on_dev,
    a, inp->on_dev, *range, *range);
}

void
gkyl_array_set_range_to_range_cu(struct gkyl_array *out, double a,
  const struct gkyl_array* inp, const struct gkyl_range *out_range, const struct gkyl_range *inp_range)
{
  if (inp_range->volume > 0) {
    dim3 dimGrid, dimBlock;
    gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, *inp_range, min(out->ncomp, inp->ncomp));

    gkyl_array_set_range_cu_kernel<<<dimGrid, dimBlock>>>(out->on_dev,
      a, inp->on_dev, *out_range, *inp_range);
  }
}

void
gkyl_array_set_offset_range_cu(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, int coff, const struct gkyl_range *range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, *range, min(out->ncomp, inp->ncomp));

  gkyl_array_set_offset_range_cu_kernel<<<dimGrid, dimBlock>>>(out->on_dev,
    a, inp->on_dev, coff, *range);
}

void
gkyl_array_scale_range_cu(struct gkyl_array *out,
  double a, const struct gkyl_range *range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, *range, out->ncomp);

  gkyl_array_set_range_cu_kernel<<<dimGrid, dimBlock>>>(out->on_dev,
    a, out->on_dev, *range, *range);
}

void
gkyl_array_shiftc_range_cu(struct gkyl_array* out, double a, unsigned k, const struct gkyl_range *range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, *range, 1);

  gkyl_array_shiftc_range_cu_kernel<<<dimGrid, dimBlock>>>(out->on_dev, a, k, *range);
}

void
gkyl_array_min_by_cell_range_cu(struct gkyl_array* out, const struct gkyl_array* inp,
  double a, const struct gkyl_range *range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, *range, out->ncomp);

  gkyl_array_min_by_cell_range_cu_kernel<<<dimGrid, dimBlock>>>(out->on_dev, inp->on_dev, a, *range);
}

void
gkyl_array_copy_range_cu(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, *range, out->ncomp);

  gkyl_array_copy_range_cu_kernel<<<dimGrid, dimBlock>>>(out->on_dev,
    inp->on_dev, *range, *range);
}

void
gkyl_array_copy_range_to_range_cu(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *out_range, const struct gkyl_range *inp_range)
{
  if (inp_range->volume > 0) {
    dim3 dimGrid, dimBlock;
    gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, *inp_range, inp->ncomp);

    gkyl_array_copy_range_cu_kernel<<<dimGrid, dimBlock>>>(out->on_dev,
      inp->on_dev, *out_range, *inp_range);
  }
}

void 
gkyl_array_copy_to_buffer_cu(void *data, 
  const struct gkyl_array *arr, const struct gkyl_range *range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, *range, arr->ncomp);

  gkyl_array_copy_to_buffer_cu_kernel<<<dimGrid, dimBlock>>>(data,
    arr->on_dev, *range);
}

void 
gkyl_array_copy_from_buffer_cu(struct gkyl_array *arr,
  const void *data, const struct gkyl_range *range)
{
  int nelem = range->volume*arr->ncomp;
  int nthreads = GKYL_DEFAULT_NUM_THREADS;
  int nblocks = gkyl_int_div_up(nelem, nthreads);
  gkyl_array_copy_from_buffer_cu_kernel<<<nblocks, nthreads>>>(arr->on_dev,
    data, *range);
}

void
gkyl_array_copy_to_buffer_fn_cu(void *data, const struct gkyl_array *arr,
  const struct gkyl_range *range, struct gkyl_array_copy_func *cf)
{
  if (range->volume > 0) {
    int nblocks = range->nblocks;
    int nthreads = range->nthreads;

    gkyl_array_copy_to_buffer_fn_cu_kernel<<<nblocks, nthreads>>>(
      data, arr->on_dev, *range, cf);
  }
}

void
gkyl_array_flip_copy_to_buffer_fn_cu(void *data, const struct gkyl_array *arr,
  int dir, const struct gkyl_range *range, struct gkyl_array_copy_func *cf)
{
  if (range->volume > 0) {
    int nblocks = range->nblocks;
    int nthreads = range->nthreads;

    struct gkyl_range buff_range;
    gkyl_range_init(&buff_range, range->ndim, range->lower, range->upper);
  
    gkyl_array_flip_copy_to_buffer_fn_cu_kernel<<<nblocks, nthreads>>>(data,
      arr->on_dev, dir, *range, buff_range, cf);
  }
}

void
gkyl_array_max_by_cell_per_cell_avg_range_cu(struct gkyl_array *out,
  const struct gkyl_array* inp, const struct gkyl_range *range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, *range, 1);
  gkyl_array_max_by_cell_per_cell_avg_range_cu_kernel<<<dimGrid, dimBlock>>>(out->on_dev, inp->on_dev, *range);
}
