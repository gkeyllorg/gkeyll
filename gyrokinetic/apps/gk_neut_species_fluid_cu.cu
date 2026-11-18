/* -*- c++ -*- */
extern "C" {
#include <float.h>
#include <math.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_cart_modal_serendip_priv.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gk_neut_species_fluid_cu.h>
}

__global__ static void
dirichlet_ghost_cu_ker(
  int cdim, int dir, int skin_idx_dir, double mass, double number_density,
  struct gkyl_range local_ext, struct gkyl_range ghost, struct gkyl_array *density,
  struct gkyl_array *tensor
)
{
  int gidx[GKYL_MAX_CDIM], sidx[GKYL_MAX_CDIM];
  const int nb = 1 << cdim, nK = cdim * cdim;
  const double rho_b_modal = pow(sqrt(2.0), cdim) * mass * number_density;
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < ghost.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&ghost, tid, gidx);
    for (int d = 0; d < cdim; ++d) {
      sidx[d] = gidx[d];
    }
    sidx[dir] = skin_idx_dir;
    const double *rho_skin =
      (const double *)gkyl_array_cfetch(density, gkyl_range_idx(&local_ext, sidx));
    double *rho_ghost = (double *)gkyl_array_fetch(density, gkyl_range_idx(&local_ext, gidx));
    fos_list[cdim].fs[1](dir, rho_skin, rho_ghost);
    for (int k = 0; k < nb; ++k) {
      rho_ghost[k] = -rho_ghost[k];
    }
    rho_ghost[0] += 2.0 * rho_b_modal;

    const double *Kskin =
      (const double *)gkyl_array_cfetch(tensor, gkyl_range_idx(&local_ext, sidx));
    double *Kghost = (double *)gkyl_array_fetch(tensor, gkyl_range_idx(&local_ext, gidx));
    for (int k = 0; k < nK; ++k) {
      fos_list[cdim].fs[1](dir, &Kskin[k * nb], &Kghost[k * nb]);
    }
  }
}

void
gk_neut_fluid_diffusion_dirichlet_ghost_cu(
  int cdim, int dir, int skin_idx_dir, double mass, double number_density,
  const struct gkyl_range *local_ext, const struct gkyl_range *ghost, struct gkyl_array *density,
  struct gkyl_array *tensor
)
{
  dirichlet_ghost_cu_ker<<<ghost->nblocks, ghost->nthreads>>>(
    cdim, dir, skin_idx_dir, mass, number_density, *local_ext, *ghost, density->on_dev,
    tensor->on_dev
  );
}

__global__ static void
zero_flux_ghost_cu_ker(
  int cdim, int dir, int skin_idx_dir, struct gkyl_range local_ext, struct gkyl_range ghost,
  struct gkyl_array *density, struct gkyl_array *tensor
)
{
  int gidx[GKYL_MAX_CDIM], sidx[GKYL_MAX_CDIM];
  const int nb = 1 << cdim, nK = cdim * cdim;
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < ghost.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&ghost, tid, gidx);
    for (int d = 0; d < cdim; ++d) {
      sidx[d] = gidx[d];
    }
    sidx[dir] = skin_idx_dir;
    long sloc = gkyl_range_idx(&local_ext, sidx), gloc = gkyl_range_idx(&local_ext, gidx);
    fos_list[cdim].fs[1](
      dir, (const double *)gkyl_array_cfetch(density, sloc),
      (double *)gkyl_array_fetch(density, gloc)
    );
    const double *Kskin = (const double *)gkyl_array_cfetch(tensor, sloc);
    double *Kghost = (double *)gkyl_array_fetch(tensor, gloc);
    for (int k = 0; k < nK; ++k) {
      fes_list[cdim].fs[1](dir, &Kskin[k * nb], &Kghost[k * nb]);
    }
  }
}

void
gk_neut_fluid_diffusion_zero_flux_ghost_cu(
  int cdim, int dir, int skin_idx_dir, const struct gkyl_range *local_ext,
  const struct gkyl_range *ghost, struct gkyl_array *density, struct gkyl_array *tensor
)
{
  zero_flux_ghost_cu_ker<<<ghost->nblocks, ghost->nthreads>>>(
    cdim, dir, skin_idx_dir, *local_ext, *ghost, density->on_dev, tensor->on_dev
  );
}

__global__ static void
recycling_source_cu_ker(
  int cdim, int dir, int skin_idx_dir, struct gkyl_range local_ext, struct gkyl_range ghost,
  const struct gkyl_array *flux, struct gkyl_array *rhs
)
{
  int gidx[GKYL_MAX_CDIM], sidx[GKYL_MAX_CDIM];
  const int nb = 1 << cdim;
  double src[8];
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < ghost.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&ghost, tid, gidx);
    for (int d = 0; d < cdim; ++d) {
      sidx[d] = gidx[d];
    }
    sidx[dir] = skin_idx_dir;
    fos_list[cdim].fs[1](
      dir, (const double *)gkyl_array_cfetch(flux, gkyl_range_idx(&local_ext, gidx)), src
    );
    // Each ghost cell maps to a distinct skin cell, so no atomics are needed.
    double *rhs_skin = (double *)gkyl_array_fetch(rhs, gkyl_range_idx(&local_ext, sidx));
    for (int k = 0; k < nb; ++k) {
      rhs_skin[k] += src[k];
    }
  }
}

void
gk_neut_fluid_recycling_source_cu(
  int cdim, int dir, int skin_idx_dir, const struct gkyl_range *local_ext,
  const struct gkyl_range *ghost, const struct gkyl_array *flux, struct gkyl_array *rhs
)
{
  recycling_source_cu_ker<<<ghost->nblocks, ghost->nthreads>>>(
    cdim, dir, skin_idx_dir, *local_ext, *ghost, flux->on_dev, rhs->on_dev
  );
}

__global__ static void
diffusion_floor_cu_ker(
  int cdim, double nu_min, struct gkyl_range range, const struct gkyl_array *vt_sq,
  const struct gkyl_array *nu, struct gkyl_array *D
)
{
  int idx[GKYL_MAX_CDIM];
  const int nb = 1 << cdim;
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < range.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&range, tid, idx);
    long loc = gkyl_range_idx(&range, idx);
    const double *vt = (const double *)gkyl_array_cfetch(vt_sq, loc);
    const double *nup = (const double *)gkyl_array_cfetch(nu, loc);
    double *Dp = (double *)gkyl_array_fetch(D, loc);
    double nodes[8 * GKYL_MAX_CDIM] = {0.0}, Dnodal[8] = {0.0};
    nl_list[cdim].nl[1](nodes);
    for (int q = 0; q < nb; ++q) {
      double vtq = eve_list[cdim].ev[1](&nodes[q * cdim], vt);
      double nuq = eve_list[cdim].ev[1](&nodes[q * cdim], nup);
      Dnodal[q] = fmax(vtq, 0.0) / fmax(nuq, nu_min);
    }
    n2m_list[cdim].n2m[1](Dnodal, Dp);
  }
}

void
gk_neut_fluid_diffusion_floor_cu(
  int cdim, double nu_min, const struct gkyl_range *range, const struct gkyl_array *vt_sq,
  const struct gkyl_array *nu, struct gkyl_array *D
)
{
  diffusion_floor_cu_ker<<<range->nblocks, range->nthreads>>>(
    cdim, nu_min, *range, vt_sq->on_dev, nu->on_dev, D->on_dev
  );
}

__global__ static void
vec_set_cu_ker(
  struct gkyl_range range, struct gkyl_array *out, double a, const struct gkyl_array *x, double b,
  const struct gkyl_array *y
)
{
  int idx[GKYL_MAX_CDIM];
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < range.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&range, tid, idx);
    long loc = gkyl_range_idx(&range, idx);
    double *op = (double *)gkyl_array_fetch(out, loc);
    const double *xp = (const double *)gkyl_array_cfetch(x, loc);
    const double *yp = y ? (const double *)gkyl_array_cfetch(y, loc) : 0;
    for (unsigned k = 0; k < out->ncomp; ++k) {
      op[k] = a * xp[k] + (yp ? b * yp[k] : 0.0);
    }
  }
}

void
gk_neut_fluid_diffusion_vec_set_cu(
  const struct gkyl_range *range, struct gkyl_array *out, double a, const struct gkyl_array *x,
  double b, const struct gkyl_array *y
)
{
  vec_set_cu_ker<<<range->nblocks, range->nthreads>>>(
    *range, out->on_dev, a, x->on_dev, b, y ? y->on_dev : 0
  );
}

__global__ static void
p_update_cu_ker(
  struct gkyl_range range, struct gkyl_array *p, const struct gkyl_array *r,
  const struct gkyl_array *v, double beta, double omega
)
{
  int idx[GKYL_MAX_CDIM];
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < range.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&range, tid, idx);
    long loc = gkyl_range_idx(&range, idx);
    double *pp = (double *)gkyl_array_fetch(p, loc);
    const double *rp = (const double *)gkyl_array_cfetch(r, loc);
    const double *vp = (const double *)gkyl_array_cfetch(v, loc);
    for (unsigned k = 0; k < p->ncomp; ++k) {
      pp[k] = rp[k] + beta * (pp[k] - omega * vp[k]);
    }
  }
}

void
gk_neut_fluid_diffusion_p_update_cu(
  const struct gkyl_range *range, struct gkyl_array *p, const struct gkyl_array *r,
  const struct gkyl_array *v, double beta, double omega
)
{
  p_update_cu_ker<<<range->nblocks, range->nthreads>>>(
    *range, p->on_dev, r->on_dev, v->on_dev, beta, omega
  );
}

__global__ static void
dot_cu_ker(
  struct gkyl_range range, const struct gkyl_array *a, const struct gkyl_array *b, double *sum
)
{
  int idx[GKYL_MAX_CDIM];
  double thread_sum = 0.0;
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < range.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&range, tid, idx);
    long loc = gkyl_range_idx(&range, idx);
    const double *ap = (const double *)gkyl_array_cfetch(a, loc);
    const double *bp = (const double *)gkyl_array_cfetch(b, loc);
    for (unsigned k = 0; k < a->ncomp; ++k) {
      thread_sum += ap[k] * bp[k];
    }
  }
  atomicAdd(sum, thread_sum);
}

double
gk_neut_fluid_diffusion_dot_cu(
  const struct gkyl_range *range, const struct gkyl_array *a, const struct gkyl_array *b,
  double *dot_cu
)
{
  gkyl_cu_memset(dot_cu, 0, sizeof(double));
  dot_cu_ker<<<range->nblocks, range->nthreads>>>(*range, a->on_dev, b->on_dev, dot_cu);
  double sum;
  gkyl_cu_memcpy(&sum, dot_cu, sizeof(double), GKYL_CU_MEMCPY_D2H);
  return sum;
}

__global__ static void
positivity_cu_ker(int cdim, struct gkyl_range range, struct gkyl_array *rho)
{
  int idx[GKYL_MAX_CDIM];
  const int nb = 1 << cdim;
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < range.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&range, tid, idx);
    double *rhop = (double *)gkyl_array_fetch(rho, gkyl_range_idx(&range, idx));
    const double rho_avg = rhop[0] / pow(sqrt(2.0), cdim);
    if (rho_avg <= 0.0) {
      for (int k = 0; k < nb; ++k) {
        rhop[k] = 0.0;
      }
      continue;
    }
    double nodes[8 * GKYL_MAX_CDIM] = {0.0};
    nl_list[cdim].nl[1](nodes);
    double rho_min = DBL_MAX;
    for (int q = 0; q < nb; ++q) {
      rho_min = fmin(rho_min, eve_list[cdim].ev[1](&nodes[q * cdim], rhop));
    }
    if (rho_min < 0.0) {
      double theta = fmin(1.0, rho_avg / (rho_avg - rho_min)) * (1.0 - 1.0e-14);
      for (int k = 1; k < nb; ++k) {
        rhop[k] *= theta;
      }
    }
  }
}

void
gk_neut_fluid_diffusion_positivity_cu(
  int cdim, const struct gkyl_range *range, struct gkyl_array *rho
)
{
  positivity_cu_ker<<<range->nblocks, range->nthreads>>>(cdim, *range, rho->on_dev);
}
