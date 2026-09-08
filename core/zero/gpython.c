/* gpython.c — implementation of the postgkyl support API (gkyl_gpython.h).
 *
 * Compiled into libg0core.so alongside the code it wraps, so the shim and
 * the library can never drift apart. Every struct field access, signature,
 * and by-value convention below is checked by the C compiler against the
 * headers in this tree: a core API change that affects the shim fails this
 * build, here at the producer, instead of corrupting postgkyl's data at
 * runtime. Consumers see only gkyl_gpython.h (opaque handles + scalars).
 */
#include <gkyl_gpython.h>

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_array.h>
#include <gkyl_array_average.h>
#include <gkyl_array_dg_reduce.h>
#include <gkyl_array_integrate.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_format_desc.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_differentiate.h>
#include <gkyl_dg_eval_at_coord_proj.h>
#include <gkyl_dynvec.h>
#include <gkyl_proj_powsqrt_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

_Static_assert(gpython_MAX_DIM == GKYL_MAX_DIM,
    "gpython_MAX_DIM out of sync with GKYL_MAX_DIM");

int
gpython_api_version(void)
{
  return GPYTHON_API_VERSION;
}

/* The handles are the gkyl objects themselves; the opaque gpython types exist
 * so that no gkyl type is ever named in gpython.h. */
#define ARR(a) ((struct gkyl_array *)(a))
#define CARR(a) ((const struct gkyl_array *)(a))
#define BAS(b) ((struct gkyl_basis *)(b))
#define CBAS(b) ((const struct gkyl_basis *)(b))

/* ---- arrays ------------------------------------------------------------ */
gpython_array *
gpython_array_new(size_t ncomp, size_t size)
{
  return (gpython_array *)gkyl_array_new(GKYL_DOUBLE, ncomp, size);
}

gpython_array *
gpython_array_from_buff(size_t ncomp, size_t size, double *data)
{
  return (gpython_array *)gkyl_array_new_from_buff(GKYL_DOUBLE, ncomp, size, data);
}

gpython_array *
gpython_array_clone(const gpython_array *a)
{
  return (gpython_array *)gkyl_array_clone(CARR(a));
}

void
gpython_array_release(gpython_array *a)
{
  gkyl_array_release(CARR(a));
}

size_t
gpython_array_ncomp(const gpython_array *a)
{
  return CARR(a)->ncomp;
}

size_t
gpython_array_size(const gpython_array *a)
{
  return CARR(a)->size;
}

double *
gpython_array_data(gpython_array *a)
{
  return (double *)ARR(a)->data;
}

/* ---- file I/O ----------------------------------------------------------- */
static void
grid_out(const struct gkyl_rect_grid *grid, int *ndim, double *lower,
    double *upper, int *cells)
{
  *ndim = grid->ndim;
  for (int d = 0; d < grid->ndim; ++d) {
    lower[d] = grid->lower[d];
    upper[d] = grid->upper[d];
    cells[d] = grid->cells[d];
  }
}

int
gpython_file_type(const char *fname)
{
  return gkyl_get_gkyl_file_type(fname);
}

int
gpython_read_header(const char *fname, int *ndim, double *lower, double *upper,
    int *cells, int *file_type, size_t *esznc, size_t *tot_cells, char **meta,
    size_t *meta_sz)
{
  struct gkyl_rect_grid grid;
  struct gkyl_array_header_info hdr;
  enum gkyl_array_rio_status status =
      gkyl_grid_sub_array_header_read(&grid, &hdr, fname);
  if (status != GKYL_ARRAY_RIO_SUCCESS)
    return (int)status;

  grid_out(&grid, ndim, lower, upper, cells);
  *file_type = (int)hdr.file_type;
  *esznc = (size_t)hdr.esznc;
  *tot_cells = (size_t)hdr.tot_cells;
  *meta = NULL;
  *meta_sz = (size_t)hdr.meta_size;
  if (hdr.meta_size > 0) {
    *meta = malloc(hdr.meta_size);
    memcpy(*meta, hdr.meta, hdr.meta_size);
  }
  gkyl_array_header_info_release(&hdr);
  return 0;
}

void
gpython_meta_release(char *meta)
{
  free(meta);
}

gpython_array *
gpython_read_field(const char *fname, int *ndim, double *lower, double *upper,
    int *cells)
{
  struct gkyl_rect_grid grid;
  struct gkyl_array *arr = gkyl_grid_array_new_from_file(&grid, fname);
  if (!arr)
    return NULL;
  grid_out(&grid, ndim, lower, upper, cells);
  return (gpython_array *)arr;
}

const char *
gpython_status_msg(int status)
{
  return gkyl_array_rio_status_msg((enum gkyl_array_rio_status)status);
}

/* ---- basis --------------------------------------------------------------
 * The heap-allocated struct gkyl_basis is plain data filled by the in-place
 * initializers; the function pointers it carries are dispatched HERE, in
 * compiled code, never from the interpreter. */
gpython_basis *
gpython_basis_new(const char *type, int ndim, int poly_order)
{
  struct gkyl_basis *b = malloc(sizeof(struct gkyl_basis));
  if (strcmp(type, "serendipity") == 0)
    gkyl_cart_modal_serendip(b, ndim, poly_order);
  else if (strcmp(type, "tensor") == 0)
    gkyl_cart_modal_tensor(b, ndim, poly_order);
  else {
    free(b);
    return NULL;
  }
  return (gpython_basis *)b;
}

gpython_basis *
gpython_basis_new_hybrid(const char *type, int cdim, int vdim)
{
  struct gkyl_basis *b = malloc(sizeof(struct gkyl_basis));
  if (strcmp(type, "hybrid") == 0)
    gkyl_cart_modal_hybrid(b, cdim, vdim);
  else if (strcmp(type, "gkhybrid") == 0)
    gkyl_cart_modal_gkhybrid(b, cdim, vdim);
  else {
    free(b);
    return NULL;
  }
  return (gpython_basis *)b;
}

void
gpython_basis_release(gpython_basis *b)
{
  free(b);
}

int
gpython_basis_ndim(const gpython_basis *b)
{
  return (int)CBAS(b)->ndim;
}

int
gpython_basis_poly_order(const gpython_basis *b)
{
  return (int)CBAS(b)->poly_order;
}

int
gpython_basis_num_basis(const gpython_basis *b)
{
  return (int)CBAS(b)->num_basis;
}

const char *
gpython_basis_id(const gpython_basis *b)
{
  return CBAS(b)->id;
}

void
gpython_basis_eval(const gpython_basis *b, const double *z, double *bvals)
{
  CBAS(b)->eval(z, bvals);
}

void
gpython_basis_node_list(const gpython_basis *b, double *coords)
{
  CBAS(b)->node_list(coords);
}

void
gpython_basis_nodal_to_modal(const gpython_basis *b, const double *fnodal,
    double *fmodal)
{
  CBAS(b)->nodal_to_modal(fnodal, fmodal);
}

/* ---- weak (DG) algebra --------------------------------------------------
 * Arrays hold nfields = ncomp/num_basis DG fields side by side; the c_oop
 * arguments of the gkyl entry points are FIELD indices. The basis struct is
 * passed BY VALUE, exactly as the header declares — the compiler, not a
 * hand-written mirror, guarantees the ABI. */
static int
nfields_of(const gpython_basis *b, const gpython_array *a)
{
  size_t nb = (size_t)CBAS(b)->num_basis;
  if (CARR(a)->ncomp % nb != 0)
    return -1;
  return (int)(CARR(a)->ncomp / nb);
}

static int
weak_shapes_ok(const gpython_array *out, const gpython_array *a1, const gpython_array *a2)
{
  return CARR(out)->ncomp == CARR(a1)->ncomp &&
         CARR(out)->size == CARR(a1)->size &&
         (!a2 || (CARR(a1)->ncomp == CARR(a2)->ncomp &&
                  CARR(a1)->size == CARR(a2)->size));
}

int
gpython_dg_mul(const gpython_basis *b, gpython_array *out, const gpython_array *a1,
    const gpython_array *a2)
{
  int nf = nfields_of(b, a1);
  if (nf < 0 || !weak_shapes_ok(out, a1, a2))
    return 1;
  for (int f = 0; f < nf; ++f)
    gkyl_dg_mul_op(CBAS(b), f, ARR(out), f, CARR(a1), f, CARR(a2));
  return 0;
}

int
gpython_dg_div(const gpython_basis *b, gpython_array *out, const gpython_array *a1,
    const gpython_array *a2)
{
  int nf = nfields_of(b, a1);
  if (nf < 0 || !weak_shapes_ok(out, a1, a2))
    return 1;
  gkyl_dg_bin_op_mem *mem =
      gkyl_dg_bin_op_mem_new(CARR(a1)->size, CBAS(b)->num_basis);
  for (int f = 0; f < nf; ++f)
    gkyl_dg_div_op(mem, CBAS(b), f, ARR(out), f, CARR(a1), f, CARR(a2));
  gkyl_dg_bin_op_mem_release(mem);
  return 0;
}

int
gpython_dg_inv(const gpython_basis *b, gpython_array *out, const gpython_array *a1)
{
  int nf = nfields_of(b, a1);
  if (nf < 0 || !weak_shapes_ok(out, a1, NULL))
    return 1;
  for (int f = 0; f < nf; ++f)
    gkyl_dg_inv_op(CBAS(b), f, ARR(out), f, CARR(a1));
  return 0;
}

/* conf*phase: single-field only, so shapes are just num_basis equality plus
 * the range volumes built from the caller's cell counts matching the
 * arrays' sizes exactly (mirrors gpython_array_integrate's range convention:
 * 1-indexed, lower=1, upper=cells). */
int
gpython_dg_mul_conf_phase(const gpython_basis *cbasis, const gpython_basis *pbasis,
    gpython_array *pout, const gpython_array *cop, const gpython_array *pop,
    const int *conf_cells, const int *phase_cells)
{
  int cdim = CBAS(cbasis)->ndim, pdim = CBAS(pbasis)->ndim;
  size_t cnb = (size_t)CBAS(cbasis)->num_basis, pnb = (size_t)CBAS(pbasis)->num_basis;
  if (cnb >= pnb || CARR(cop)->ncomp != cnb || CARR(pop)->ncomp != pnb ||
      CARR(pout)->ncomp != pnb || CARR(pout)->size != CARR(pop)->size)
    return 1;

  int ones[GKYL_MAX_DIM];
  for (int d = 0; d < pdim; ++d) ones[d] = 1;

  struct gkyl_range crange, prange;
  gkyl_range_init(&crange, cdim, ones, conf_cells);
  gkyl_range_init(&prange, pdim, ones, phase_cells);
  if ((size_t)crange.volume != CARR(cop)->size ||
      (size_t)prange.volume != CARR(pop)->size)
    return 1;

  gkyl_dg_mul_conf_phase_op_range(CBAS(cbasis), CBAS(pbasis), ARR(pout),
      CARR(cop), CARR(pop), &crange, &prange);
  return 0;
}

/* Local DG derivative: like gpython_dg_mul, the per-field loop lives here since
 * gkyl_dg_differentiate_op_local's c_oop/c_iop are field indices, not a
 * whole-array operation. */
int
gpython_dg_differentiate(const gpython_basis *b, int dir, int diff_order, double dx,
    gpython_array *out, const gpython_array *in)
{
  int nf = nfields_of(b, in);
  if (nf < 0 || !weak_shapes_ok(out, in, NULL))
    return 1;
  int ndim = CBAS(b)->ndim;
  if (dir < 0 || dir >= ndim || diff_order < 1 || diff_order > 2)
    return 1;
  for (int f = 0; f < nf; ++f)
    gkyl_dg_differentiate_op_local(CBAS(b), dir, diff_order, dx, f, ARR(out),
        f, CARR(in));
  return 0;
}

/* ---- linear coefficient ops / reductions -------------------------------- */
void
gpython_array_set(gpython_array *out, double c, const gpython_array *a)
{
  gkyl_array_set(ARR(out), c, CARR(a));
}

void
gpython_array_accumulate(gpython_array *out, double c, const gpython_array *a)
{
  gkyl_array_accumulate(ARR(out), c, CARR(a));
}

void
gpython_array_scale(gpython_array *a, double c)
{
  gkyl_array_scale(ARR(a), c);
}

void
gpython_array_shiftc(gpython_array *a, double val, unsigned comp)
{
  gkyl_array_shiftc(ARR(a), val, comp);
}

void
gpython_array_reduce(double *out, const gpython_array *a, int op)
{
  static const enum gkyl_array_op ops[] = { GKYL_MIN, GKYL_MAX, GKYL_SUM };
  gkyl_array_reduce(out, CARR(a), ops[op]);
}

int
gpython_array_dg_reduce(double *out, const gpython_basis *b, const gpython_array *a,
    int comp, int op)
{
  size_t nb = (size_t)CBAS(b)->num_basis;
  if (comp < 0 || CARR(a)->ncomp < (size_t)(comp + 1) * nb)
    return 1;
  static const enum gkyl_array_op ops[] = { GKYL_MIN, GKYL_MAX, GKYL_SUM };
  gkyl_array_dg_reducec(out, CARR(a), comp, ops[op], CBAS(b));
  return 0;
}

/* ---- integration --------------------------------------------------------- */
int
gpython_array_integrate(int ndim, const double *lower, const double *upper,
    const int *cells, const gpython_basis *b, int nfields, int op, double factor,
    const gpython_array *a, double *out)
{
  static const enum gkyl_array_integrate_op ops[] = {
    GKYL_ARRAY_INTEGRATE_OP_NONE,
    GKYL_ARRAY_INTEGRATE_OP_ABS,
    GKYL_ARRAY_INTEGRATE_OP_SQ,
  };

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int ones[GKYL_MAX_DIM];
  for (int d = 0; d < ndim; ++d)
    ones[d] = 1;
  struct gkyl_range range;
  gkyl_range_init(&range, ndim, ones, cells);
  if ((size_t)range.volume != CARR(a)->size)
    return 1;

  struct gkyl_array_integrate *up =
      gkyl_array_integrate_new(&grid, CBAS(b), nfields, ops[op], false);
  gkyl_array_integrate_advance(up, CARR(a), factor, NULL, &range, NULL, out);
  gkyl_array_integrate_release(up);
  return 0;
}

/* ---- averaging ------------------------------------------------------------
 * Single-field only (like gpython_dg_mul_conf_phase): `a`/`weight`/`out` must
 * each carry exactly one basis's worth of coefficients per cell -- the
 * kernel gkyl_array_average_choose_kernel dispatches has no field-index
 * parameter at all, so a multi-field caller must loop, one field at a time,
 * in Python (dg/modal.py). */
int
gpython_array_average(int ndim, const double *lower, const double *upper,
    const int *cells, const gpython_basis *b, const gpython_basis *b_avg,
    int ndim_avg, const int *cells_avg, const int *avg_dim,
    const gpython_array *weight, const gpython_array *a, gpython_array *out)
{
  size_t nb = (size_t)CBAS(b)->num_basis, nb_avg = (size_t)CBAS(b_avg)->num_basis;
  if (CARR(a)->ncomp != nb || CARR(out)->ncomp != nb_avg ||
      (weight && (CARR(weight)->ncomp != nb || CARR(weight)->size != CARR(a)->size)))
    return 1;

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int ones[GKYL_MAX_DIM];
  for (int d = 0; d < ndim; ++d)
    ones[d] = 1;
  struct gkyl_range range;
  gkyl_range_init(&range, ndim, ones, cells);
  if ((size_t)range.volume != CARR(a)->size)
    return 1;

  int ones_avg[GKYL_MAX_DIM];
  for (int d = 0; d < ndim_avg; ++d)
    ones_avg[d] = 1;
  struct gkyl_range range_avg;
  gkyl_range_init(&range_avg, ndim_avg, ones_avg, cells_avg);
  if ((size_t)range_avg.volume != CARR(out)->size)
    return 1;

  /* range_avg doubles as local_avg_ext (mirrors gpython_dg_mul_conf_phase's
   * range convention: 1-indexed, lower=1, upper=cells -- there is no ghost
   * layer anywhere in the gpython boundary, so the "extended" reduced range
   * gkyl_array_average_new only reads to size the integrated weight is
   * identical to the reduced range itself). */
  struct gkyl_array_average *up = gkyl_array_average_new(&grid, CBAS(b),
      CBAS(b_avg), &range, &range_avg, &range_avg,
      weight ? CARR(weight) : NULL, avg_dim, false);
  gkyl_array_average_advance(up, CARR(a), ARR(out));
  gkyl_array_average_release(up);
  return 0;
}

/* ---- evaluate-and-project ------------------------------------------------- */
gpython_array *
gpython_eval_at_coord_proj(const gpython_basis *b, int cdim_do, int ndim,
    const double *lower, const double *upper, const int *cells,
    int num_eval, const int *eval_dirs, const double *eval_coords,
    int ndim_tar, const int *cells_tar, const gpython_array *in,
    int *out_btype, int *out_poly_order, int *out_cdim, int *out_vdim)
{
  int nf = nfields_of(b, in);
  if (nf < 0 || num_eval < 1 || num_eval > ndim)
    return NULL;

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int ones[GKYL_MAX_DIM];
  for (int d = 0; d < ndim; ++d)
    ones[d] = 1;
  struct gkyl_range rng_do;
  gkyl_range_init(&rng_do, ndim, ones, cells);
  if ((size_t)rng_do.volume != CARR(in)->size)
    return NULL;

  int ones_tar[GKYL_MAX_DIM];
  for (int d = 0; d < ndim_tar; ++d)
    ones_tar[d] = 1;
  struct gkyl_range rng_tar;
  gkyl_range_init(&rng_tar, ndim_tar, ones_tar, cells_tar);

  struct gkyl_dg_eval_at_coord_proj *up =
      gkyl_dg_eval_at_coord_proj_new(cdim_do, CBAS(b), num_eval, eval_dirs,
          false);

  enum gkyl_basis_type btype_tar;
  int cdim_tar, ndim_tar_full, poly_order_tar, num_basis_tar;
  gkyl_dg_eval_at_coord_proj_target_basis(up, &cdim_tar, &ndim_tar_full,
      &btype_tar, &poly_order_tar, &num_basis_tar);

  struct gkyl_array *out = gkyl_array_new(GKYL_DOUBLE,
      (size_t)nf * (size_t)num_basis_tar, (size_t)rng_tar.volume);

  bool pick_lower[GKYL_MAX_DIM];
  int known_index[GKYL_MAX_DIM];
  for (int d = 0; d < ndim; ++d)
    known_index[d] = -1;
  for (int i = 0; i < num_eval; ++i)
    pick_lower[i] = false;

  gkyl_dg_eval_at_coord_proj_advance(up, eval_coords, &grid, pick_lower,
      known_index, &rng_do, &rng_tar, CARR(in), out);
  gkyl_dg_eval_at_coord_proj_release(up);

  *out_btype = (int)btype_tar;
  *out_poly_order = poly_order_tar;
  *out_cdim = cdim_tar;
  *out_vdim = ndim_tar_full - cdim_tar;

  return (gpython_array *)out;
}

/* ---- pow(sqrt) projection -------------------------------------------------- */
int
gpython_powsqrt(const gpython_basis *b, int num_quad, double exponent,
    int ndim, const int *cells, gpython_array *out, const gpython_array *in)
{
  size_t nb = (size_t)CBAS(b)->num_basis;
  if (CARR(in)->ncomp != nb || CARR(out)->ncomp != nb)
    return 1;

  int ones[GKYL_MAX_DIM];
  for (int d = 0; d < ndim; ++d)
    ones[d] = 1;
  struct gkyl_range range;
  gkyl_range_init(&range, ndim, ones, cells);
  if ((size_t)range.volume != CARR(in)->size || CARR(out)->size != CARR(in)->size)
    return 1;

  struct gkyl_proj_powsqrt_on_basis *up =
      gkyl_proj_powsqrt_on_basis_new(CBAS(b), num_quad, false);
  gkyl_proj_powsqrt_on_basis_advance(up, &range, exponent, CARR(in), ARR(out));
  gkyl_proj_powsqrt_on_basis_release(up);
  return 0;
}

/* ---- writing -------------------------------------------------------------- */
int
gpython_write_field(const char *fname, int ndim, const double *lower,
    const double *upper, const int *cells, const char *meta, size_t meta_sz,
    const gpython_array *a)
{
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int ones[GKYL_MAX_DIM];
  for (int d = 0; d < ndim; ++d)
    ones[d] = 1;
  struct gkyl_range range;
  gkyl_range_init(&range, ndim, ones, cells);
  if ((size_t)range.volume != CARR(a)->size)
    return -1; /* sentinel: never a valid gkyl_array_rio_status value */

  struct gkyl_msgpack_data md = { .meta_sz = meta_sz, .meta = (char *)meta };
  return (int)gkyl_grid_sub_array_write(&grid, &range, meta_sz ? &md : NULL,
      CARR(a), fname);
}

/* ---- dynvector (time-series) I/O ------------------------------------------
 * The dynvec object itself never crosses gkyl_gpython.h: it is created, filled,
 * and released entirely inside this function, and its contents move to the
 * caller only as gpython_arrays (already opaque, already RAII'd). */
int
gpython_dynvec_read(const char *fname, size_t *ncomp, gpython_array **tm,
    gpython_array **data)
{
  struct gkyl_dynvec_etype_ncomp info = gkyl_dynvec_read_ncomp(fname);
  if (info.ncomp == 0)
    return 1; /* missing file or empty/unrecognized header */
  if (info.type != GKYL_DOUBLE)
    return 2; /* postgkyl only ever writes/reads double dynvectors */

  gkyl_dynvec vec = gkyl_dynvec_new(GKYL_DOUBLE, info.ncomp);
  bool ok = gkyl_dynvec_read(vec, fname);
  if (!ok) {
    gkyl_dynvec_release(vec);
    return 3;
  }

  size_t n = gkyl_dynvec_size(vec);
  struct gkyl_array *tm_arr = gkyl_array_new(GKYL_DOUBLE, 1, n);
  struct gkyl_array *data_arr = gkyl_array_new(GKYL_DOUBLE, info.ncomp, n);
  gkyl_dynvec_to_array(vec, tm_arr, data_arr);
  gkyl_dynvec_release(vec);

  *ncomp = info.ncomp;
  *tm = (gpython_array *)tm_arr;
  *data = (gpython_array *)data_arr;
  return 0;
}

int
gpython_dynvec_write(const char *fname, size_t ncomp, size_t n, const double *tm,
    const double *data)
{
  gkyl_dynvec vec = gkyl_dynvec_new(GKYL_DOUBLE, ncomp);
  for (size_t i = 0; i < n; ++i)
    gkyl_dynvec_append(vec, tm[i], data + i * ncomp);
  int status = gkyl_dynvec_write(vec, fname);
  gkyl_dynvec_release(vec);
  return status;
}
