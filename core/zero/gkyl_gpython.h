/* gkyl_gpython.h — the postgkyl support API (the "gpython shim").
 *
 * A small, stable C surface over the core library for the postgkyl
 * post-processing tool. It deliberately exposes NO gkyl types: only opaque
 * handles, scalars, and plain buffers, so a consumer compiled against this
 * header holds no knowledge of struct layouts or calling conventions.
 * Everything layout- and ABI-sensitive lives in zero/gpython.c, compiled into
 * libg0core.so — so shim and library can never drift apart, and any core
 * API change that affects the shim fails THIS build, at the producer.
 *
 * Consumers (postgkyl's _gpython CPython extension) include only this header
 * and link -lg0core. Bump GPYTHON_API_VERSION when this header changes shape;
 * consumers check it at load time to reject a stale pairing.
 */
#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Bump when this header changes shape; checked by consumers at load time.
 * v2: added gpython_array_dg_reduce, gpython_write_field, gpython_dynvec_{read,write}.
 * v3: added gpython_basis_new_hybrid (hybrid/gkhybrid bases).
 * v4: added gpython_dg_mul_conf_phase (conf-space x phase-space weak product).
 * v5: added gpython_array_average (weighted/plain average over a dim subset).
 * v6: added gpython_dg_differentiate (local DG derivative) and
 *     gpython_eval_at_coord_proj (evaluate at coords + project onto a
 *     lower-dimensional target basis).
 * v7: added gpython_powsqrt (pow(sqrt(.), exponent) via
 *     gkyl_proj_powsqrt_on_basis). */
#define GPYTHON_API_VERSION 7
int gpython_api_version(void);

/* Enough for GKYL_MAX_DIM; fixed here so callers can size buffers without
 * seeing gkyl headers. gpython.c statically asserts it matches GKYL_MAX_DIM. */
#define gpython_MAX_DIM 7

typedef struct gpython_array gpython_array; /* wraps a gkyl_array* (double data)   */
typedef struct gpython_basis gpython_basis; /* wraps a struct gkyl_basis           */

/* ---- arrays: create / copy / introspect / release ---------------------- */
gpython_array *gpython_array_new(size_t ncomp, size_t size); /* zero-initialized   */
/* Zero-copy view of caller-owned memory laid out (size x ncomp) doubles.
 * The buffer must outlive the returned handle. */
gpython_array *gpython_array_from_buff(size_t ncomp, size_t size, double *data);
gpython_array *gpython_array_clone(const gpython_array *a);
void gpython_array_release(gpython_array *a);
size_t gpython_array_ncomp(const gpython_array *a);
size_t gpython_array_size(const gpython_array *a);
double *gpython_array_data(gpython_array *a); /* contiguous (size x ncomp) buffer  */

/* ---- file I/O (gkyl_array_rio) ------------------------------------------
 * Grid out-parameters are caller-allocated buffers of length gpython_MAX_DIM.
 * Status codes are gkyl_array_rio_status values; 0 == success.            */
int gpython_file_type(const char *fname); /* gkyl file type, or -1 */
int gpython_read_header(const char *fname, int *ndim, double *lower,
    double *upper, int *cells, int *file_type, size_t *esznc,
    size_t *tot_cells, char **meta, size_t *meta_sz); /* *meta: malloc'd   */
void gpython_meta_release(char *meta);
gpython_array *gpython_read_field(const char *fname, int *ndim, double *lower,
    double *upper, int *cells); /* NULL on failure */
const char *gpython_status_msg(int status);

/* ---- basis objects ------------------------------------------------------
 * type: "serendipity" | "tensor" (the long names stored in .gkyl metadata).
 * Returns NULL for unknown types. Use gpython_basis_new_hybrid for the two
 * mixed-order bases below.                                                 */
gpython_basis *gpython_basis_new(const char *type, int ndim, int poly_order);
/* type: "hybrid" | "gkhybrid" — mixed-order bases with a fixed poly_order
 * of 1, parameterized by configuration/velocity dimension (ndim = cdim +
 * vdim) rather than by an explicit poly_order, mirroring Gkeyll's own
 * gkyl_cart_modal_hybrid/gkhybrid constructors exactly. Returns NULL for
 * unknown types; out-of-range cdim/vdim is a process abort in the
 * underlying constructor (see gkyl_basis.h), not a clean failure — callers
 * must respect its bounds.                                                 */
gpython_basis *gpython_basis_new_hybrid(const char *type, int cdim, int vdim);
void gpython_basis_release(gpython_basis *b);
int gpython_basis_ndim(const gpython_basis *b);
int gpython_basis_poly_order(const gpython_basis *b);
int gpython_basis_num_basis(const gpython_basis *b);
const char *gpython_basis_id(const gpython_basis *b);
/* Evaluate all num_basis functions at reference point z[ndim].            */
void gpython_basis_eval(const gpython_basis *b, const double *z, double *bvals);
/* Node coordinates: coords[num_basis * ndim].                              */
void gpython_basis_node_list(const gpython_basis *b, double *coords);
/* Exact nodal -> modal change of basis on one cell's num_basis values.     */
void gpython_basis_nodal_to_modal(const gpython_basis *b, const double *fnodal,
    double *fmodal);

/* ---- weak (DG) algebra --------------------------------------------------
 * Operands must have ncomp == nfields * num_basis; the per-field loop runs
 * here. Returns 0 on success, nonzero on shape mismatch.                   */
int gpython_dg_mul(const gpython_basis *b, gpython_array *out, const gpython_array *a1,
    const gpython_array *a2);
int gpython_dg_div(const gpython_basis *b, gpython_array *out, const gpython_array *a1,
    const gpython_array *a2);
int gpython_dg_inv(const gpython_basis *b, gpython_array *out, const gpython_array *a1);
/* Conf-space x phase-space weak product: pout = cop * pop, where cop lives
 * on a conf-space grid of cbasis->ndim dimensions and pop/pout live on a
 * phase-space grid of pbasis->ndim (> cbasis->ndim) dimensions sharing the
 * conf-space grid's leading dimensions (gkyl_dg_mul_conf_phase_op_range).
 * Single-field only: cop.ncomp == cbasis num_basis, pop/pout.ncomp == pbasis
 * num_basis -- the underlying kernel takes no field-index arguments.
 * conf_cells/phase_cells give each grid's per-dimension cell count (length
 * cbasis->ndim / pbasis->ndim), used to build the index ranges Gkeyll maps
 * each phase cell's conf cell through. Returns 0 on success, nonzero if the
 * arrays don't match the bases or cells don't cover the arrays.             */
int gpython_dg_mul_conf_phase(const gpython_basis *cbasis, const gpython_basis *pbasis,
    gpython_array *pout, const gpython_array *cop, const gpython_array *pop,
    const int *conf_cells, const int *phase_cells);

/* Local DG derivative (gkyl_dg_differentiate_op_local): differentiates the
 * DG expansion independently in every cell (no inter-cell stencil), field
 * by field -- serendipity/tensor only (a Gkeyll limit; hybrid/gkhybrid have
 * no differentiate kernels and hit an unconditional `assert(false)` in the
 * underlying dispatch, so this must never be called with those bases).
 * `dir` is 0-based, `diff_order` is 1 or 2. `out`/`in` must each carry
 * ncomp == nfields * num_basis coefficients (the same field-count
 * convention as gpython_dg_mul); `dx` is the cell length in `dir`. Returns 0 on
 * success, nonzero on a field-count/shape mismatch or an out-of-range
 * `dir`/`diff_order` -- the caller (postgkyl's dg/modal.py) is responsible
 * for staying within Gkeyll's compiled kernel-table bounds (ndim/poly_order
 * per basis_type), which this function does NOT check: an out-of-table
 * combination is a process abort in the kernel dispatch, not a clean
 * failure.                                                                 */
int gpython_dg_differentiate(const gpython_basis *b, int dir, int diff_order,
    double dx, gpython_array *out, const gpython_array *in);

/* ---- linear coefficient ops / reductions ------------------------------- */
void gpython_array_set(gpython_array *out, double c, const gpython_array *a);
void gpython_array_accumulate(gpython_array *out, double c, const gpython_array *a);
void gpython_array_scale(gpython_array *a, double c);
void gpython_array_shiftc(gpython_array *a, double val, unsigned comp);
/* op: 0 min, 1 max, 2 sum; out[ncomp]. Reduces the raw DG COEFFICIENTS —
 * exact for "sum" but not the field's true min/max; see
 * gpython_array_dg_reduce for the field-aware version.                        */
void gpython_array_reduce(double *out, const gpython_array *a, int op);

/* Field-aware reduction (gkyl_array_dg_reducec): evaluates the DG expansion
 * of field `comp` (0-based, ncomp must be a multiple of num_basis) at each
 * cell's Gauss-Legendre quadrature nodes and reduces THOSE — the true
 * min/max/sum of the represented field, not the coefficients. op: 0/1/2
 * min/max/sum. Returns 0 on success, nonzero if comp is out of range.      */
int gpython_array_dg_reduce(double *out, const gpython_basis *b, const gpython_array *a,
    int comp, int op);

/* ---- integration (gkyl_array_integrate) ---------------------------------
 * op: 0 none, 1 abs, 2 sq. out[nfields]. Returns 0 on success, nonzero if
 * the grid does not cover the array.                                        */
int gpython_array_integrate(int ndim, const double *lower, const double *upper,
    const int *cells, const gpython_basis *b, int nfields, int op, double factor,
    const gpython_array *a, double *out);

/* ---- averaging (gkyl_array_average) -------------------------------------
 * Single-field weighted (or plain) average of `a` over the donor dims
 * flagged in `avg_dim` (length ndim, 1 = averaged, 0 = kept):
 *   weight != NULL:  int f(x) w(x) dx^avg / int w(x) dx^avg
 *   weight == NULL:  int f(x) dx^avg / int dx^avg
 * `b` is the full donor basis (ndim/poly_order); `b_avg` the reduced-
 * dimension (same poly_order) target basis for the surviving `ndim_avg`
 * dims -- `ndim_avg` is 1 with `cells_avg = {1}` for a full reduction
 * (Gkeyll keeps a single-cell target rather than a true 0-dim basis, same
 * convention `gpython_array_integrate`'s caller uses at the Python layer).
 * `a`/`weight`/`out` are single-field: their ncomp must equal the relevant
 * basis's num_basis exactly (gkyl_array_average takes no field-index
 * argument, like gpython_dg_mul_conf_phase; a multi-field caller loops in
 * Python). Returns 0 on success, nonzero if `a`/`weight`/`out` don't match
 * the bases or the donor/target cell counts don't cover them.               */
int gpython_array_average(int ndim, const double *lower, const double *upper,
    const int *cells, const gpython_basis *b, const gpython_basis *b_avg,
    int ndim_avg, const int *cells_avg, const int *avg_dim,
    const gpython_array *weight, const gpython_array *a, gpython_array *out);

/* ---- evaluate-and-project (gkyl_dg_eval_at_coord_proj) --------------------
 * Evaluate a donor DG field at physical coordinates `eval_coords` in the
 * `num_eval` directions `eval_dirs` (0-based, distinct, < ndim), and project
 * the result onto the lower-dimensional target basis Gkeyll picks for that
 * elimination -- e.g. a 3-D field evaluated at one z1 coordinate becomes a
 * 2-D field over the surviving (z0, z2) directions. Allocates and returns
 * the target array (ncomp = nfields * num_basis_tar, nfields = a.ncomp /
 * donor num_basis, same field-count convention as gpython_dg_mul); NULL on a
 * shape mismatch or invalid eval_dirs.
 *
 * `cdim_do` is the donor's configuration-space dimension count (equal to
 * `ndim` for pure serendipity/tensor fields; less than `ndim` for a
 * gkhybrid phase-space field, where it marks the conf/velocity split --
 * see gkeyll's own gkyl_cart_modal_gkhybrid convention). `ndim_tar`/
 * `cells_tar` describe the target's rectangular index range with the same
 * convention gpython_array_average's ndim_avg/cells_avg use: the surviving
 * donor dims' cell counts in donor order, or ndim_tar=1/cells_tar={1} when
 * every donor direction is evaluated away (Gkeyll always keeps at least one
 * target dimension; there is no true 0-dim basis).
 *
 * The target basis can differ in TYPE from the donor (e.g. eliminating a
 * gkhybrid velocity direction can yield a plain serendipity target), so its
 * metadata is queried from the updater only after construction and reported
 * via the out-params:
 *   *out_btype      -- gkyl_basis_type ordinal (gkyl_basis.h); translate to
 *                       the pgkyl string vocabulary in Python
 *   *out_poly_order, *out_cdim, *out_vdim -- target basis parameters
 *     (out_vdim = target ndim - out_cdim)
 *
 * Guarded to Gkeyll's compiled kernel-table coverage the same way
 * gpython_dg_differentiate is: an out-of-table (basis_type, ndim, poly_order,
 * eval_dirs) combination is a process abort in the kernel dispatch, not a
 * clean failure, so the caller (postgkyl's dg/modal.py) must stay within
 * the documented coverage before calling this function.                   */
gpython_array *gpython_eval_at_coord_proj(const gpython_basis *b, int cdim_do, int ndim,
    const double *lower, const double *upper, const int *cells,
    int num_eval, const int *eval_dirs, const double *eval_coords,
    int ndim_tar, const int *cells_tar, const gpython_array *in,
    int *out_btype, int *out_poly_order, int *out_cdim, int *out_vdim);

/* ---- pow(sqrt) projection (gkyl_proj_powsqrt_on_basis) -------------------
 * Single-field ``out = pow(sqrt(in), exponent)``, i.e. ``in ** (exponent/2)``,
 * via Gauss-Legendre quadrature: evaluate at each cell's quadrature nodes
 * (clamping a negative node value to 1e-40, the updater's own convention,
 * not checked against the DG coefficients) and project back onto the basis.
 * Unlike every weak (DG) op above, this is NOT a fixed per-(basis_type,
 * ndim, poly_order) compiled kernel table -- the updater works off `basis`'s
 * own `eval` callback, so it has no coverage restriction at all. `num_quad`
 * is the number of quadrature nodes per dimension (Gkeyll's own gyrokinetic
 * app convention: poly_order + 1). `out`/`in` must each carry exactly one
 * basis's worth of coefficients (ncomp == b->num_basis) -- the underlying
 * updater takes no field-index argument, like gpython_dg_mul_conf_phase; a
 * multi-field caller loops in Python. No physical grid is needed (only cell
 * indexing), unlike gpython_array_integrate/gpython_array_average, so there
 * is no lower/upper here. Returns 0 on success, nonzero if the operands
 * don't match the basis or `cells` does not cover them. */
int gpython_powsqrt(const gpython_basis *b, int num_quad, double exponent,
    int ndim, const int *cells, gpython_array *out, const gpython_array *in);

/* ---- writing (gkyl_array_rio) --------------------------------------------
 * Mirrors gpython_read_field: writes the FULL array over a uniform grid built
 * from lower/upper/cells (the same convention gpython_array_integrate uses).
 * meta/meta_sz: a raw msgpack byte blob, or meta_sz == 0 for none — encoding
 * policy lives in postgkyl's io/ layer, exactly as decoding does for reads.
 * Returns -1 if the grid does not cover the array (checked here to avoid an
 * out-of-bounds C loop); otherwise a gkyl_array_rio_status (0 == success,
 * rendered by gpython_status_msg).                                             */
int gpython_write_field(const char *fname, int ndim, const double *lower,
    const double *upper, const int *cells, const char *meta, size_t meta_sz,
    const gpython_array *a);

/* ---- dynvector (time-series) I/O (gkyl_dynvec) ---------------------------
 * Double-precision dynvectors only (the only kind postgkyl ever produces).
 * gpython_dynvec_read allocates *tm (ncomp=1) and *data (ncomp=*ncomp), both
 * sized to the number of time samples found in the file — the caller
 * releases them like any other gpython_array. Returns 0 on success; nonzero if
 * the file is missing/unreadable or not a double-typed dynvector.          */
int gpython_dynvec_read(const char *fname, size_t *ncomp, gpython_array **tm,
    gpython_array **data);
/* Writes a dynvector built by zipping tm[n] with data[n*ncomp:(n+1)*ncomp].
 * Returns 0 on success.                                                    */
int gpython_dynvec_write(const char *fname, size_t ncomp, size_t n,
    const double *tm, const double *data);

#ifdef __cplusplus
}
#endif
