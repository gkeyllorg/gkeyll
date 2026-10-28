/** 
 *
 * NVECTOR wrapper of the Gkeyll data structure.
 *
 */

#include <gkyl_sundials.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <string.h> // for memcpy.

#ifdef GKYL_HAVE_SUNDIALS

#include <gkyl_sundials_priv.h>

/**
 * Create a new Gkeyll array.
 *
 * @param use_gpu Whether to allocate it on the GPU.
 * @param nc Number of components per cell.
 * @param size Number of cells.
 * @return A new gkyl_array pointer.
 */
static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (use_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

/**
 * Make an empty copy of a given SUNDIALS NVECTOR.
 *
 * @param nvin An input NVECTOR.
 * @return A new nvector.
 */
static N_Vector
snvec_clone_empty(N_Vector nvin)
{
  N_Vector nvout;
  N_VectorContent_Gkeyll content;
  if (nvin == NULL)
    return (NULL);

  // Create vector.
  nvout = NULL;
  nvout = N_VNewEmpty(nvin->sunctx);
  if (nvout == NULL)
    return (NULL);

  // Attach operations.
  if (N_VCopyOps(nvin, nvout)) {
    N_VDestroy(nvout);
    return (NULL);
  }

  // Create content.
  content = NULL;
  content = (N_VectorContent_Gkeyll) gkyl_malloc(sizeof *content);
  if (content == NULL) {
    N_VDestroy(nvout);
    return (NULL);
  }

  // Attach content.
  nvout->content = content;

  // Initialize content.
  content->own_vector = SUNFALSE;
  // gkyl_sundials_nvec_clone_empty must be called only in
  // gkyl_sundials_nvec_clone to ensure use_gpu has the correct flag.
  // Otherwise, use_gpu flag will be false even if it must be true.

  // TO DO: Check to verify if this function is called separately
  content->use_gpu       = SUNFALSE;
  content->arr           = NULL;
  content->comm          = NULL;
  content->local_range   = NULL;
  content->no_error_norm = SUNFALSE;
  content->not_stepped   = SUNFALSE;
  content->red_mem       = NULL;

  return nvout;
}

/**
 * Make an copy of a given SUNDIALS NVECTOR.
 *
 * @param nvin An input NVECTOR.
 * @return A new nvector.
 */
static N_Vector
snvec_clone(N_Vector nvin)
{
  N_Vector nvout;

  struct gkyl_array* nvout_arr;
  struct gkyl_array* nvin_arr = NV_CONTENT_GKZ(nvin)->arr;

  nvout = NULL;
  nvout = snvec_clone_empty(nvin);
  if (nvout == NULL)
    return (NULL);

  nvout_arr = mkarr(NV_CONTENT_GKZ(nvin)->use_gpu, nvin_arr->ncomp, nvin_arr->size);

  NV_CONTENT_GKZ(nvout)->arr           = nvout_arr;
  NV_CONTENT_GKZ(nvout)->use_gpu       = NV_CONTENT_GKZ(nvin)->use_gpu;
  NV_CONTENT_GKZ(nvout)->comm          = NV_CONTENT_GKZ(nvin)->comm;
  NV_CONTENT_GKZ(nvout)->local_range   = NV_CONTENT_GKZ(nvin)->local_range;
  NV_CONTENT_GKZ(nvout)->global_volume = NV_CONTENT_GKZ(nvin)->global_volume;
  NV_CONTENT_GKZ(nvout)->no_error_norm = NV_CONTENT_GKZ(nvin)->no_error_norm;
  NV_CONTENT_GKZ(nvout)->not_stepped   = NV_CONTENT_GKZ(nvin)->not_stepped;
  NV_CONTENT_GKZ(nvout)->red_mem       = NV_CONTENT_GKZ(nvin)->red_mem;
  NV_CONTENT_GKZ(nvout)->own_vector    = SUNTRUE;

  return nvout;
}

/**
 * Destroy SUNDIALS NVECTOR.
 *
 * @param nvin NVECTOR to be destroyed.
 */
static void
snvec_destroy(N_Vector nvin)
{
  if (nvin == 0)
    return;

  // Free content.
  if (nvin->content != 0) {
    // Free the gkeyll parvector if it's owned by the vector wrapper.
    if (NV_CONTENT_GKZ(nvin)->own_vector && NV_CONTENT_GKZ(nvin)->arr != 0) {
      gkyl_array_release(NV_CONTENT_GKZ(nvin)->arr);
      NV_CONTENT_GKZ(nvin)->arr = 0;
    }
    free(nvin->content);
    nvin->content = 0;
  }

  // Free ops and vector.
  if (nvin->ops != 0) {
    free(nvin->ops);
    nvin->ops = 0;
  }
  free(nvin);
  nvin = 0;
}

/**
 * Set some entries 0. Needed by SUNDIALS but may be removed in
 * future SUNDIALS releases.
 *
 * @param v Input Nvector.
 * @param x Input index.
 * @param y Input index.
 */
static void
snvec_space(N_Vector v, sunindextype* x, sunindextype* y)
{
  *x = 0;
  *y = 0;
}

/**
 * Return the number of DOF wrapped by the nvector, excluding ghost cells.
 * future SUNDIALS releases.
 *
 * @param v Input Nvector.
 * @return Number of DOF wrapped by v (excluding ghosts).
 */
static sunindextype
snvec_get_length(N_Vector x)
{
  struct gkyl_array *xarr = NV_CONTENT_GKZ(x)->arr;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(x)->local_range;
  return xarr->ncomp * local_range->volume;
}

/**
 * Perform the linear sum z = a*x + b*y where z, x and y are SUNDIALS NVECTORS,
 * and a and b are scalars.
 *
 * @param a Scalar to multiply `x`.
 * @param x Nvector to multiply by `a` and add to `z`.
 * @param b Scalar to multiply `y`.
 * @param y Nvector to multiply by `b` and add to `z`.
 * @param z Output Nvector.
 */
static void
snvec_linear_sum(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y, N_Vector z)
{
  bool not_stepped = NV_CONTENT_GKZ(z)->not_stepped;
  if (not_stepped)
    return;

  struct gkyl_array *xarr = NV_CONTENT_GKZ(x)->arr;
  struct gkyl_array *yarr = NV_CONTENT_GKZ(y)->arr;
  struct gkyl_array *zarr = NV_CONTENT_GKZ(z)->arr;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(x)->local_range;

  gkyl_array_comp_op_range(zarr, GKYL_AXPBY, a, xarr, b, yarr, local_range);
}

/**
 * Set the whole SUNDIALS NVECTOR to a constant value.
 *
 * @param c Scalar value.
 * @param z Nvector to assign `c` value to.
 */
static void
snvec_const(sunrealtype c, N_Vector z)
{
  bool not_stepped = NV_CONTENT_GKZ(z)->not_stepped;
  if (not_stepped)
    return;

  struct gkyl_array *zarr = NV_CONTENT_GKZ(z)->arr;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(z)->local_range;

  gkyl_array_clear_range(zarr, c, local_range);
}

/**
 * Scale a given SUNDIALS NVECTOR x and assign it to another
 * NVECTOR z: z = c * x.
 *
 * @param c Scalar value.
 * @param x Nvector to multiply by c.
 * @param z Output Nvector.
 */
static void
snvec_scale(sunrealtype c, N_Vector x, N_Vector z)
{
  // NOTE: we allow vectors not stepped to scale because I think this method is
  // used for copying (e.g. ycur -> yn), which is needed for some df/dt diagnostics
  // to come out right.
  struct gkyl_array *xarr = NV_CONTENT_GKZ(x)->arr;
  struct gkyl_array *zarr = NV_CONTENT_GKZ(z)->arr;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(x)->local_range;

  gkyl_array_set_range(zarr, c, xarr, local_range);
}

/**
 * Component-wise division of two SUNDIALS NVECTORS.
 *   w = u / v.
 *
 * @param u Numerator Nvector.
 * @param v Denominator Nvector.
 * @param w Output Nvector.
 */
static void
snvec_div(N_Vector u, N_Vector v, N_Vector w)
{
  bool not_stepped = NV_CONTENT_GKZ(w)->not_stepped;
  if (not_stepped)
    return;

  struct gkyl_array *u_arr = NV_CONTENT_GKZ(u)->arr;
  struct gkyl_array *v_arr = NV_CONTENT_GKZ(v)->arr;
  struct gkyl_array *w_arr = NV_CONTENT_GKZ(w)->arr;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(u)->local_range;

  gkyl_array_comp_op_range(w_arr, GKYL_DIV, SUN_RCONST(1.0), u_arr,
    SUN_RCONST(0.0), v_arr, local_range);
}

/**
 * Component-wise absolute value.
 *   v_k = |u_k|.
 *
 * @param u Input Nvector.
 * @param v Output Nvector.
 */
static void
snvec_abs(N_Vector u, N_Vector v)
{
  bool not_stepped = NV_CONTENT_GKZ(v)->not_stepped;
  if (not_stepped)
    return;

  struct gkyl_array *u_arr = NV_CONTENT_GKZ(u)->arr;
  struct gkyl_array *v_arr = NV_CONTENT_GKZ(v)->arr;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(u)->local_range;

  gkyl_array_comp_op_range(v_arr, GKYL_ABS, SUN_RCONST(1.0), u_arr,
    SUN_RCONST(1.0), v_arr, local_range);
}

/**
 * Component-wise reciprocal:
 *   v_k = 1/u_k.
 *
 * @param u Input Nvector.
 * @param v Output Nvector.
 */
static void
snvec_inv(N_Vector u, N_Vector v)
{
  bool not_stepped = NV_CONTENT_GKZ(v)->not_stepped;
  if (not_stepped)
    return;

  struct gkyl_array *u_arr       = NV_CONTENT_GKZ(u)->arr;
  struct gkyl_array *v_arr       = NV_CONTENT_GKZ(v)->arr;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(u)->local_range;

  gkyl_array_comp_op_range(v_arr, GKYL_INV, SUN_RCONST(1.0), u_arr,
    SUN_RCONST(1.0), v_arr, local_range);
}

/**
 * Add a constant component-wise:
 *   v_k = u_k + x.
 *
 * @param u Input Nvector.
 * @param x Scala value to add.
 * @param v Output Nvector.
 */
static void
snvec_add_const(N_Vector u, sunrealtype x, N_Vector v)
{
  bool not_stepped = NV_CONTENT_GKZ(v)->not_stepped;
  if (not_stepped)
    return;

  struct gkyl_array *u_arr       = NV_CONTENT_GKZ(u)->arr;
  struct gkyl_array *v_arr       = NV_CONTENT_GKZ(v)->arr;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(u)->local_range;

  gkyl_array_copy_range(v_arr, u_arr, local_range);

  sunindextype ncomp = u_arr->ncomp;

  for (sunindextype i = 0; i < ncomp; ++i)
    gkyl_array_shiftc_range(v_arr, x, i, local_range);
}

/**
 * Compute the component-wise dot product between two SUNDIALS NVECTORS.
 *
 * @param x Input Nvector.
 * @param y Input Nvector.
 * @return Dot product of x and y.
 */
static sunrealtype
snvec_dot_product(N_Vector x, N_Vector y)
{
  bool no_error_norm = NV_CONTENT_GKZ(x)->no_error_norm;
  if (no_error_norm)
    return 0.0;

  struct gkyl_array *x_arr = NV_CONTENT_GKZ(x)->arr;
  struct gkyl_array *y_arr = NV_CONTENT_GKZ(y)->arr;
  struct gkyl_comm *comm = NV_CONTENT_GKZ(y)->comm;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(x)->local_range;
  struct gkyl_sundials_reduce_mem *red_mem = NV_CONTENT_GKZ(x)->red_mem;
  bool use_gpu = NV_CONTENT_GKZ(y)->use_gpu;

  // Sum reduce (component-wise) x_i^{(k)} * y_i^{(k)}.
  int ncomp = x_arr->ncomp;
  gkyl_array_reduce_weighted_range(red_mem->red_local, x_arr, y_arr, GKYL_SUM, local_range);
  gkyl_comm_allreduce(comm, GKYL_DOUBLE, GKYL_SUM, ncomp, red_mem->red_local, red_mem->red_global);

  if (use_gpu)
    gkyl_cu_memcpy(red_mem->red_global_ho, red_mem->red_global, ncomp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(red_mem->red_global_ho, red_mem->red_global, ncomp*sizeof(double));

  sunrealtype dot_prod = 0.0;
  for (sunindextype i = 0; i < ncomp; ++i)
    dot_prod += red_mem->red_global_ho[i];

  return dot_prod;
}

/**
 * Compute the absolute max norm of a SUNDIALS NVECTOR (component-wise).
 *
 * @param u Input Nvector.
 * @return Absolute max norm.
 */
static sunrealtype
snvec_abs_max_norm(N_Vector u)
{
  bool no_error_norm = NV_CONTENT_GKZ(u)->no_error_norm;
  if (no_error_norm)
    return 0.0;

  struct gkyl_array *u_arr = NV_CONTENT_GKZ(u)->arr;
  struct gkyl_comm *comm = NV_CONTENT_GKZ(u)->comm;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(u)->local_range;
  struct gkyl_sundials_reduce_mem *red_mem = NV_CONTENT_GKZ(u)->red_mem;
  bool use_gpu = NV_CONTENT_GKZ(u)->use_gpu;

  int ncomp = u_arr->ncomp;

  gkyl_array_reduce_range(red_mem->red_local, u_arr, GKYL_ABS_MAX, local_range);
  gkyl_comm_allreduce(comm, GKYL_DOUBLE, GKYL_MAX, ncomp, red_mem->red_local, red_mem->red_global);

  if (use_gpu)
    gkyl_cu_memcpy(red_mem->red_global_ho, red_mem->red_global, ncomp * sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(red_mem->red_global_ho, red_mem->red_global, ncomp * sizeof(double));

  sunrealtype u_abs_max = -1.0;
  for (sunindextype i = 0; i < ncomp; ++i)
    u_abs_max = fmax(u_abs_max, red_mem->red_global_ho[i]);

  return u_abs_max;
}

/**
 * Compute the WRMS norm of a Nvector x (typically the error in the solution
 * vector) given the weight w. The norm is
 *   out = sqrt( (1/N) * sum_i^N w_i |x_i|_C^2 ),
 * where
 *   |x_i|_C^2 = (1/N_b) * sum_{k=1}^{N_b} x_i^{(k)},
 *   N = number of cells,
 *   N_b = number of basis functions per cell,
 *   x_i^{(k)} = k-th coefficient in the i-th cell.
 * The weight w is computed externally, e.g.
 *   w_i = 1/[eps_rel*|y_i|_C+eps_abs]
 * where eps_rel and eps_abs are the relative and absolute tolerance,
 * respectively, and y is the solution vector at some stage.
 *
 * @param c Scalar value.
 * @param x Nvector to multiply by c.
 * @param z Output Nvector.
 */
static sunrealtype
snvec_wrms_norm(N_Vector nvx, N_Vector nvwgt)
{
  bool no_error_norm = NV_CONTENT_GKZ(nvx)->no_error_norm;
  if (no_error_norm)
    return 0.0;

  struct gkyl_array *nvx_arr = NV_CONTENT_GKZ(nvx)->arr;
  struct gkyl_array *nvwgt_arr = NV_CONTENT_GKZ(nvwgt)->arr;
  struct gkyl_comm *comm = NV_CONTENT_GKZ(nvwgt)->comm;
  struct gkyl_range *local_range = NV_CONTENT_GKZ(nvx)->local_range;
  struct gkyl_sundials_reduce_mem *red_mem = NV_CONTENT_GKZ(nvx)->red_mem;
  bool use_gpu = NV_CONTENT_GKZ(nvx)->use_gpu;

  int ncomp = nvx_arr->ncomp;

  // Reduce over cells.
  gkyl_array_reduce_weighted_range(red_mem->red_local, nvx_arr, nvwgt_arr, GKYL_RMS, local_range);
  gkyl_comm_allreduce(comm, GKYL_DOUBLE, GKYL_SUM, ncomp, red_mem->red_local, red_mem->red_global);

  if (use_gpu)
    gkyl_cu_memcpy(red_mem->red_global_ho, red_mem->red_global, ncomp * sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(red_mem->red_global_ho, red_mem->red_global, ncomp * sizeof(double));

  // Reduce over components.
  //  sunrealtype red_out = 0.0;
  //  for (sunindextype i = 0; i < ncomp; ++i) red_out += red_ho[i];
  // Use the 0th component because each component should have the same result.
  sunrealtype red_out = red_mem->red_global_ho[0];

  // Divide by the global number of cells: the sum above was reduced over all ranks.
  red_out = SUNRsqrt(red_out / NV_CONTENT_GKZ(nvx)->global_volume);

  return red_out;
}

/**
 * Construct a new SUNDIALS NVECTOR without the underlying Gkeyll array.
 *
 * @param sunctx SUNDIALS context object.
 * @return A new nvector.
 */
static N_Vector
snvec_new_empty(SUNContext sunctx)
{
  N_Vector nvec;
  N_VectorContent_Gkeyll content;

  // Create an empty vector object.
  nvec = 0;
  nvec = N_VNewEmpty(sunctx);
  if (nvec == 0) {
    fprintf(stderr, "Error: Creating empty nvector.\n");
    assert(false);
  }

  // Attach operations.
  // Constructors, destructors, and utility operations.
  nvec->ops->nvclone      = snvec_clone;
  nvec->ops->nvcloneempty = snvec_clone_empty;
  nvec->ops->nvdestroy    = snvec_destroy;

  // Data operations.
  nvec->ops->nvspace     = snvec_space;
  nvec->ops->nvgetlength = snvec_get_length;
  nvec->ops->nvlinearsum = snvec_linear_sum;
  nvec->ops->nvconst     = snvec_const;
  nvec->ops->nvscale     = snvec_scale;
  nvec->ops->nvdiv       = snvec_div;
  nvec->ops->nvabs       = snvec_abs;
  nvec->ops->nvinv       = snvec_inv;
  nvec->ops->nvaddconst  = snvec_add_const;
  nvec->ops->nvdotprod   = snvec_dot_product;
  nvec->ops->nvmaxnorm   = snvec_abs_max_norm;
  nvec->ops->nvwrmsnorm  = snvec_wrms_norm;

  // Create content.
  content = 0;
  content = (N_VectorContent_Gkeyll) gkyl_malloc(sizeof *content);
  if (content == 0) {
    fprintf(stderr, "Error: Allocating content for nvector.\n");
    assert(false);
  }

  // Attach content.
  nvec->content = content;

  // Attach gkyl_array.
  content->own_vector    = SUNFALSE;
  content->use_gpu       = SUNFALSE;
  content->arr           = 0;
  content->comm          = 0;
  content->local_range   = 0;
  content->no_error_norm = SUNFALSE;
  content->not_stepped   = SUNFALSE;
  content->red_mem       = 0;

  return nvec;
}

struct gkyl_array*
snvec_get_array(N_Vector nvin)
{
  return NV_CONTENT_GKZ(nvin)->arr;
}

struct gkyl_array*
smanynvec_get_array(N_Vector manynvin, int nvidx)
{
  N_Vector nvin = N_VGetSubvector_ManyVector(manynvin, nvidx);
  return snvec_get_array(nvin);
}

static void
smanynvec_clear(N_Vector manynvin, double cval)
{
  // Set all DOF of input ManyNVector 'manynvin' to the value 'cval'.
  int num_subvec = N_VGetNumSubvectors_ManyVector(manynvin);
  for (int i=0; i<num_subvec; i++) {
    struct gkyl_array *garr = smanynvec_get_array(manynvin, i);
    gkyl_array_clear(garr, cval);
  }
}

static void 
smanynvec_randomize(N_Vector manynvin, double lo, double up)
{
  // Set all DOF of input ManyNVector 'manynvin' different random numbers.
  int num_subvec = N_VGetNumSubvectors_ManyVector(manynvin);
  for (int i=0; i<num_subvec; i++) {
    struct gkyl_array *garr = smanynvec_get_array(manynvin, i);
    gkyl_array_randomize(garr, lo, up);
  }
}

struct gkyl_sundials *
gkyl_sundials_new(bool use_gpu)
{
  struct gkyl_sundials *gksun = gkyl_malloc(sizeof(*gksun));

  gksun->use_gpu = use_gpu;

  // Create the SUNDIALS context object.
  SUNContext_Create(SUN_COMM_NULL, &gksun->sunctx);

  return gksun;
}

/**
 * Check if rk_method equals:
 *   - opt_A if opt_B=0.
 *   - opt_B if opt_A=0.
 *   - opt_A | opt_B.
 *
 * @param rk_method Input method option(s).
 * @param opt_A First option to check for.
 * @param opt_B Second option to check for.
 * @return Whether input method meets desired condition.
 */
bool
gkyl_sundials_check_rk_method(enum gkyl_sundials_rk_method rk_method,
  enum gkyl_sundials_rk_method opt_A, enum gkyl_sundials_rk_method opt_B)
{
  bool out;
  if (opt_A == 0) {
    out = rk_method & opt_B;
  }
  else if (opt_B == 0) {
    out = rk_method & opt_A;
  }
  else if ((opt_A != 0) && (opt_B != 0)) {
    out = (rk_method & (opt_A | opt_B)) == (opt_A | opt_B);
  }
  else {
    printf("sundials.c gkyl_sundials_check_rk_method: cannot pass 0 for both opt_A and opt_B.\n");
    assert(false);
  }
  return out;
}

/**
 * Translate the enum used to indicate LSRK SSP method in Gkeyll,
 * to the LSRK SSP method types in SUNDIALS.
 *
 * @param gk_lsrk_method LSRK method Gkeyll enum.
 * @return SUNDIALS LSRK method flag.
 */
static ARKODE_LSRKMethodType
gs_translate_gk_to_sundials_method_ssprk(enum gkyl_sundials_rk_method gk_rk_method)
{
  if (gkyl_sundials_check_rk_method(gk_rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_SSP_S_2))
    return ARKODE_LSRK_SSP_S_2; // Optimal 2nd order s-stage SSP RK method.
  else if (gkyl_sundials_check_rk_method(gk_rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_SSP_S_3))
    return ARKODE_LSRK_SSP_S_3; // Optimal 3rd order s-stage SSP RK method.
  else if (gkyl_sundials_check_rk_method(gk_rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_SSP_10_4))
    return ARKODE_LSRK_SSP_10_4; // Optimal 4th order 10-stage SSP RK method.
  else {
    fprintf(stderr, "gs_translate_gk_to_sundials_method_ssprk: wrong input\n");
    assert(false);
  }
  return 0;
}

/**
 * Translate the enum used to indicate LSRK STS method in Gkeyll,
 * to the LSRK STS method types in SUNDIALS.
 *
 * @param gk_lsrk_method LSRK method Gkeyll enum.
 * @return SUNDIALS LSRK method flag.
 */
static ARKODE_LSRKMethodType
gs_translate_gk_to_sundials_method_sts(enum gkyl_sundials_rk_method gk_rk_method)
{
  if (gkyl_sundials_check_rk_method(gk_rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_RKC_2))
    return ARKODE_LSRK_RKC_2; // 2nd order Runge-Kutta-Chebyshev (RKC).
  else if (gkyl_sundials_check_rk_method(gk_rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_RKL_2))
    return ARKODE_LSRK_RKL_2; // 2nd order Runge-Kutta-Legendre (RKL).
  else { 
    fprintf(stderr, "gs_translate_gk_to_sundials_method_sts: wrong input\n");
    assert(false);
  }
  return 0;
}

/**
 * Translate the enum used to indicate operator splitting method in Gkeyll,
 * to the operator splitting method types in SUNDIALS. See
 *   https://sundials.readthedocs.io/en/latest/arkode/Usage/SplittingStep/SplittingStepCoefficients.html#arkode-usage-splittingstep-splittingstepcoefficients
 * for more info.
 *
 * @param gk_opsplit_method LSRK method Gkeyll enum.
 * @return SUNDIALS LSRK method flag.
 */
enum ARKODE_SplittingCoefficientsID
gs_translate_gk_to_sundials_method_opsplit(enum gkyl_sundials_opsplit_method gk_opsplit_method)
{
  if (gk_opsplit_method == GKYL_SUNDIALS_OPSPLIT_METHOD_LIE_TROTTER_1_1_2)
    return ARKODE_SPLITTING_LIE_TROTTER_1_1_2;
  else if (gk_opsplit_method == GKYL_SUNDIALS_OPSPLIT_METHOD_STRANG_2_2_2)
    return ARKODE_SPLITTING_STRANG_2_2_2;
  else if (gk_opsplit_method == GKYL_SUNDIALS_OPSPLIT_METHOD_BEST_2_2_2)
    return ARKODE_SPLITTING_BEST_2_2_2;
  else if (gk_opsplit_method == GKYL_SUNDIALS_OPSPLIT_METHOD_SUZUKI_3_3_2)
    return ARKODE_SPLITTING_SUZUKI_3_3_2;
  else if (gk_opsplit_method == GKYL_SUNDIALS_OPSPLIT_METHOD_RUTH_3_3_2)
    return ARKODE_SPLITTING_RUTH_3_3_2;
  else if (gk_opsplit_method == GKYL_SUNDIALS_OPSPLIT_METHOD_YOSHIDA_4_4_2)
    return ARKODE_SPLITTING_YOSHIDA_4_4_2;
  else if (gk_opsplit_method == GKYL_SUNDIALS_OPSPLIT_METHOD_YOSHIDA_8_6_2)
    return ARKODE_SPLITTING_YOSHIDA_8_6_2;
  else {
    fprintf(stderr, "gs_translate_gk_to_sundials_method_opsplit: wrong input\n");
    assert(false);
  }
  return 0;
}

static void
gkyl_sundials_stepper_init_ssp_rk33(struct gkyl_sundials *gksun,
  struct gkyl_sundials_stepper_inp *inp)
{
  // Initialize Gkeyll's native 3rd order 3-stage SSP RK but through SUNDIALS.
  int flag;
  N_Vector nvin = inp->gsnv->nvec;
  if (nvin == 0) {
    fprintf(stderr, "\ngkyl_sundials_stepper_init_ssp_rk33: Error, gsnv has null N_Vector.\n");
    assert(false);
  }

  // Call ERKStepCreate to initialize the ARK timestepper module and
  // specify the right-hand side function in y'=f(t,y), the initial time
  // T0, and the initial dependent variable vector y.
  gksun->arkode_mem_ssprk = ERKStepCreate(gksun->has_sts? gksun->dfdt_ssprk_func : gksun->dfdt_func,
    inp->t_curr, nvin, nvin->sunctx);
  sundials_check_flag((void*)gksun->arkode_mem_ssprk, "ERKStepCreate", 0);

  // Set routines.
  inp->app_ctx->arkode_mem_ssprk = gksun->arkode_mem_ssprk;
  flag = ARKodeSetUserData(gksun->arkode_mem_ssprk, inp->app_ctx);
  sundials_check_flag(&flag, "ARKodeSetUserData", 1);

  // Set the initial step size.
  sunrealtype dt_stab = 0.1; // Recalculated and passed with _reset method.
  flag = ARKodeSetInitStep(gksun->arkode_mem_ssprk, dt_stab);
  sundials_check_flag(&flag, "ARKodeSetInitStep", 1);

  if (!gksun->has_sts) {
    // Set CFL stable time step size function.
    // Don't do it if we are doing operator split as in that case the time step
    // is set in the update function in gyrokinetic.c.
    flag = ARKodeSetStabilityFn(gksun->arkode_mem_ssprk, gksun->cfl_stable_dt_func, inp->app_ctx);
    sundials_check_flag(&flag, "ARKodeSetStabilityFn", 1);
  }

  // Set CFL safety factor one to ensure exact use of the stable time step size.
  flag = ARKodeSetCFLFraction(gksun->arkode_mem_ssprk, SUN_RCONST(1.0));
  sundials_check_flag(&flag, "ARKodeSetCFLFraction", 1);

  // Specify max number of steps allowed.
  flag = ARKodeSetMaxNumSteps(gksun->arkode_mem_ssprk, inp->max_steps);
  sundials_check_flag(&flag, "ARKodeSetMaxNumSteps", 1);

  // Insert the SSP33 Butcher tableau
  ARKodeButcherTable B_ssp33;
  B_ssp33 = ARKodeButcherTable_Alloc(3, SUNTRUE);
  sundials_check_flag((void*)(B_ssp33), "ARKodeButcherTable_Alloc", 0);
  B_ssp33->A[1][0] = SUN_RCONST(1.0);
  B_ssp33->A[2][0] = SUN_RCONST(0.25);
  B_ssp33->A[2][1] = SUN_RCONST(0.25);
  B_ssp33->b[0]    = SUN_RCONST(1.0/6.0);
  B_ssp33->b[1]    = SUN_RCONST(1.0/6.0);
  B_ssp33->b[2]    = SUN_RCONST(2.0/3.0);
  B_ssp33->c[1]    = SUN_RCONST(1.0);
  B_ssp33->c[2]    = SUN_RCONST(0.5);
  B_ssp33->q       = 3;
  B_ssp33->p       = 3;    // dummy embedding order
  B_ssp33->d[0]    = SUN_RCONST(1.0/6.0); // dummy embedding
  B_ssp33->d[1]    = SUN_RCONST(1.0/6.0);
  B_ssp33->d[2]    = SUN_RCONST(2.0/3.0);

  flag = ERKStepSetTable(gksun->arkode_mem_ssprk, B_ssp33);
  sundials_check_flag(&flag, "ERKStepSetTable", 1);

  // Set pre/post processing methods in arkode mem.
  flag = ARKodeSetPreRhsFn(gksun->arkode_mem_ssprk, gksun->pre_process_rk_stage_ssprk_func);
  sundials_check_flag(&flag, "ARKodeSetPreRhsFn", 1);
  flag = ARKodeSetPostprocessStageFn(gksun->arkode_mem_ssprk, gksun->post_process_rk_stage_ssprk_func);
  sundials_check_flag(&flag, "ARKodeSetPostprocessStageFn", 1);
  flag = ARKodeSetPreStepFn(gksun->arkode_mem_ssprk, gksun->pre_process_step_ssprk_func);
  sundials_check_flag(&flag, "ARKodeSetPreStepFn", 1);
  flag = ARKodeSetPostStepFn(gksun->arkode_mem_ssprk, gksun->post_process_step_ssprk_func);
  sundials_check_flag(&flag, "ARKodeSetPostStepFn", 1);

  // Free the Butcher tableau.
  ARKodeButcherTable_Free(B_ssp33);
}

static void
gkyl_sundials_stepper_init_ssp_rk(struct gkyl_sundials *gksun,
  struct gkyl_sundials_stepper_inp *inp)
{
  // Initialize a SSP RK time stepper (not 3-stage SSP RK3).
  int flag;
  N_Vector nvin = inp->gsnv->nvec;
  if (nvin == 0) {
    fprintf(stderr, "\ngkyl_sundials_stepper_init_ssp_rk: Error, gsnv has null N_Vector.\n");
    assert(false);
  }

  // Call LSRKStepCreateSSP to initialize the ARK timestepper module and
  // specify the right-hand side function in dfdt, the initial time
  // t_curr, and the initial dependent variable vector nvin.
  gksun->arkode_mem_ssprk = LSRKStepCreateSSP(gksun->has_sts? gksun->dfdt_ssprk_func : gksun->dfdt_func,
    inp->t_curr, nvin, nvin->sunctx);
  sundials_check_flag((void*)gksun->arkode_mem_ssprk, "LSRKStepCreateSSP", 0);

  // Set user data (app pointer).
  inp->app_ctx->arkode_mem_ssprk = gksun->arkode_mem_ssprk;
  flag = ARKodeSetUserData(gksun->arkode_mem_ssprk, inp->app_ctx);
  sundials_check_flag(&flag, "ARKodeSetUserData", 1);

  // Specify tolerances.
  flag = ARKodeSStolerances(gksun->arkode_mem_ssprk, inp->rel_tol, inp->abs_tol);
  sundials_check_flag(&flag, "ARKStepSStolerances", 1);

  // Specify max number of steps allowed.
  flag = ARKodeSetMaxNumSteps(gksun->arkode_mem_ssprk, inp->max_steps);
  sundials_check_flag(&flag, "ARKodeSetMaxNumSteps", 1);

  // Specify the SSP method.
  flag = LSRKStepSetSSPMethod(gksun->arkode_mem_ssprk, gs_translate_gk_to_sundials_method_ssprk(inp->rk_method));
  sundials_check_flag(&flag, "LSRKStepSetSSPMethod", 1);

  if ( !(gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_SSP_10_4)) ) {
    // Specify the number of SSP stages.
    flag = LSRKStepSetNumSSPStages(gksun->arkode_mem_ssprk, inp->num_stages);
    sundials_check_flag(&flag, "LSRKStepSetNumSSPStages", 1);
  }

  // Attach the error function.
  flag = ARKodeWFtolerances(gksun->arkode_mem_ssprk, gksun->snvec_efun_cell_norm_func);
  sundials_check_flag(&flag, "ARKodeWFtolerances", 1);

  // Set pre/post processing methods in arkode mem.
  flag = ARKodeSetPreRhsFn(gksun->arkode_mem_ssprk, gksun->pre_process_rk_stage_ssprk_func);
  sundials_check_flag(&flag, "ARKodeSetPreRhsFn", 1);
  flag = ARKodeSetPostprocessStageFn(gksun->arkode_mem_ssprk, gksun->post_process_rk_stage_ssprk_func);
  sundials_check_flag(&flag, "ARKodeSetPostprocessStageFn", 1);
  flag = ARKodeSetPreStepFn(gksun->arkode_mem_ssprk, gksun->pre_process_step_ssprk_func);
  sundials_check_flag(&flag, "ARKodeSetPreStepFn", 1);
  flag = ARKodeSetPostStepFn(gksun->arkode_mem_ssprk, gksun->post_process_step_ssprk_func);
  sundials_check_flag(&flag, "ARKodeSetPostStepFn", 1);
}

static void
gkyl_sundials_sts_set_defaults(struct gkyl_sundials *gksun,
  struct gkyl_sundials_stepper_inp *inp)
{
  // Set default values of STS inputs (shared by LSRKStep STS and ExtSTS) if
  // the user didn't provide a value. The default maximum number of stages is
  // set by each method.
  // A tolerance <= 0 means not set (inputs are zero-initialized, and a zero
  // tolerance would make the error weights infinite where the solution is 0).
  if (inp->rel_tol <= 0.0)
    inp->rel_tol = 1e-7;
  
  if (inp->abs_tol <= 0.0)
    inp->abs_tol = 1e-12;

  if (inp->max_steps == 0)
    inp->max_steps = 1e5;

  if (inp->gsnv == NULL) {
    fprintf(stderr, "\ngkyl_sundials_sts_set_defaults: Error, gsnv is null.\n");
    assert(false);
  }

  if (inp->app_ctx == NULL) {
    fprintf(stderr, "\ngkyl_sundials_sts_set_defaults: Error, app_ctx is null.\n");
    assert(false);
  }

  gksun->dom_eig_est = NULL;

  if (inp->dee_max_iter <= 0)
    inp->dee_max_iter = 1e3;
  
  if (fabs(inp->dee_rel_tol) <= 0.0)
    inp->dee_rel_tol = 0.01;
  
  if (inp->dee_num_init_warmups < 0)
    inp->dee_num_init_warmups = 0;
  
  if (inp->dee_num_succ_warmups < 0)
    inp->dee_num_succ_warmups = 0;
  
  if (inp->dee_frequency < 0)
    inp->dee_frequency = 10;

  if (inp->dee_safety_fac <= 1.0)
    inp->dee_safety_fac = 1.1;
}

static void
gkyl_sundials_sts_dom_eig_init(struct gkyl_sundials *gksun,
  struct gkyl_sundials_stepper_inp *inp, void *arkode_mem_sts, ARKDomEigFn gk_dom_eig_func)
{
  // Set up the dominant eigenvalue estimation (DEE) of an LSRKStep STS
  // integrator, either with Gkeyll's estimate (gk_dom_eig_func) or with
  // SUNDIALS' power iteration estimator.
  int flag;
  N_Vector nvin = inp->gsnv->nvec;

  // Randomize the input vector (needed by DEE).
  smanynvec_randomize(nvin, 0.5, 1.5);

  if (inp->dee_by_gkeyll) {
    // Gkeyll estimates the dominant eigenvalue.
    flag = LSRKStepSetDomEigFn(arkode_mem_sts, gk_dom_eig_func);
    sundials_check_flag(&flag, "LSRKStepSetDomEigFn", 1);
  }
  else {
    // Set the initial eigenvector for DEE. Should reset after ICs are set.
    gksun->dom_eig_est = SUNDomEigEstimator_Power(nvin, inp->dee_max_iter,
      inp->dee_rel_tol, nvin->sunctx);
    sundials_check_flag((void*)gksun->dom_eig_est, "SUNDomEigEstimator_Power", 0);
  
    flag = LSRKStepSetDomEigEstimator(arkode_mem_sts, gksun->dom_eig_est);
    sundials_check_flag(&flag, "LSRKStepSetDomEigEstimator", 1);
  
    flag = LSRKStepSetNumDomEigEstInitPreprocessIters(arkode_mem_sts, inp->dee_num_init_warmups);
    sundials_check_flag(&flag, "LSRKStepSetNumDomEigEstInitPreprocessIters", 1);
  
    flag = LSRKStepSetNumDomEigEstPreprocessIters(arkode_mem_sts, inp->dee_num_succ_warmups);
    sundials_check_flag(&flag, "LSRKStepSetNumDomEigEstPreprocessIters", 1);
  }

  // Zero out input vector (that's what's expected/assumed during initialization).
  smanynvec_clear(nvin, 0.0);

  // Specify after how many successful steps dom_eig is recomputed.
  // Note that nsteps = 0 refers to constant dominant eigenvalue.
  flag = LSRKStepSetDomEigFrequency(arkode_mem_sts, inp->dee_frequency);
  sundials_check_flag(&flag, "LSRKStepSetDomEigFrequency", 1);

  // Specify safety factor for user provided dom_eig.
  flag = LSRKStepSetDomEigSafetyFactor(arkode_mem_sts, inp->dee_safety_fac);
  sundials_check_flag(&flag, "LSRKStepSetDomEigSafetyFactor", 1);
}

static void
gkyl_sundials_stepper_init_sts(struct gkyl_sundials *gksun,
  struct gkyl_sundials_stepper_inp *inp)
{
  // Initialize the Super Time Stepping integrator.

  // Set default values if user didn't provide a value.
  gkyl_sundials_sts_set_defaults(gksun, inp);

  if (inp->max_num_stages <= 0)
    inp->max_num_stages = 200;
  // Finished setting default values.

  int flag;
  N_Vector nvin = inp->gsnv->nvec;
  if (nvin == 0) {
    fprintf(stderr, "\ngkyl_sundials_stepper_init_sts: Error, gsnv has null N_Vector.\n");
    assert(false);
  }

  // Call LSRKStepCreateSTS to initialize the ARK timestepper module and
  // specify the right-hand side function in dfdt, the initial time
  // t_curr, and the initial dependent variable vector nvin.
  gksun->arkode_mem_sts = LSRKStepCreateSTS(gksun->has_ssprk? gksun->dfdt_sts_func : gksun->dfdt_func, inp->t_curr, nvin, nvin->sunctx);
  sundials_check_flag((void*)gksun->arkode_mem_sts, "LSRKStepCreateSTS", 0);

  // Set user data (app pointer).
  inp->app_ctx->arkode_mem_sts = gksun->arkode_mem_sts;
  flag = ARKodeSetUserData(gksun->arkode_mem_sts, inp->app_ctx);
  sundials_check_flag(&flag, "ARKodeSetUserData", 1);

  // Specify tolerances.
  flag = ARKodeSStolerances(gksun->arkode_mem_sts, inp->rel_tol, inp->abs_tol);
  sundials_check_flag(&flag, "ARKStepSStolerances", 1);

  // Dominant eigenvalue estimation.
  gkyl_sundials_sts_dom_eig_init(gksun, inp, gksun->arkode_mem_sts,
    gksun->has_ssprk? gksun->gk_dom_eig_sts_func : gksun->gk_dom_eig_func);

  // Specify max number of stages allowed.
  flag = LSRKStepSetMaxNumStages(gksun->arkode_mem_sts, inp->max_num_stages);
  sundials_check_flag(&flag, "LSRKStepSetMaxNumStages", 1);

  // Specify max number of steps allowed.
  flag = ARKodeSetMaxNumSteps(gksun->arkode_mem_sts, inp->max_steps);
  sundials_check_flag(&flag, "ARKodeSetMaxNumSteps", 1);

  // Specify the STS method.
  flag = LSRKStepSetSTSMethod(gksun->arkode_mem_sts, gs_translate_gk_to_sundials_method_sts(inp->rk_method));
  sundials_check_flag(&flag, "LSRKStepSetSTSMethod", 1);

  // Attach the error function.
  flag = ARKodeWFtolerances(gksun->arkode_mem_sts, gksun->snvec_efun_cell_norm_func);
  sundials_check_flag(&flag, "ARKodeWFtolerances (sts)", 1);

  // Set pre/post processing methods in arkode mem.
  flag = ARKodeSetPreRhsFn(gksun->arkode_mem_sts, gksun->pre_process_rk_stage_sts_func);
  sundials_check_flag(&flag, "ARKodeSetPreRhsFn", 1);
  flag = ARKodeSetPostprocessStageFn(gksun->arkode_mem_sts, gksun->post_process_rk_stage_sts_func);
  sundials_check_flag(&flag, "ARKodeSetPostprocessStageFn", 1);
  flag = ARKodeSetPreStepFn(gksun->arkode_mem_sts, gksun->pre_process_step_sts_func);
  sundials_check_flag(&flag, "ARKodeSetPreStepFn", 1);
  flag = ARKodeSetPostStepFn(gksun->arkode_mem_sts, gksun->post_process_step_sts_func);
  sundials_check_flag(&flag, "ARKodeSetPostStepFn", 1);
}

static void
gkyl_sundials_stepper_init_mri_extsts(struct gkyl_sundials *gksun,
  struct gkyl_sundials_stepper_inp *inp)
{
  // Initialize an extended super time stepping (ExtSTS) method: an MRI method
  // whose outer (explicit) stages evolve the terms otherwise stepped with
  // SSP-RK, and whose inner STS integrator (one step per outer stage) evolves
  // the terms otherwise stepped with STS.

  // Set default values if user didn't provide a value. If max_num_stages is
  // not given we keep the ExtSTS default (10000), because a single STS step
  // spans a whole outer stage.
  gkyl_sundials_sts_set_defaults(gksun, inp);
  // Finished setting default values.

  int flag;
  N_Vector nvin = inp->gsnv->nvec;
  if (nvin == 0) {
    fprintf(stderr, "\ngkyl_sundials_stepper_init_mri_extsts: Error, gsnv has null N_Vector.\n");
    assert(false);
  }

  // Create the MRIStep ExtSTS integrator, with the STS terms as the diffusive
  // function fd, the SSP-RK terms as the explicit function fe, and no
  // implicit function fi.
  gksun->arkode_mem_outer = MRIStepCreateExtSTS(gksun->dfdt_sts_func, gksun->dfdt_ssprk_func, NULL,
    inp->t_curr, nvin, nvin->sunctx);
  sundials_check_flag((void*)gksun->arkode_mem_outer, "MRIStepCreateExtSTS", 0);

  // Access the inner LSRKStep STS integrator (owned by the MRIStep integrator).
  flag = MRIStepGetSTS(gksun->arkode_mem_outer, &gksun->arkode_mem_sts);
  sundials_check_flag(&flag, "MRIStepGetSTS", 1);
  sundials_check_flag(gksun->arkode_mem_sts, "MRIStepGetSTS", 0);

  // The MRIStep integrator also plays the role of the SSP-RK stepper: the
  // SSP-RK df/dt and the SSP-RK step hooks query the stage index, current
  // state and step size of arkode_mem_ssprk.
  gksun->arkode_mem_ssprk = gksun->arkode_mem_outer;
  inp->app_ctx->arkode_mem_ssprk = gksun->arkode_mem_ssprk;
  inp->app_ctx->arkode_mem_sts = gksun->arkode_mem_sts;

  // Set user data (app pointer). This is passed on to the inner STS integrator.
  flag = ARKodeSetUserData(gksun->arkode_mem_outer, inp->app_ctx);
  sundials_check_flag(&flag, "ARKodeSetUserData", 1);

  // Specify tolerances and attach the error function (used with adaptive steps).
  flag = ARKodeSStolerances(gksun->arkode_mem_outer, inp->rel_tol, inp->abs_tol);
  sundials_check_flag(&flag, "ARKodeSStolerances", 1);
  flag = ARKodeWFtolerances(gksun->arkode_mem_outer, gksun->snvec_efun_cell_norm_func);
  sundials_check_flag(&flag, "ARKodeWFtolerances (extsts)", 1);

  // Specify the STS method of the inner integrator.
  flag = LSRKStepSetSTSMethod(gksun->arkode_mem_sts, gs_translate_gk_to_sundials_method_sts(inp->rk_method));
  sundials_check_flag(&flag, "LSRKStepSetSTSMethod", 1);

  // Dominant eigenvalue estimation for the STS terms.
  gkyl_sundials_sts_dom_eig_init(gksun, inp, gksun->arkode_mem_sts, gksun->gk_dom_eig_sts_func);

  if (inp->max_num_stages > 0) {
    // Specify max number of STS stages allowed.
    flag = LSRKStepSetMaxNumStages(gksun->arkode_mem_sts, inp->max_num_stages);
    sundials_check_flag(&flag, "LSRKStepSetMaxNumStages", 1);
  }

  if (inp->mri_coupling != 0) {
    // Set the outer MRI method.
    MRIStepCoupling mri_coupling = MRIStepCoupling_LoadTable((ARKODE_MRITableID) inp->mri_coupling);
    sundials_check_flag((void*)mri_coupling, "MRIStepCoupling_LoadTable", 0);
    flag = MRIStepSetCoupling(gksun->arkode_mem_outer, mri_coupling);
    sundials_check_flag(&flag, "MRIStepSetCoupling", 1);
    MRIStepCoupling_Free(mri_coupling);
  }

  // Specify max number of steps allowed.
  flag = ARKodeSetMaxNumSteps(gksun->arkode_mem_outer, inp->max_steps);
  sundials_check_flag(&flag, "ARKodeSetMaxNumSteps", 1);

  if (!inp->fixed_step) {
    // Adaptive outer steps (temporal error control), limited by the CFL stable
    // time step of the SSP-RK terms. A fixed step is instead set from the same
    // CFL limit in the reset and update functions.
    flag = ARKodeSetStabilityFn(gksun->arkode_mem_outer, gksun->cfl_stable_dt_ssprk_func, inp->app_ctx);
    sundials_check_flag(&flag, "ARKodeSetStabilityFn", 1);

    // Set CFL safety factor one to ensure exact use of the stable time step size.
    flag = ARKodeSetCFLFraction(gksun->arkode_mem_outer, SUN_RCONST(1.0));
    sundials_check_flag(&flag, "ARKodeSetCFLFraction", 1);

    // Time step adaptivity controller and parameters. Compared with the default
    // controller, these reduce rejected steps by ~4-5x when the step is limited
    // by the stability of the outer method rather than by the CFL cap, with
    // similar accepted step sizes (rt_gk_sheath_2x2v_p1).
    gksun->adapt_controller = SUNAdaptController_Soderlind(nvin->sunctx);
    sundials_check_flag((void*)gksun->adapt_controller, "SUNAdaptController_Soderlind", 0);
    flag = SUNAdaptController_SetParams_Soderlind(gksun->adapt_controller,
      SUN_RCONST(0.875), SUN_RCONST(-0.5), SUN_RCONST(0.0), SUN_RCONST(0.0), SUN_RCONST(0.0));
    sundials_check_flag(&flag, "SUNAdaptController_SetParams_Soderlind", 1);
    flag = ARKodeSetAdaptController(gksun->arkode_mem_outer, gksun->adapt_controller);
    sundials_check_flag(&flag, "ARKodeSetAdaptController", 1);
    flag = ARKodeSetSafetyFactor(gksun->arkode_mem_outer, SUN_RCONST(0.9));
    sundials_check_flag(&flag, "ARKodeSetSafetyFactor", 1);
    flag = ARKodeSetMaxGrowth(gksun->arkode_mem_outer, SUN_RCONST(1.5));
    sundials_check_flag(&flag, "ARKodeSetMaxGrowth", 1);
    flag = ARKodeSetMaxFirstGrowth(gksun->arkode_mem_outer, SUN_RCONST(1.5));
    sundials_check_flag(&flag, "ARKodeSetMaxFirstGrowth", 1);
    flag = ARKodeSetMinReduction(gksun->arkode_mem_outer, SUN_RCONST(0.25));
    sundials_check_flag(&flag, "ARKodeSetMinReduction", 1);
    flag = ARKodeSetMaxEFailGrowth(gksun->arkode_mem_outer, SUN_RCONST(0.2));
    sundials_check_flag(&flag, "ARKodeSetMaxEFailGrowth", 1);
  }

  // Set pre/post processing methods in arkode mem. The SSP-RK step hooks are
  // used at the beginning/end of each outer step, and the outer stages are
  // post-processed so that fields and boundary conditions are consistent with
  // each stage before the SSP-RK terms are evaluated on it.
  flag = ARKodeSetPreRhsFn(gksun->arkode_mem_outer, gksun->pre_process_rk_stage_extsts_func);
  sundials_check_flag(&flag, "ARKodeSetPreRhsFn", 1);
  flag = ARKodeSetPostprocessStageFn(gksun->arkode_mem_outer, gksun->post_process_stage_extsts_func);
  sundials_check_flag(&flag, "ARKodeSetPostprocessStageFn", 1);
  flag = ARKodeSetPreStepFn(gksun->arkode_mem_outer, gksun->pre_process_step_ssprk_func);
  sundials_check_flag(&flag, "ARKodeSetPreStepFn", 1);
  flag = ARKodeSetPostStepFn(gksun->arkode_mem_outer, gksun->post_process_step_ssprk_func);
  sundials_check_flag(&flag, "ARKodeSetPostStepFn", 1);
}

void
gkyl_sundials_stepper_init(struct gkyl_sundials *gksun,
  struct gkyl_sundials_stepper_inp *inp)
{
  gksun->is_extsts = false;
  gksun->adapt_controller = NULL;
  if (gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_MRI_EXTSTS)) {
    // Extended super time stepping (ExtSTS).
    if ( gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_RK_SSP_3_3) ||
         gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_SSP_S_2) ||
         gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_SSP_S_3) ||
         gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_SSP_10_4) ) {
      fprintf(stderr, "\ngkyl_sundials_stepper_init: Error, GKYL_SUNDIALS_METHOD_MRI_EXTSTS cannot be combined with an SSP-RK method.\n");
      assert(false);
    }
    if ( !(gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_RKC_2) ||
           gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_RKL_2)) ) {
      printf("gkyl_sundials_stepper_init: no STS method given with GKYL_SUNDIALS_METHOD_MRI_EXTSTS, using GKYL_SUNDIALS_METHOD_LSRK_RKL_2.\n");
      inp->rk_method |= GKYL_SUNDIALS_METHOD_LSRK_RKL_2;
    }

    // Maximum number steps.
    if (inp->max_steps == 0)
      inp->max_steps = 100000;

    gksun->stepper_inp = inp; // Store stepper inputs.
    gksun->app_ctx = inp->app_ctx; // Copy pointer to app pointer.

    gksun->has_ssprk = false;
    gksun->has_sts = true;
    gksun->is_extsts = true;
    gksun->is_opsplit = false;

    gkyl_sundials_stepper_init_mri_extsts(gksun, inp);

    // Store outermost stepper mem in context so it can be accessed in methods.
    inp->app_ctx->arkode_mem_outer = gksun->arkode_mem_outer;
    return;
  }

  // Set default values if user didn't provide a value.
  bool is_lsrk_ssp_s_3 = gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_SSP_S_3);
  if ( is_lsrk_ssp_s_3 && (inp->num_stages == 3) ) {
    // Use 3rd order 3-stage SSP RK, without embedding, adapting dt using Gkeyll's CFL constraint.
    if (gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_RKC_2))
      inp->rk_method = GKYL_SUNDIALS_METHOD_RK_SSP_3_3 | GKYL_SUNDIALS_METHOD_LSRK_RKC_2;
    else if (gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_RKL_2))
      inp->rk_method = GKYL_SUNDIALS_METHOD_RK_SSP_3_3 | GKYL_SUNDIALS_METHOD_LSRK_RKL_2;
    else
      inp->rk_method = GKYL_SUNDIALS_METHOD_RK_SSP_3_3;
  }
  
  // Maximum number steps (not used if using ARK_ONE_STEP, unless using operator-split).
  if (inp->max_steps == 0)
    inp->max_steps = 100000;
  // Finished setting default values.

  gksun->stepper_inp = inp; // Store stepper inputs.
  gksun->app_ctx = inp->app_ctx; // Copy pointer to app pointer.

  gksun->has_ssprk = gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_RK_SSP_3_3) ||
                     gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_SSP_S_2) ||
                     gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_SSP_S_3) ||
                     gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_SSP_10_4);
  gksun->has_sts = gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_RKC_2) ||
                   gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_RKL_2);

  if (inp->fixed_step && !gksun->has_ssprk) {
    // Fixed steps are set by the CFL limit of the SSP-RK terms.
    fprintf(stderr, "\ngkyl_sundials_stepper_init: Error, fixed_step requires a method with SSP-RK terms.\n");
    assert(false);
  }

  if (gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_RK_SSP_3_3)) {
    // Gkeyll's native 3rd order 3-stage SSP RK.
    gkyl_sundials_stepper_init_ssp_rk33(gksun, inp);
  }
  else if ( gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_SSP_S_2) ||
            gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_SSP_S_3) ||
            gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_SSP_10_4) ) {
    // Other SSP RK methods (with embedding to adapt dt).
    gkyl_sundials_stepper_init_ssp_rk(gksun, inp);
  }

  if ( gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_RKC_2) ||
       gkyl_sundials_check_rk_method(inp->rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_RKL_2) ) {
    // Super time stepping.
    gkyl_sundials_stepper_init_sts(gksun, inp);
  }

  gksun->is_opsplit = false;
  if (gksun->has_ssprk && gksun->has_sts) {
    // Create operator splitting integrator.
    gksun->is_opsplit = true;
    int flag;
    flag = ARKodeCreateSUNStepper(gksun->arkode_mem_ssprk, &gksun->stepper_ssprk);
    sundials_check_flag(&flag, "ARKodeCreateSUNStepper", 1);

    flag = ARKodeCreateSUNStepper(gksun->arkode_mem_sts, &gksun->stepper_sts);
    sundials_check_flag(&flag, "ARKodeCreateSUNStepper", 1);

    SUNStepper steppers[] = {gksun->stepper_ssprk, gksun->stepper_sts};
    int num_partitions = sizeof(steppers)/sizeof(steppers[0]);

    N_Vector nvin = inp->gsnv->nvec;
    if (nvin == 0) {
      fprintf(stderr, "\ngkyl_sundials_stepper_init: Error, gsnv has null N_Vector.\n");
      assert(false);
    }

    gksun->arkode_mem_outer = SplittingStepCreate(steppers, num_partitions, inp->t_curr, nvin, nvin->sunctx);
    sundials_check_flag(gksun->arkode_mem_outer, "SplittingStepCreate", 0);

    flag = ARKodeSetUserData(gksun->arkode_mem_outer, inp->app_ctx);

    if (inp->opsplit_method != GKYL_SUNDIALS_OPSPLIT_METHOD_NONE) {
      // Set operator split method.
      SplittingStepCoefficients opsplit_coeffs = SplittingStepCoefficients_LoadCoefficients(gs_translate_gk_to_sundials_method_opsplit(inp->opsplit_method));
      sundials_check_flag(opsplit_coeffs, "SplittingStepCoefficients_LoadCoefficients", 0);
      flag = SplittingStepSetCoefficients(gksun->arkode_mem_outer, opsplit_coeffs);
      sundials_check_flag(&flag, "SplittingStepSetCoefficients", 1);
      SplittingStepCoefficients_Destroy(&opsplit_coeffs);
    }
  
    // Set pre/post processing methods in arkode mem.
    flag = ARKodeSetPreStepFn(gksun->arkode_mem_outer, gksun->pre_process_step_opsplit_func);
    sundials_check_flag(&flag, "ARKodeSetPreStepFn", 1);
    flag = ARKodeSetPostStepFn(gksun->arkode_mem_outer, gksun->post_process_step_opsplit_func);
    sundials_check_flag(&flag, "ARKodeSetPostStepFn", 1);
  }
  else {
    gksun->arkode_mem_outer = gksun->has_ssprk? gksun->arkode_mem_ssprk : gksun->arkode_mem_sts;
  }

  // Store outermost stepper mem in context so it can be accessed in methods.
  inp->app_ctx->arkode_mem_outer = gksun->arkode_mem_outer;

}

void
gkyl_sundials_arkode_reset(struct gkyl_sundials *gksun, double time,
  struct gkyl_sundials_nvec *gsmanynv, struct gkyl_sundials_nvec *gsmanynv_buff)
{
  N_Vector manynvin = gsmanynv->nvec;

  int flag;

  flag = ARKodeReset(gksun->arkode_mem_outer, time, manynvin);
  sundials_check_flag(&flag, "ARKodeReset", 1);

  if (gksun->has_sts) {
    if ( !(gksun->stepper_inp->dee_by_gkeyll) ) {
      // Pass ICs to the eigenvalue estimate.
      flag = SUNDomEigEstimator_SetInitialGuess(gksun->dom_eig_est, manynvin);
      sundials_check_flag(&flag, "SUNDomEigEstimator_SetInitialGuess", 1);
    }
  }

  if (gksun->has_ssprk || gksun->is_extsts) {
    if (gksun->is_opsplit || gksun->is_extsts || gksun->stepper_inp->fixed_step ||
        gkyl_sundials_check_rk_method(gksun->stepper_inp->rk_method, 0, GKYL_SUNDIALS_METHOD_RK_SSP_3_3)) {
      // Estimate the CFL stable dt.
      // Use the temporary NVector.
      N_Vector manynvbuff = gsmanynv_buff->nvec;
  
      // Compute dt.
      if (gksun->is_opsplit || gksun->is_extsts) {
        // For ExtSTS arkode_mem_ssprk is the MRIStep (outer) memory.
        double dt_ssprk;
        flag = gksun->dfdt_ssprk_func(time, manynvin, manynvbuff, gksun->app_ctx);
        flag = gksun->cfl_stable_dt_ssprk_func(manynvin, time, &dt_ssprk, gksun->app_ctx);

        // Set the initial step size of the SSP-RK stepper.
        flag = ARKodeSetInitStep(gksun->arkode_mem_ssprk, dt_ssprk);
        sundials_check_flag(&flag, "ARKodeSetInitStep", 1);

        if (gksun->is_opsplit || gksun->stepper_inp->fixed_step) {
          // Set outer dt to the SSP-RK dt.
          gkyl_sundials_set_fixed_step(gksun, dt_ssprk);
        }
      }
      else {
        double dt_init;
        flag = gksun->dfdt_func(time, manynvin, manynvbuff, gksun->app_ctx);
        flag = gksun->cfl_stable_dt_func(manynvin, time, &dt_init, gksun->app_ctx);
        // Set the initial step size.
        flag = ARKodeSetInitStep(gksun->arkode_mem_ssprk, dt_init);
        sundials_check_flag(&flag, "ARKodeSetInitStep", 1);

        if (gksun->stepper_inp->fixed_step) {
          // Set the (outer) dt to the CFL stable dt.
          gkyl_sundials_set_fixed_step(gksun, dt_init);
        }
      }
    }
  }

}

void
gkyl_sundials_set_fixed_step(struct gkyl_sundials *gksun, double dt)
{
  int flag = ARKodeSetFixedStep(gksun->arkode_mem_outer, dt);
  sundials_check_flag(&flag, "ARKodeSetFixedStep", 1);
}

void
gkyl_sundials_set_fixed_step_ssprk(struct gkyl_sundials *gksun, double dt)
{
  int flag = ARKodeSetFixedStep(gksun->arkode_mem_ssprk, dt);
  sundials_check_flag(&flag, "ARKodeSetFixedStep", 1);
}

void
gkyl_sundials_set_fixed_step_sts(struct gkyl_sundials *gksun, double dt)
{
  int flag = ARKodeSetFixedStep(gksun->arkode_mem_sts, dt);
  sundials_check_flag(&flag, "ARKodeSetFixedStep", 1);
}

double
gkyl_sundials_get_cfl_dt(struct gkyl_sundials *gksun)
{
  double dt;
  int flag = gksun->cfl_stable_dt_func(0, 0, &dt, gksun->app_ctx);
  return dt;
}

double
gkyl_sundials_get_cfl_dt_ssprk(struct gkyl_sundials *gksun)
{
  double dt;
  int flag = gksun->cfl_stable_dt_ssprk_func(0, 0, &dt, gksun->app_ctx);
  return dt;
}

double
gkyl_sundials_get_cfl_dt_sts(struct gkyl_sundials *gksun)
{
  double dt;
  int flag = gksun->cfl_stable_dt_sts_func(0, 0, &dt, gksun->app_ctx);
  return dt;
}

int
gkyl_sundials_evolve(struct gkyl_sundials *gksun, double t_new,
  struct gkyl_sundials_nvec *gsnv, double *t_curr)
{
  N_Vector nvin = gsnv->nvec;

  // Call integrator to evolve the solution to time t_new.
  int flag = ARKodeEvolve(gksun->arkode_mem_outer, t_new, nvin, t_curr, ARK_ONE_STEP);
  sundials_check_flag(&flag, "ARKodeEvolve", 1);
  if (flag != ARK_SUCCESS) {
    fprintf(stderr, "\ngkyl_sundials_evolve: Error, ARKodeEvolve failed with flag = %d.\n", flag);
    return flag;
  }

  if (gksun->stepper_inp->print_arkode_stats) {
    printf("\n\n\nStep successful: advanced to time %g\n", *t_curr);
    flag = ARKodePrintAllStats(gksun->arkode_mem_outer, stdout, SUN_OUTPUTFORMAT_TABLE);
    sundials_check_flag(&flag, "ARKodePrintAllStats", 1);
    if (flag != ARK_SUCCESS) {
      fprintf(stderr, "\ngkyl_sundials_evolve: Error, ARKodePrintAllStats failed with flag = %d.\n", flag);
      return flag;
    }
    if (gksun->dom_eig_est != NULL) {
      flag = SUNDomEigEstimator_Write(gksun->dom_eig_est, stdout);
      sundials_check_flag(&flag, "SUNDomEigEstimator_Write", 1);
      if (flag != ARK_SUCCESS) {
        fprintf(stderr, "\ngkyl_sundials_evolve: Error, SUNDomEigEstimator_Write failed with flag = %d.\n", flag);
        return flag;
      }
    }
  }

  return flag;
}

long
gkyl_sundials_get_num_error_test_failures(struct gkyl_sundials *gksun)
{
  long num_fail = 0;
  if (gksun->has_ssprk) {
    long num_fail_ssprk;
    int flag = ARKodeGetNumErrTestFails(gksun->arkode_mem_outer, &num_fail_ssprk); 
    sundials_check_flag(&flag, "ARKodeGetNumErrTestFails", 1);
    num_fail += num_fail_ssprk;
  }

  if (gksun->has_sts) {
    long num_fail_sts;
    int flag = ARKodeGetNumErrTestFails(gksun->arkode_mem_outer, &num_fail_sts); 
    sundials_check_flag(&flag, "ARKodeGetNumErrTestFails", 1);
    num_fail += num_fail_sts;
  }

  return num_fail;
}

long
gkyl_sundials_get_num_rhs_evals(struct gkyl_sundials *gksun)
{
  long num_evals = 0;
  if (gksun->has_ssprk || gksun->is_extsts) {
    // For ExtSTS arkode_mem_ssprk is the MRIStep memory, whose explicit (fe)
    // evaluations are the SSP-RK terms.
    long num_evals_ssprk = 0;
    int flag = ARKodeGetNumRhsEvals(gksun->arkode_mem_ssprk, 0, &num_evals_ssprk);
    sundials_check_flag(&flag, "ARKodeGetNumRhsEvals", 1);
    num_evals += num_evals_ssprk;
  }

  if (gksun->has_sts) {
    long num_evals_sts = 0;
    int flag = ARKodeGetNumRhsEvals(gksun->arkode_mem_sts, 0, &num_evals_sts); 
    sundials_check_flag(&flag, "ARKodeGetNumRhsEvals", 1);
    num_evals += num_evals_sts;
  }

  return num_evals;
}

double
gkyl_sundials_get_last_dt(struct gkyl_sundials *gksun)
{
  double dt_out;
  int flag = ARKodeGetLastStep(gksun->arkode_mem_outer, &dt_out);
  sundials_check_flag(&flag, "ARKodeGetLastStep", 1);
  return dt_out;
}

double
gkyl_sundials_get_last_dt_ssprk(struct gkyl_sundials *gksun)
{
  double dt_out;
  int flag = ARKodeGetLastStep(gksun->arkode_mem_ssprk, &dt_out);
  sundials_check_flag(&flag, "ARKodeGetLastStep", 1);
  return dt_out;
}

double
gkyl_sundials_get_last_dt_sts(struct gkyl_sundials *gksun)
{
  double dt_out;
  int flag = ARKodeGetLastStep(gksun->arkode_mem_sts, &dt_out);
  sundials_check_flag(&flag, "ARKodeGetLastStep", 1);
  return dt_out;
}

double
gkyl_sundials_get_current_dt(struct gkyl_sundials *gksun)
{
  double dt_out;
  int flag = ARKodeGetCurrentStep(gksun->arkode_mem_outer, &dt_out);
  sundials_check_flag(&flag, "ARKodeGetCurrentStep", 1);
  return dt_out;
}

double
gkyl_sundials_get_current_dt_ssprk(struct gkyl_sundials *gksun)
{
  double dt_out;
  int flag = ARKodeGetCurrentStep(gksun->arkode_mem_ssprk, &dt_out);
  sundials_check_flag(&flag, "ARKodeGetCurrentStep", 1);
  return dt_out;
}

double
gkyl_sundials_get_current_dt_sts(struct gkyl_sundials *gksun)
{
  double dt_out;
  int flag = ARKodeGetCurrentStep(gksun->arkode_mem_sts, &dt_out);
  sundials_check_flag(&flag, "ARKodeGetCurrentStep", 1);
  return dt_out;
}

bool
gkyl_sundials_use_operator_split(struct gkyl_sundials *gksun)
{
  return gksun->is_opsplit;
}

bool
gkyl_sundials_operator_split_in_method(enum gkyl_sundials_rk_method rk_method)
{
  bool has_ssprk = gkyl_sundials_check_rk_method(rk_method, 0, GKYL_SUNDIALS_METHOD_RK_SSP_3_3) ||
                   gkyl_sundials_check_rk_method(rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_SSP_S_2) ||
                   gkyl_sundials_check_rk_method(rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_SSP_S_3) ||
                   gkyl_sundials_check_rk_method(rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_SSP_10_4);
  bool has_sts = gkyl_sundials_check_rk_method(rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_RKC_2) ||
                 gkyl_sundials_check_rk_method(rk_method, 0, GKYL_SUNDIALS_METHOD_LSRK_RKL_2);

  bool is_opsplit = false;
  if (has_ssprk && has_sts)
    is_opsplit = true;

  return is_opsplit;
}

void
gkyl_sundials_release(struct gkyl_sundials *gksun)
{
  if (gksun->is_extsts) {
    // Freeing the MRIStep memory also frees its inner STS memory.
    ARKodeFree(&gksun->arkode_mem_outer);
    // The time step controller is not owned by ARKODE.
    if (gksun->adapt_controller != NULL)
      SUNAdaptController_Destroy(gksun->adapt_controller);
  }
  else {
    if (gksun->has_ssprk)
      ARKodeFree(&gksun->arkode_mem_ssprk);

    if (gksun->has_sts)
      ARKodeFree(&gksun->arkode_mem_sts);

    if (gksun->is_opsplit) {
      SUNStepper_Destroy(&gksun->stepper_ssprk);
      SUNStepper_Destroy(&gksun->stepper_sts);
      ARKodeFree(&gksun->arkode_mem_outer);
    }
  }

  SUNContext_Free(&gksun->sunctx);

  gkyl_free(gksun);
}

struct gkyl_sundials_nvec*
gkyl_sundials_nvec_new(struct gkyl_sundials *gksun, struct gkyl_array *arr,
  struct gkyl_comm *comm, struct gkyl_range *local_range, bool no_error_norm, bool not_stepped)
{
  struct gkyl_sundials_nvec *gsnv = gkyl_malloc(sizeof(*gsnv));

  // Allocate memory needed for reductions
  bool use_gpu = gkyl_array_is_cu_dev(arr);
  int ncomp = arr->ncomp;
  if (use_gpu) {
    gsnv->red_mem.red_local  = gkyl_cu_malloc(ncomp * sizeof(double));
    gsnv->red_mem.red_global = gkyl_cu_malloc(ncomp * sizeof(double));
  }
  else {
    gsnv->red_mem.red_local  = gkyl_malloc(ncomp * sizeof(double));
    gsnv->red_mem.red_global = gkyl_malloc(ncomp * sizeof(double));
  }
  gsnv->red_mem.red_global_ho = gkyl_malloc(ncomp * sizeof(double));

  // Allocate the Sundials NVector.
  N_Vector nvout;
  nvout = 0;
  nvout = snvec_new_empty(gksun->sunctx);
  if (nvout == 0) {
    fprintf(stderr, "Error: Creating new empty nvector.\n");
    assert(false);
  }

  // If a quantity is not stepped, it shouldn't be included in error calculation.
  if (not_stepped)
    assert(no_error_norm);

  NV_CONTENT_GKZ(nvout)->own_vector = SUNFALSE;
  NV_CONTENT_GKZ(nvout)->use_gpu = gkyl_array_is_cu_dev(arr);
  NV_CONTENT_GKZ(nvout)->comm = comm;
  NV_CONTENT_GKZ(nvout)->local_range = local_range;
  // Global number of cells, needed to normalize the WRMS norm consistently across ranks.
  double local_volume = local_range->volume, global_volume = local_volume;
  if (comm) {
    gkyl_comm_allreduce_host(comm, GKYL_DOUBLE, GKYL_SUM, 1, &local_volume, &global_volume);
  }
  NV_CONTENT_GKZ(nvout)->global_volume = global_volume;
  NV_CONTENT_GKZ(nvout)->no_error_norm = no_error_norm;
  NV_CONTENT_GKZ(nvout)->not_stepped = not_stepped;
  NV_CONTENT_GKZ(nvout)->arr = arr;
  NV_CONTENT_GKZ(nvout)->red_mem = &gsnv->red_mem;

  gsnv->nvec = nvout;
  return gsnv;
}

struct gkyl_array*
gkyl_sundials_nvec_get_array(struct gkyl_sundials_nvec *gsnv)
{
  N_Vector nvin = gsnv->nvec;
  return snvec_get_array(nvin);
}

void
gkyl_sundials_nvec_release(struct gkyl_sundials_nvec *gsnv)
{
  N_Vector nvin = gsnv->nvec;
  struct gkyl_array *arr = snvec_get_array(nvin);
  bool use_gpu = gkyl_array_is_cu_dev(arr);

  if (use_gpu) {
    gkyl_cu_free(gsnv->red_mem.red_local);
    gkyl_cu_free(gsnv->red_mem.red_global);
  }
  else {
    gkyl_free(gsnv->red_mem.red_local);
    gkyl_free(gsnv->red_mem.red_global);
  }
  gkyl_free(gsnv->red_mem.red_global_ho);

  snvec_destroy(nvin);
  gkyl_free(gsnv);
}

struct gkyl_sundials_nvec*
gkyl_sundials_many_nvec_new(struct gkyl_sundials *gksun, int num_nvector,
  struct gkyl_sundials_nvec *gsnv_arr[])
{
  struct gkyl_sundials_nvec *gsmanynv = gkyl_malloc(sizeof(*gsmanynv));
  
  // Create array of NVectors.
  N_Vector nvarr[num_nvector];
  for (int i=0; i<num_nvector; i++)
    nvarr[i] = gsnv_arr[i]->nvec;

  // Allocate the Sundials ManyNVector.
  N_Vector nvout = 0;
  nvout = N_VNew_ManyVector(num_nvector, nvarr, gksun->sunctx);
  if (nvout == 0) {
    fprintf(stderr, "Error: Creating new many nvector.\n");
    assert(false);
  }

  gsmanynv->nvec = nvout;

  return gsmanynv;
}

int 
gkyl_sundials_many_nvec_get_num_subvec(struct gkyl_sundials *gksun, struct gkyl_sundials_nvec* gsmanynv)
{
  return N_VGetNumSubvectors_ManyVector(gsmanynv->nvec);
}

void
gkyl_sundials_many_nvec_release(struct gkyl_sundials_nvec *gsmanynv)
{
  N_Vector nvin = gsmanynv->nvec;
  N_VDestroy(nvin);
  gkyl_free(gsmanynv);
}

#else

struct gkyl_sundials*
gkyl_sundials_new(bool use_gpu)
{
  fprintf(stderr, "\nGkeyll was not built with SUNDIALS. Use Gkeyll-native steppers or build Gkeyll with SUNDIALS.\n");
  return 0;
}

bool
gkyl_sundials_check_rk_method(enum gkyl_sundials_rk_method rk_method,
  enum gkyl_sundials_rk_method opt_A, enum gkyl_sundials_rk_method opt_B)
{
  return false;
};

void
gkyl_sundials_stepper_init(struct gkyl_sundials *gksun,
  struct gkyl_sundials_stepper_inp *inp)
{
  // Do nothing.
}

void
gkyl_sundials_arkode_reset(struct gkyl_sundials *gksun, double time,
  struct gkyl_sundials_nvec *gsmanynv, struct gkyl_sundials_nvec *gsmanynv_buff)
{
  // Do nothing.
}

void
gkyl_sundials_set_fixed_step(struct gkyl_sundials *gksun, double dt)
{
  // Do nothing.
}

void
gkyl_sundials_set_fixed_step_ssprk(struct gkyl_sundials *gksun, double dt)
{
  // Do nothing.
}

void
gkyl_sundials_set_fixed_step_sts(struct gkyl_sundials *gksun, double dt)
{
  // Do nothing.
}

double
gkyl_sundials_get_cfl_dt(struct gkyl_sundials *gksun)
{
  // Do nothing.
  return 0.0;
}

double
gkyl_sundials_get_cfl_dt_ssprk(struct gkyl_sundials *gksun)
{
  // Do nothing.
  return 0.0;
}

double
gkyl_sundials_get_cfl_dt_sts(struct gkyl_sundials *gksun)
{
  // Do nothing.
  return 0.0;
}

int
gkyl_sundials_evolve(struct gkyl_sundials *gksun, double t_new,
  struct gkyl_sundials_nvec *gsnv, double *t_curr)
{
  // Do nothing.
  return 0;
}

void
gkyl_sundials_release(struct gkyl_sundials *gksun)
{
  // Do nothing.
}

struct gkyl_sundials_nvec*
gkyl_sundials_nvec_new(struct gkyl_sundials *gksun, struct gkyl_array *arr,
  struct gkyl_comm *comm, struct gkyl_range *local_range, bool no_error_norm, bool not_stepped)
{
  // Do nothing.
  return 0;
}

struct gkyl_array*
gkyl_sundials_nvec_get_array(struct gkyl_sundials_nvec *gsnv)
{
  // Do nothing.
  return 0;
}

long
gkyl_sundials_get_num_error_test_failures(struct gkyl_sundials *gksun)
{
  // Do nothing.
  return 0;
}

long
gkyl_sundials_get_num_rhs_evals(struct gkyl_sundials *gksun)
{
  // Do nothing.
  return 0;
}

double
gkyl_sundials_get_last_dt(struct gkyl_sundials *gksun)
{
  // Do nothing.
  return 0;
}

double
gkyl_sundials_get_last_dt_ssprk(struct gkyl_sundials *gksun)
{
  // Do nothing.
  return 0;
}

double
gkyl_sundials_get_last_dt_sts(struct gkyl_sundials *gksun)
{
  // Do nothing.
  return 0;
}

double
gkyl_sundials_get_current_dt(struct gkyl_sundials *gksun)
{
  // Do nothing.
  return 0;
}

double
gkyl_sundials_get_current_dt_ssprk(struct gkyl_sundials *gksun)
{
  // Do nothing.
  return 0;
}

double
gkyl_sundials_get_current_dt_sts(struct gkyl_sundials *gksun)
{
  // Do nothing.
  return 0;
}

bool
gkyl_sundials_use_operator_split(struct gkyl_sundials *gksun)
{
  // Do nothing.
  return false;
}

bool
gkyl_sundials_operator_split_in_method(enum gkyl_sundials_rk_method rk_method)
{
  // Do nothing.
  return false;
}

void
gkyl_sundials_nvec_release(struct gkyl_sundials_nvec *gsnv)
{
  // Do nothing.
}

struct gkyl_sundials_nvec*
gkyl_sundials_many_nvec_new(struct gkyl_sundials *gksun, int num_nvector,
  struct gkyl_sundials_nvec *gsnv_arr[])
{
  // Do nothing.
  return 0;
}


int 
gkyl_sundials_many_nvec_get_num_subvec(struct gkyl_sundials *gksun, struct gkyl_sundials_nvec* gsmanynv)
{
  // Do nothing.
  return 0;
}

void
gkyl_sundials_many_nvec_release(struct gkyl_sundials_nvec *gsmanynv)
{
  // Do nothing.
}
#endif
