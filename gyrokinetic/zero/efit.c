#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_efit.h>
#include <gkyl_efit_priv.h>

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_nodal_ops.h>
#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <stdint.h>

// Keep malformed-file rejection effective under -ffast-math.
static bool efit_finite(double value)
{
  uint64_t bits;
  memcpy(&bits,&value,sizeof bits);
  return (bits & UINT64_C(0x7ff0000000000000)) != UINT64_C(0x7ff0000000000000);
}


static double
efit_xpt_value(const double *arr, int n, int idx)
{
  return idx < n ? arr[idx] : DBL_MAX;
}

static void
write_xpt_diag_once(const struct gkyl_efit *up)
{
  static char last_filepath[1024] = { 0 };
  if (strncmp(last_filepath, up->filepath, sizeof(last_filepath)) == 0)
    return;

  strncpy(last_filepath, up->filepath, sizeof(last_filepath)-1);
  last_filepath[sizeof(last_filepath)-1] = '\0';

  double cubic_R0 = efit_xpt_value(up->Rxpt_cubic, up->num_xpts_cubic, 0);
  double cubic_Z0 = efit_xpt_value(up->Zxpt_cubic, up->num_xpts_cubic, 0);
  double cubic_R1 = efit_xpt_value(up->Rxpt_cubic, up->num_xpts_cubic, 1);
  double cubic_Z1 = efit_xpt_value(up->Zxpt_cubic, up->num_xpts_cubic, 1);
  double quad_R0 = efit_xpt_value(up->Rxpt, up->num_xpts, 0);
  double quad_Z0 = efit_xpt_value(up->Zxpt, up->num_xpts, 0);
  double quad_R1 = efit_xpt_value(up->Rxpt, up->num_xpts, 1);
  double quad_Z1 = efit_xpt_value(up->Zxpt, up->num_xpts, 1);

  fprintf(stderr,
    "XPT_DIAG name=%s filepath=%s reflect=%d cubic_n=%d cubic_psisep=%.16e "
    "cubic0=(%.16e,%.16e) cubic1=(%.16e,%.16e) quad_n=%d quad_psisep=%.16e "
    "quad0=(%.16e,%.16e) quad1=(%.16e,%.16e) raw_quad_found=%d "
    "raw_quad=(%.16e,%.16e) raw_quad_psi=%.16e raw_quad_dist_cell=%.16e "
    "fallback_to_cubic=%d\n",
    up->name, up->filepath, up->reflect, up->num_xpts_cubic, up->psisep_cubic,
    cubic_R0, cubic_Z0, cubic_R1, cubic_Z1, up->num_xpts, up->psisep,
    quad_R0, quad_Z0, quad_R1, quad_Z1, up->xpt_diag_quadratic_found,
    up->xpt_diag_quad_R, up->xpt_diag_quad_Z, up->xpt_diag_quad_psi,
    up->xpt_diag_quad_dist_cell, up->xpt_diag_fallback_to_cubic);

  const char *csv_path = getenv("GKYL_EFIT_XPT_DIAG_CSV");
  if (!csv_path || !csv_path[0])
    return;

  FILE *csv = fopen(csv_path, "a");
  if (!csv) {
    fprintf(stderr, "XPT_DIAG failed_to_open_csv=%s\n", csv_path);
    return;
  }

  fprintf(csv,
    "%s,%s,%d,%d,%.16e,%.16e,%.16e,%.16e,%.16e,%d,%.16e,%.16e,%.16e,%.16e,%.16e,%d,%.16e,%.16e,%.16e,%.16e,%d\n",
    up->name, up->filepath, up->reflect, up->num_xpts_cubic, up->psisep_cubic,
    cubic_R0, cubic_Z0, cubic_R1, cubic_Z1, up->num_xpts, up->psisep,
    quad_R0, quad_Z0, quad_R1, quad_Z1, up->xpt_diag_quadratic_found,
    up->xpt_diag_quad_R, up->xpt_diag_quad_Z, up->xpt_diag_quad_psi,
    up->xpt_diag_quad_dist_cell, up->xpt_diag_fallback_to_cubic);
  fclose(csv);
}

// Do any two non-adjacent edges of the closed outline properly cross?
//
// O(n^2) over n = 31..87 vertices, run once when an equilibrium is read, so the
// quadratic cost is irrelevant. Adjacent edges share an endpoint by
// construction and are skipped; zero-length edges (a repeated vertex, which is
// common and benign -- all 450 NSTX-U outlines have two) are skipped as well,
// since they cannot cross anything.
//
// A properly crossing pair means the crossing-number interior test is being
// asked a question that has no single right answer.
int
gkyl_efit_limiter_self_intersections(const struct gkyl_efit *e)
{
  if (!e || e->limiter_status != 1 || e->limiter_n < 4) return 0;
  const int n = e->limiter_n;
  int crossings = 0, overlaps = 0;
  for (int i=0; i<n; ++i) {
    double ax=e->limiter_R[i], ay=e->limiter_Z[i];
    double bx=e->limiter_R[(i+1)%n], by=e->limiter_Z[(i+1)%n];
    if (ax==bx && ay==by) continue;                    // zero-length edge
    for (int j=i+1; j<n; ++j) {
      // Skip the two edges that share a vertex with edge i, including the
      // wrap-around pair (i=0, j=n-1).
      if (j==i+1 || (i==0 && j==n-1)) continue;
      double cx=e->limiter_R[j], cy=e->limiter_Z[j];
      double dx=e->limiter_R[(j+1)%n], dy=e->limiter_Z[(j+1)%n];
      if (cx==dx && cy==dy) continue;
      double r1=(bx-ax)*(cy-ay)-(by-ay)*(cx-ax);
      double r2=(bx-ax)*(dy-ay)-(by-ay)*(dx-ax);
      double r3=(dx-cx)*(ay-cy)-(dy-cy)*(ax-cx);
      double r4=(dx-cx)*(by-cy)-(dy-cy)*(bx-cx);
      // A proper crossing: strict sign changes on both tests. The interior is
      // genuinely ambiguous here.
      if (((r1>0)!=(r2>0)) && ((r3>0)!=(r4>0))) { crossings++; continue; }
      // Collinear edges that OVERLAP over a positive length: the outline
      // retraces part of itself. This does not change the crossing-number
      // interior (retracing a line preserves parity), which is why the one
      // real instance measures as zero area and zero containment difference --
      // but it is still malformed input worth naming. Touching at a single
      // point is NOT counted: that is what a duplicate vertex produces, and
      // every one of the 450 NSTX-U outlines has two of those.
      if (r1==0.0 && r2==0.0 && r3==0.0 && r4==0.0) {
        double ux=bx-ax, uy=by-ay, uu=ux*ux+uy*uy;
        if (uu>0.0) {
          double tc=((cx-ax)*ux+(cy-ay)*uy)/uu;
          double td=((dx-ax)*ux+(dy-ay)*uy)/uu;
          double t0=tc<td?tc:td, t1=tc<td?td:tc;
          double lo=t0>0.0?t0:0.0, hi=t1<1.0?t1:1.0;
          if (hi-lo > 0.0) overlaps++;
        }
      }
    }
  }
  return crossings + overlaps;
}

gkyl_efit* gkyl_efit_new(const struct gkyl_efit_inp *inp)
{
  gkyl_efit *up = gkyl_calloc(1, sizeof(struct gkyl_efit));

  up->reflect = inp->reflect;
  up->use_gpu = inp->use_gpu;
  up->filepath = inp->filepath;

  // Copy the optional X-point search outline; the caller's arrays need not outlive us.
  up->xpt_bound_n = 0;
  up->xpt_bound_R = 0;
  up->xpt_bound_Z = 0;
  if (inp->xpt_bound_n > 0 && inp->xpt_bound_R && inp->xpt_bound_Z) {
    up->xpt_bound_n = inp->xpt_bound_n;
    up->xpt_bound_R = gkyl_malloc(sizeof(double)*inp->xpt_bound_n);
    up->xpt_bound_Z = gkyl_malloc(sizeof(double)*inp->xpt_bound_n);
    for (int i=0; i<inp->xpt_bound_n; ++i) {
      up->xpt_bound_R[i] = inp->xpt_bound_R[i];
      up->xpt_bound_Z[i] = inp->xpt_bound_Z[i];
    }
  }
  get_stripped_filename(up->filepath, up->name);

  gkyl_cart_modal_tensor(&up->rzbasis_cubic, 2, 3);
  gkyl_cart_modal_serendip(&up->fluxbasis, 1, inp->flux_poly_order);
  gkyl_cart_modal_tensor(&up->rzbasis, 2, inp->rz_poly_order);
  
  // Check if file exists using gkyl_check_file_exists and handle error only on rank 0.
  if (!gkyl_check_file_exists(up->filepath)) {
    fprintf(stderr, "efit.c: Failed to open the eqdsk file: %s\n", up->filepath);
    assert(false);
  }

  FILE *ptr = fopen(up->filepath,"r"); 

  // Read the last two ints in the first line, assuming they are N_R and N_Z.
  int MAX_LINE_LENGTH = 256;
  char first_line[MAX_LINE_LENGTH];
  char *token;
  if (fgets(first_line, sizeof(first_line), ptr) != NULL) {
    // Remove potential newline character from the end of the line.
    first_line[strcspn(first_line, "\n")] = 0;

    // Tokenize the string based on whitespace.
    char temp_line[MAX_LINE_LENGTH];
    strcpy(temp_line, first_line);
    token = strtok(temp_line, " \t\r\n");
    while (token != NULL) {
      // Check if the token is a number.
      int is_number = 1;
      for (int i = 0; i < strlen(token); i++) {
        if (!isdigit((unsigned char)token[i])) {
          is_number = 0;
          break;
        }
      }

      if (is_number) {
        // Store the last two integers.
        up->nr = up->nz;
        up->nz = atoi(token);
      }

      // Get the next token.
      token = strtok(NULL, " \t\r\n");
    }
  }

  // Read the non-array parameters, all are doubles:
  // rdim,zdim,rcentr,rleft,zmid;
  // rmaxis,zmaxis,simag,sibry,bcentr;
  // current,simag,xdum,rmaxis,xdum;
  // zmaxis,xdum,sibry,xdum,xdum;
  size_t status = fscanf(ptr,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
    &up->rdim, &up->zdim, &up->rcentr, &up->rleft, &up->zmid, &up-> rmaxis, &up->zmaxis, 
    &up->simag, &up->sibry, &up->bcentr, &up-> current, &up->simag, &up->xdum, &up->rmaxis, 
    &up->xdum, &up-> zmaxis, &up->xdum, &up->sibry, &up->xdum, &up->xdum);

  up->psisep = up->sibry;
  up->psisep_cubic = up->sibry;

  // Set zmid to 0 for double null
  if (up->reflect) {
    up->zmid = 0.0;
    up->zmaxis = 0.0;
  }


  // Now we need to make the grid
  up->zmin = up->zmid - up->zdim/2;
  up->zmax = up->zmid + up->zdim/2;
  up->rmin = up->rleft;
  up->rmax = up->rleft+up->rdim;

  double rzlower[2] = {up->rmin, up->zmin };
  double rzupper[2] = {up->rmax, up->zmax};
  int rzcells[2] = {0};
  int rzghost[2] = {1,1};
  if(up->rzbasis.poly_order==1){
    rzcells[0] = up->nr-1;
    rzcells[1]= up->nz-1;
  }
  if(up->rzbasis.poly_order==2){
    rzcells[0] = (up->nr-1)/2;
    rzcells[1] = (up->nz-1)/2;
  }
  gkyl_rect_grid_init(&up->rzgrid, 2, rzlower, rzupper, rzcells);
  gkyl_create_grid_ranges(&up->rzgrid, rzghost, &up->rzlocal_ext, &up->rzlocal);

  int cells_cubic[2] = {up->nr-1, up->nz-1};
  int rzghost_cubic[2] = {0,0};
  gkyl_rect_grid_init(&up->rzgrid_cubic, 2, rzlower, rzupper, cells_cubic);
  gkyl_create_grid_ranges(&up->rzgrid_cubic, rzghost_cubic, &up->rzlocal_cubic_ext, &up->rzlocal_cubic);

  double fluxlower[1];
  double fluxupper[1];
  bool step_convention; // True if psi increases toward magnetic axis
  if (up->simag > up->sibry) {
    step_convention = true;
    fluxlower[0] = up->sibry;
    fluxupper[0] = up->simag;
  }
  else {
    step_convention = false;
    fluxlower[0] = up->simag;
    fluxupper[0] = up->sibry;
  }

  int fluxcells[1] = {0};
  int fluxghost[2] = {1,1};
  if (up->fluxbasis.poly_order==1){
    fluxcells[0] = up->nr-1;
  }
  if (up->fluxbasis.poly_order==2){
    fluxcells[0] = (up->nr-1)/2;
  }

  gkyl_rect_grid_init(&up->fluxgrid, 1, fluxlower, fluxupper, fluxcells);
  gkyl_create_grid_ranges(&up->fluxgrid, fluxghost, &up->fluxlocal_ext, &up->fluxlocal);

  // allocate the necessary arrays
  up->psizr = gkyl_array_new(GKYL_DOUBLE, up->rzbasis.num_basis, up->rzlocal_ext.volume);
  up->psizr_cubic = gkyl_array_new(GKYL_DOUBLE, up->rzbasis_cubic.num_basis, up->rzlocal_cubic_ext.volume);
  up->bmagzr = gkyl_array_new(GKYL_DOUBLE, up->rzbasis.num_basis, up->rzlocal_ext.volume);
  up->fpolflux = gkyl_array_new(GKYL_DOUBLE, up->fluxbasis.num_basis, up->fluxlocal_ext.volume);
  up->fpolprimeflux = gkyl_array_new(GKYL_DOUBLE, up->fluxbasis.num_basis, up->fluxlocal_ext.volume);
  up->qflux = gkyl_array_new(GKYL_DOUBLE, up->fluxbasis.num_basis, up->fluxlocal_ext.volume);

  // Read fpol because we do want that
  int flux_node_nums[1] = {up->nr};
  struct gkyl_range flux_nrange;
  gkyl_range_init_from_shape(&flux_nrange, 1, flux_node_nums);
  struct gkyl_array *fpolflux_n = gkyl_array_new(GKYL_DOUBLE, 1, flux_nrange.volume);
  int fidx[1];
  // fpol is given on a uniform flux grid from the magnetic axis to plasma boundary
  if (step_convention) {
    for (int i = up->nr-1; i>=0; i--){
      fidx[0] = i;
      double *fpol_n= gkyl_array_fetch(fpolflux_n, gkyl_range_idx(&flux_nrange, fidx));
      status = fscanf(ptr,"%lf", fpol_n);
    }
  }
  else {
    for(int i = 0; i<up->nr; i++){
      fidx[0] = i;
      double *fpol_n= gkyl_array_fetch(fpolflux_n, gkyl_range_idx(&flux_nrange, fidx));
      status = fscanf(ptr,"%lf", fpol_n);
    }
  }

  struct gkyl_nodal_ops *n2m_flux = gkyl_nodal_ops_new(&up->fluxbasis, &up->fluxgrid, false);
  gkyl_nodal_ops_n2m(n2m_flux, &up->fluxbasis, &up->fluxgrid, 
    &flux_nrange, &up->fluxlocal, 1, fpolflux_n, up->fpolflux, false);

  // Now we have 3 of the 1d arrays, all of length nr :
  // pres, ffprim, pprime
  // I don't actually care about pres or pprime, so skip those

  //skip pres
  for(int i = 0; i<up->nr; i++){
    status = fscanf(ptr, "%lf", &up->xdum);
  }

  // read ffprime and divide out f
  struct gkyl_array *fpolprimeflux_n = gkyl_array_new(GKYL_DOUBLE, 1, flux_nrange.volume);
  // fpol*fpolprime is given on a uniform flux grid from the magnetic axis to plasma boundary
  if (step_convention) {
    for (int i = up->nr-1; i>=0; i--){
      fidx[0] = i;
      double *fpolprime_n = gkyl_array_fetch(fpolprimeflux_n, gkyl_range_idx(&flux_nrange, fidx));
      status = fscanf(ptr,"%lf", fpolprime_n);
      double *fpol_n = gkyl_array_fetch(fpolflux_n, gkyl_range_idx(&flux_nrange, fidx));
      fpolprime_n[0] = fpolprime_n[0]/fpol_n[0]; // divide out fpol
    }
  }
  else {
    for(int i = 0; i<up->nr; i++){
      fidx[0] = i;
      double *fpolprime_n= gkyl_array_fetch(fpolprimeflux_n, gkyl_range_idx(&flux_nrange, fidx));
      status = fscanf(ptr,"%lf", fpolprime_n);
      double *fpol_n = gkyl_array_fetch(fpolflux_n, gkyl_range_idx(&flux_nrange, fidx));
      fpolprime_n[0] = fpolprime_n[0]/fpol_n[0]; // divide out fpol
    }
  }
  gkyl_nodal_ops_n2m(n2m_flux, &up->fluxbasis, &up->fluxgrid, 
    &flux_nrange, &up->fluxlocal, 1, fpolprimeflux_n, up->fpolprimeflux, false);

  // skip pprime
  for(int i = 0; i<up->nr; i++){
    status = fscanf(ptr, "%lf", &up->xdum);
  }

  // Now we are gonna wanna read psi
  int node_nums[2] = {up->nr, up->nz};
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, up->rzgrid.ndim, node_nums);
  struct gkyl_array *psizr_n = gkyl_array_new(GKYL_DOUBLE, 1, nrange.volume);

  // Now lets loop through
  // Not only do we want psi at the nodes, we also want psi/R and psi/R^2 so we can use them for the magnetc field
  double R = up->rmin;
  double dR = up->rdim/(up->nr-1);
  int idx[2];
  for(int iz = 0; iz < up->nz; iz++){
    idx[1] = iz;
    for(int ir = 0; ir < up->nr; ir++){
      R = up->rmin+ir*dR;
      idx[0] = ir;
      // set psi
      double *psi_n = gkyl_array_fetch(psizr_n, gkyl_range_idx(&nrange, idx));
      status = fscanf(ptr,"%lf", psi_n);
    }
  }

  // We filled psizr_nodal
  struct gkyl_nodal_ops *n2m_rz = gkyl_nodal_ops_new(&up->rzbasis, &up->rzgrid, false);
  gkyl_nodal_ops_n2m(n2m_rz, &up->rzbasis, &up->rzgrid, &nrange, &up->rzlocal, 1, psizr_n, up->psizr, false);

  // Reflect psi for double null
  // Reflect DG coeffs rather than nodal data to avoid symmetry errors in n2m conversion
  if (up->reflect) {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &up->rzlocal);
    while (gkyl_range_iter_next(&iter)) {
      if (iter.idx[1] < gkyl_range_shape(&up->rzlocal,1)/2 +1 ) {
        int idx_change[2] = {iter.idx[0], gkyl_range_shape(&up->rzlocal, 1) - iter.idx[1]+1};
        const double *coeffs_ref = gkyl_array_cfetch(up->psizr, gkyl_range_idx(&up->rzlocal, iter.idx));
        double *coeffs  = gkyl_array_fetch(up->psizr, gkyl_range_idx(&up->rzlocal, idx_change));
        up->rzbasis.flip_odd_sign( 1, coeffs_ref, coeffs);
      }
    }
  }
 
  // Now lets read the q profile
  struct gkyl_array *qflux_n = gkyl_array_new(GKYL_DOUBLE, 1, flux_nrange.volume);
  int geqdsk_sign_convention = up->sibry > up->simag ? 0 : 1;
  if (geqdsk_sign_convention) {
    // psi increases toward magnetic axis.
    for (int i = up->nr-1; i>=0; i--) {
      fidx[0] = i;
      double *q_n= gkyl_array_fetch(qflux_n, gkyl_range_idx(&flux_nrange, fidx));
      status = fscanf(ptr, "%lf", q_n);
    }
  } else {
    // psi increases away from magnetic axis.
    for (int i = 0; i<up->nr; i++) {
      fidx[0] = i;
      double *q_n= gkyl_array_fetch(qflux_n, gkyl_range_idx(&flux_nrange, fidx));
      status = fscanf(ptr, "%lf", q_n);
    }
  }
  gkyl_nodal_ops_n2m(n2m_flux, &up->fluxbasis, &up->fluxgrid, 
    &flux_nrange, &up->fluxlocal, 1, qflux_n, up->qflux, false);


  // Standard EQDSK tail: plasma boundary count, limiter count, then R,Z pairs.
  // Older equilibrium files can omit the tail; only a requested extension
  // requires it. Never expose a partially read wall as usable geometry.
  int nboundary = 0, nlimiter = 0;
  int counts = fscanf(ptr, "%d %d", &nboundary, &nlimiter);
  if (counts != EOF) {
    up->limiter_status = -1;
    if (counts == 2 && nboundary >= 0 && nboundary <= 1000000 &&
        nlimiter >= 0 && nlimiter <= 1000000) {
      bool valid = true;
      double R, Z;
      for (int i=0; i<nboundary && valid; ++i)
        valid = fscanf(ptr, "%lf %lf", &R, &Z) == 2 && efit_finite(R) && efit_finite(Z);
      // Classify by vertex count exhaustively. Previously nlimiter==1 matched
      // neither branch and kept the pessimistic -1, so a readable record was
      // reported as corrupt (this is why step.geqdsk read as malformed); and
      // nlimiter==2 was accepted as status 1 even though every consumer
      // requires >=3, so the reader and its consumers disagreed.
      if (valid && nlimiter >= 1) {
        double *wall_R = gkyl_malloc(nlimiter*sizeof(double));
        double *wall_Z = gkyl_malloc(nlimiter*sizeof(double));
        for (int i=0; i<nlimiter && valid; ++i)
          valid = fscanf(ptr, "%lf %lf", &wall_R[i], &wall_Z[i]) == 2 &&
            efit_finite(wall_R[i]) && efit_finite(wall_Z[i]);
        if (valid && nlimiter >= 3) {
          up->limiter_n = nlimiter;
          up->limiter_R = wall_R; up->limiter_Z = wall_Z;
          up->limiter_status = 1;
          // Shape is diagnosed here and nowhere else. It does NOT feed
          // limiter_status: this reports, it does not decide.
          up->limiter_self_intersections = gkyl_efit_limiter_self_intersections(up);
          if (up->limiter_self_intersections > 0)
            fprintf(stderr,
              "TOK_GEO_WALL_MALFORMED_SHAPE name=%s vertices=%d self_intersections=%d\n",
              up->name, up->limiter_n, up->limiter_self_intersections);
        }
        else if (valid) {
          // Readable, but too few vertices to bound a region. Keep the count
          // for diagnostics; do not expose an unusable outline as geometry.
          up->limiter_n = nlimiter;
          gkyl_free(wall_R); gkyl_free(wall_Z);
          up->limiter_status = 2;
        }
        else { gkyl_free(wall_R); gkyl_free(wall_Z); }
      }
      else if (valid && nlimiter == 0)
        up->limiter_status = 0;
    }
  }

  // Make the cubic interpolator
  up->evf  = gkyl_dg_basis_ops_evalf_new(&up->rzgrid_cubic, psizr_n);
  gkyl_dg_basis_op_mem *mem = 0;
  mem = gkyl_dg_alloc_cubic_2d(cells_cubic);
  gkyl_dg_calc_cubic_2d_from_nodal_vals(mem, cells_cubic, up->rzgrid_cubic.dx, psizr_n, up->psizr_cubic);
  gkyl_dg_basis_op_mem_release(mem);

  // Calculate B.
  struct gkyl_array *bpolzr_n = gkyl_array_new(GKYL_DOUBLE, 1, nrange.volume);
  struct gkyl_array *bphizr_n = gkyl_array_new(GKYL_DOUBLE, 1, nrange.volume);
  struct gkyl_array *bmagzr_n = gkyl_array_new(GKYL_DOUBLE, 1, nrange.volume);
  double dZ = up->zdim/(up->nz-1);
  double scale_factorR = 2.0/(up->rzgrid_cubic.dx[0]);
  double scale_factorZ = 2.0/(up->rzgrid_cubic.dx[1]);
  for (int iz = 0; iz < up->nz; iz++){
    idx[1] = iz;
    double Z = up->zmin+iz*dZ;
    for (int ir = 0; ir < up->nr; ir++){
      R = up->rmin+ir*dR;
      idx[0] = ir;

      // Calculate Bpol.
      double xn[2] = {R, Z};
      double psi_curr, br, bz;
      if (R == 0.0) {
        double fout[4];
        up->evf->eval_cubic_wgrad2(0.0, xn, fout, up->evf->ctx);
        psi_curr = fout[0];
        br = fout[3];
        bz = -fout[1];
      }
      else {
        double fout[3];
        up->evf->eval_cubic_wgrad(0.0, xn, fout, up->evf->ctx);
        psi_curr = fout[0];
        br = 1.0/R*fout[2];
        bz = -1.0/R*fout[1];
      }
      double *bpol_n = gkyl_array_fetch(bpolzr_n, gkyl_range_idx(&nrange, idx));
      bpol_n[0] = sqrt(br*br + bz*bz);

      // Calculate Bphi.
      if (psi_curr < up->fluxgrid.lower[0] || psi_curr > up->fluxgrid.upper[0]){
        psi_curr = up->sibry;
      }
      fidx[0] = up->fluxlocal.lower[0] + (int) floor((psi_curr - up->fluxgrid.lower[0])/up->fluxgrid.dx[0]);
      fidx[0] = GKYL_MIN2(fidx[0], up->fluxlocal.upper[0]);
      fidx[0] = GKYL_MAX2(fidx[0], up->fluxlocal.lower[0]);
      long flux_loc = gkyl_range_idx(&up->fluxlocal, fidx);
      const double *coeffs = gkyl_array_cfetch(up->fpolflux, flux_loc);
      double fxc;
      gkyl_rect_grid_cell_center(&up->fluxgrid, fidx, &fxc);
      double fx = (psi_curr - fxc)/(up->fluxgrid.dx[0]*0.5);
      double fpol = up->fluxbasis.eval_expand(&fx, coeffs);
      double *bphi_n = gkyl_array_fetch(bphizr_n, gkyl_range_idx(&nrange, idx));
      if (fpol == 0.0 && R == 0.0)
        bphi_n[0] = 0.0;
      else 
        bphi_n[0] = fpol/R;

      // Calculate Bmag.
      double *bmag_n = gkyl_array_fetch(bmagzr_n, gkyl_range_idx(&nrange, idx));
      bmag_n[0] = sqrt(bpol_n[0]*bpol_n[0] + bphi_n[0]*bphi_n[0]);
    }
  }
  gkyl_nodal_ops_n2m(n2m_rz, &up->rzbasis, &up->rzgrid, &nrange, &up->rzlocal, 1, bmagzr_n, up->bmagzr, false);

  // Reflect B for double null.
  // Reflect DG coeffs rather than nodal data to avoid symmetry errors in n2m conversion.
  if (up->reflect) {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &up->rzlocal);
    while (gkyl_range_iter_next(&iter)) {
      if (iter.idx[1] < gkyl_range_shape(&up->rzlocal,1)/2 +1 ) {
        int idx_change[2] = {iter.idx[0], gkyl_range_shape(&up->rzlocal, 1) - iter.idx[1]+1};
        const double *coeffs_ref = gkyl_array_cfetch(up->bmagzr, gkyl_range_idx(&up->rzlocal, iter.idx));
        double *coeffs  = gkyl_array_fetch(up->bmagzr, gkyl_range_idx(&up->rzlocal, idx_change));
        up->rzbasis.flip_odd_sign( 1, coeffs_ref, coeffs);
      }
    }
  }
  
  // Free n2m operators
  gkyl_nodal_ops_release(n2m_flux);
  gkyl_nodal_ops_release(n2m_rz);
  // Free nodal arrays
  gkyl_array_release(fpolflux_n);
  gkyl_array_release(fpolprimeflux_n);
  gkyl_array_release(psizr_n);
  gkyl_array_release(qflux_n);
  gkyl_array_release(bpolzr_n);
  gkyl_array_release(bphizr_n);
  gkyl_array_release(bmagzr_n);
  // Done reading the equilibrium and its optional material boundary.
  
  fclose(ptr);

  int num_max_xpts = 10;
  double Rxpt[num_max_xpts];
  double Zxpt[num_max_xpts];

  up->num_xpts_cubic = find_xpts_cubic(up, Rxpt, Zxpt);
  up->Rxpt_cubic = gkyl_malloc(sizeof(double)*fmax(2, up->num_xpts_cubic));
  up->Zxpt_cubic = gkyl_malloc(sizeof(double)*fmax(2, up->num_xpts_cubic));
  for (int i = 0; i < up->num_xpts_cubic; i++) {
    up->Rxpt_cubic[i] = Rxpt[i];
    up->Zxpt_cubic[i] = Zxpt[i];
    // AS 9/24/24 This commented print statement is useful for checking the X-point Locations
    //  printf("cubic: Rxpt[%d] = %1.16f, Zxpt[%d] = %1.16f | psisep = %1.16f\n", i, up->Rxpt_cubic[i], i, up->Zxpt_cubic[i], up->psisep_cubic);
  }

  up->num_xpts = find_xpts(up, up->num_xpts_cubic, up->Rxpt_cubic, up->Zxpt_cubic, Rxpt, Zxpt);
  up->Rxpt = gkyl_malloc(sizeof(double)*fmax(2, up->num_xpts));
  up->Zxpt = gkyl_malloc(sizeof(double)*fmax(2, up->num_xpts));
  // Callers index Rxpt[0]/Zxpt[0] unconditionally.  If no X-point was found (possible
  // once an X-point search region is supplied), leave NaN rather than whatever the
  // allocation happened to contain, so the failure surfaces instead of propagating
  // silently as a garbage position.
  for (int i = 0; i < (int) fmax(2, up->num_xpts); ++i) {
    up->Rxpt[i] = NAN;
    up->Zxpt[i] = NAN;
  }
  if (up->num_xpts == 0)
    fprintf(stderr,
      "EFIT_XPT no X-point found inside the supplied search region name=%s bound_n=%d cubic_n=%d\n",
      up->name, up->xpt_bound_n, up->num_xpts_cubic);
  for (int i = 0; i < up->num_xpts; i++) {
    up->Rxpt[i] = Rxpt[i];
    up->Zxpt[i] = Zxpt[i];
    // AS 9/24/24 This commented print statement is useful for checking the X-point Locations
    //  printf("Rxpt[%d] = %1.16f, Zxpt[%d] = %1.16f | psisep = %1.16f\n", i, up->Rxpt[i], i, up->Zxpt[i], up->psisep);
  }

  write_xpt_diag_once(up);

  return up;
}

void
gkyl_efit_get_psi_bounds(const gkyl_efit *up, double *simag, double *psisep)
{
  if (simag)
    *simag = up->simag;
  if (psisep)
    *psisep = up->psisep;
}

void gkyl_efit_release(gkyl_efit* up){
  gkyl_free(up->Rxpt);
  gkyl_free(up->Zxpt);
  gkyl_free(up->Rxpt_cubic);
  gkyl_free(up->Zxpt_cubic);
  gkyl_free(up->limiter_R);
  gkyl_free(up->limiter_Z);
  gkyl_free(up->xpt_bound_R);
  gkyl_free(up->xpt_bound_Z);
  gkyl_array_release(up->psizr);
  gkyl_array_release(up->psizr_cubic);
  gkyl_array_release(up->bmagzr);
  gkyl_dg_basis_ops_evalf_release(up->evf);
  gkyl_array_release(up->fpolflux);
  gkyl_array_release(up->fpolprimeflux);
  gkyl_array_release(up->qflux);
  gkyl_free(up);
}
