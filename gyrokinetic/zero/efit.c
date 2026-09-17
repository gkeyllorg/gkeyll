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

gkyl_efit* gkyl_efit_new(const struct gkyl_efit_inp *inp)
{
  gkyl_efit *up = gkyl_calloc(1, sizeof(struct gkyl_efit));

  up->reflect = inp->reflect;
  up->use_gpu = inp->use_gpu;
  up->filepath = inp->filepath;
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
  // Done, don't care about the rest
  
  fclose(ptr);

  int num_max_xpts = 10;
  double Rxpt[num_max_xpts];
  double Zxpt[num_max_xpts];

  up->num_xpts = find_xpts(up, Rxpt, Zxpt);
  up->Rxpt = gkyl_malloc(sizeof(double)*fmax(2, up->num_xpts));
  up->Zxpt = gkyl_malloc(sizeof(double)*fmax(2, up->num_xpts));
  for (int i = 0; i < up->num_xpts; i++) {
    up->Rxpt[i] = Rxpt[i];
    up->Zxpt[i] = Zxpt[i];
    // AS 9/24/24 This commented print statement is useful for checking the X-point Locations
    //  printf("Rxpt[%d] = %1.16f, Zxpt[%d] = %1.16f | psisep = %1.16f\n", i, up->Rxpt[i], i, up->Zxpt[i], up->psisep);
  }

  up->num_xpts_cubic = find_xpts_cubic(up, Rxpt, Zxpt);
  up->Rxpt_cubic = gkyl_malloc(sizeof(double)*fmax(2, up->num_xpts_cubic));
  up->Zxpt_cubic = gkyl_malloc(sizeof(double)*fmax(2, up->num_xpts_cubic));
  for (int i = 0; i < up->num_xpts_cubic; i++) {
    up->Rxpt_cubic[i] = Rxpt[i];
    up->Zxpt_cubic[i] = Zxpt[i];
    // AS 9/24/24 This commented print statement is useful for checking the X-point Locations
    //  printf("cubic: Rxpt[%d] = %1.16f, Zxpt[%d] = %1.16f | psisep = %1.16f\n", i, up->Rxpt_cubic[i], i, up->Zxpt_cubic[i], up->psisep_cubic);
  }

  return up;
}

void gkyl_efit_release(gkyl_efit* up){
  gkyl_free(up->Rxpt);
  gkyl_free(up->Zxpt);
  gkyl_free(up->Rxpt_cubic);
  gkyl_free(up->Zxpt_cubic);
  gkyl_array_release(up->psizr);
  gkyl_array_release(up->psizr_cubic);
  gkyl_array_release(up->bmagzr);
  gkyl_dg_basis_ops_evalf_release(up->evf);
  gkyl_array_release(up->fpolflux);
  gkyl_array_release(up->fpolprimeflux);
  gkyl_array_release(up->qflux);
  gkyl_free(up);
}
