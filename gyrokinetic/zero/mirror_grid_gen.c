#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_math.h>
#include <gkyl_mirror_grid_gen.h>
#include <gkyl_rect_decomp.h>

struct gkyl_mirror_grid_gen_x {
  enum gkyl_mirror_grid_gen_field_line_coord fl_coord; // field-line coordinate to use
  bool include_axis; // add nodes on r=0 axis (the axis is assumed be psi=0)
};

// context for use in root finder
struct psirz_ctx {
  struct gkyl_basis_ops_evalf *evcub; // cubic eval functions
  double Z; // local Z value
  double psi; // psi to match
};

static inline double
floor_sqrt(double x)
{
  return sqrt( fmax(x, 1e-14) );
}

static
double psirz(double R, void *ctx)
{
  struct psirz_ctx *rctx = ctx;
  double Z = rctx->Z;
  double xn[2] = { R, Z };
  double fout[1];
  rctx->evcub->eval_cubic(0, xn, fout, rctx->evcub->ctx);
  return fout[0] - rctx->psi;
}

static double calc_running_coord(double coord_lo, int i, double dx) {
  double dels[2] = {1.0/sqrt(3), 1.0-1.0/sqrt(3) };
  double coord = coord_lo;
  for(int j = 0; j < i; j++)
    coord+=dels[j%2]*dx;
  return coord;
}

static void
curlbhat_func(double r_curr, double Z, double phi, struct gkyl_basis_ops_evalf *evcub, struct gkyl_vec3 *curlbhat)
{
  // Calculate psi's various derivatives
  double Br = 0.0, Bz = 0.0, bmag = 0.0;
  double dpsidR = 0.0, dpsidZ = 0.0;
  double d2psidR2 = 0.0, d2psidZ2 = 0.0, d2psidRdZ = 0.0;
  double dBdR = 0.0, dBdZ = 0.0;
  double dBrdR = 0.0, dBrdZ = 0.0;
  double dBzdR = 0.0, dBzdZ = 0.0;

  double xn[2] = {r_curr, Z};
  double fout[4];
  evcub->eval_cubic_wgrad(0.0, xn, fout, evcub->ctx);
  dpsidR = fout[1];
  dpsidZ = fout[2];
  evcub->eval_cubic_wgrad2(0.0, xn, fout, evcub->ctx);
  d2psidR2 = fout[1];
  d2psidZ2 = fout[2];
  d2psidRdZ = fout[3];

  Br = 1.0/r_curr*dpsidZ;
  Bz = -1.0/r_curr*dpsidR;
  bmag = sqrt(Br*Br+Bz*Bz);

  dBrdR  = 1.0/r_curr*d2psidRdZ;
  dBrdZ  = 1.0/r_curr*d2psidZ2;
  dBzdR  = -1.0/r_curr*d2psidR2;
  dBzdZ  = -1.0/r_curr*d2psidRdZ;

  dBdR = 1/bmag*(Br*dBrdR + Bz*dBzdR);
  dBdZ = 1/bmag*(Br*dBrdZ + Bz*dBzdZ);

  // Get the polar components (contravariant, upperscript components on tangent basis)
  curlbhat->x[0] =  0.0; // R component ^1
  curlbhat->x[1] = 1.0/bmag*1.0/r_curr*(dBrdZ - dBzdR) + (-dBdR*Bz/r_curr + dBdZ*Br/r_curr); // Phi component ^2
  curlbhat->x[2] = 0.0;
}


struct gkyl_mirror_grid_gen *
gkyl_mirror_grid_gen_inew(const struct gkyl_mirror_grid_gen_inp *inp)
{
  struct gkyl_mirror_grid_gen *geo = gkyl_malloc(sizeof *geo);
  geo->gg_x = gkyl_malloc(sizeof *geo->gg_x);

  geo->gg_x->fl_coord = inp->fl_coord;
  geo->gg_x->include_axis = inp->include_axis;

  int nr = inp->nrcells, nz = inp->nzcells;
  int cells[] = { nr, nz };
  double lower[2] = { inp->R[0], inp->Z[0] };
  double upper[2] = { inp->R[1], inp->Z[1] };

  struct gkyl_rect_grid gridRZ;
  gkyl_rect_grid_init(&gridRZ, 2, lower, upper, cells);

  struct gkyl_basis_ops_evalf *evcub =
    gkyl_dg_basis_ops_evalf_new(&gridRZ, inp->psiRZ);

  do {
    const char *fname = inp->psi_cubic_fname ? inp->psi_cubic_fname : "psi_cubic.gkyl";
    if (inp->write_psi_cubic)
      gkyl_dg_basis_ops_evalf_write_cubic(evcub, fname);
  } while (0);

  // Construct grid in RZ plane
  
  enum { NPSI, NAL, NZ };
  
  geo->nodes_rza = gkyl_array_new(GKYL_DOUBLE, 3, inp->nrange.volume);
  geo->nodes_psi = gkyl_array_new(GKYL_DOUBLE, 1, inp->nrange.volume);
  geo->nodes_geom = gkyl_array_new(GKYL_USER, sizeof(struct gkyl_mirror_grid_gen_geom), inp->nrange.volume);


  // Take special care for setting dpsi
  double psi_lo = inp->comp_grid->lower[0];
  double psi_up = inp->comp_grid->upper[0];
  bool inc_axis = inp->include_axis;
  // Adjust if we are using sqrt(psi) as radial coordinate
  double psic_lo = psi_lo, psic_up = psi_up;
  if (inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z) {
    // When we include axis, psi_lo is ignored
    psic_lo = inc_axis ? 0.0 : sqrt(psi_lo);
    psic_up = sqrt(psi_up);
  }
  double dpsi = (psic_up-psic_lo)/(gkyl_range_shape(&inp->global,0));

  double dz = inp->comp_grid->dx[NZ];
  double dalpha = inp->comp_grid->dx[NAL];

  double z_lo = inp->comp_grid->lower[NZ] + (inp->local.lower[NZ] - inp->global.lower[NZ])*inp->comp_grid->dx[NZ];
  double alpha_lo = inp->comp_grid->lower[NAL] + (inp->local.lower[NAL] - inp->global.lower[NAL])*inp->comp_grid->dx[NAL];
  psi_lo = inp->comp_grid->lower[NPSI] + (inp->local.lower[NPSI] - inp->global.lower[NPSI])*inp->comp_grid->dx[NPSI];
  psi_up = inp->comp_grid->lower[NPSI] + (inp->local.upper[NPSI] - inp->global.lower[NPSI])*inp->comp_grid->dx[NPSI];

  // Adjust if we are using sqrt(psi) as radial coordinate
  psic_lo = psi_lo;
  psic_up = psi_up;
  if (inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z) {
    // When we include axis, psi_lo is ignored
    psic_lo = inc_axis ? 0.0 : sqrt(psi_lo);
    psic_up = sqrt(psi_up);
  }

  double rlow = lower[0], rup = upper[0];
  double rmin = rlow + 1e-8*(rup-rlow);

  struct psirz_ctx pctx = { .evcub = evcub };

  // Compute node locations
  bool status = true;
  for (int iz=inp->nrange.lower[NZ]; iz<=inp->nrange.upper[NZ]; ++iz) {
    double zcurr = z_lo + iz*dz;

    double Zcurr;
    inp->position_map->maps[2](0.0, &zcurr,  &Zcurr,  inp->position_map->ctxs[2]);
    zcurr = Zcurr; // update zcurr to be the mapped value

    double psi_min[1], psi_max[1];
    evcub->eval_cubic(0.0, (double[2]) { rmin, zcurr }, psi_min, evcub->ctx);
    evcub->eval_cubic(0.0, (double[2]) { rup, zcurr }, psi_max, evcub->ctx);

    for (int ipsi=inp->nrange.lower[NPSI]; ipsi<=inp->nrange.upper[NPSI]; ++ipsi) {

      if (inc_axis && (ipsi == 0)) {
        for (int ia=inp->nrange.lower[NAL]; ia<=inp->nrange.upper[NAL]; ++ia){
          double alpha_curr = alpha_lo + ia*dalpha;
          int idx[3] = { ipsi, ia, iz };
          double *rz = gkyl_array_fetch(geo->nodes_rza, gkyl_range_idx(&inp->nrange, idx));
          rz[0] = 0.0; rz[1] = zcurr; rz[2] = alpha_curr;
        }
      }
      else {
        double psic_curr = psic_lo + ipsi*dpsi;
        pctx.Z = zcurr;

        // We continue to do root-finding for psi and not sqrt(psi)
        double psi_curr = psic_curr;
        if (inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z)
          psi_curr = psic_curr*psic_curr;
        
        pctx.psi =  psi_curr;
        
        struct gkyl_qr_res root = gkyl_ridders(psirz, &pctx, rmin, rup,
          psi_min[0]-psi_curr, psi_max[0]-psi_curr,
          100, 1e-10);

        if (root.status) {
          status = false;
          goto cleanup;
        }
        
        for (int ia=inp->nrange.lower[NAL]; ia<=inp->nrange.upper[NAL]; ++ia){
          double alpha_curr = alpha_lo + ia*dalpha;
          int idx[3] = { ipsi, ia, iz };
          double *rz = gkyl_array_fetch(geo->nodes_rza, gkyl_range_idx(&inp->nrange, idx));
          rz[0] = root.res; rz[1] = zcurr; rz[2] = alpha_curr;
          double *psi_coord = gkyl_array_fetch(geo->nodes_psi, gkyl_range_idx(&inp->nrange, idx));
          psi_coord[0] = psi_curr; // store psi coordinate
        }
      }
    }
  }

  enum { PSI_I, DPSI_R_I, DPSI_Z_I };
  
  // Compute geometry at nodes
  for (int iz=inp->nrange.lower[NZ]; iz<=inp->nrange.upper[NZ]; ++iz) {
    
    for (int ipsi=inp->nrange.lower[NPSI]; ipsi<=inp->nrange.upper[NPSI]; ++ipsi) {

      for (int ia=inp->nrange.lower[NAL]; ia<=inp->nrange.upper[NAL]; ++ia){
        int idx[3] = { ipsi, ia, iz };
        long loc = gkyl_range_idx(&inp->nrange, idx);
        
        const double *rzp = gkyl_array_cfetch(geo->nodes_rza, loc);
        double rz[2] = { rzp[0], rzp[1] };

        double dZ_dz = gkyl_position_map_slope(inp->position_map, 2, rzp[1],
            dz, iz, &inp->nrange); // derivative of the position map in Z direction
        
        struct gkyl_mirror_grid_gen_geom *g = gkyl_array_fetch(geo->nodes_geom, loc);
        
        if (inc_axis && (ipsi == 0)) {
          double fout2[4]; // second derivative of psi is needed
          evcub->eval_cubic_wgrad2(0.0, rz, fout2, evcub->ctx);

          // On-axis the coordinate system breaks down. Below we choose
          // some reasonable defaults for the tnagent and
          // duals. However, the Jacobians and magnetic field are
          // correct and computed using the estimated asymptotic
          // behavior of psi as r -> 0.
          
          g->dual[0].x[0] = 1.0;
          g->dual[0].x[1] = 0.0;
          g->dual[0].x[2] = 0.0;

          g->dual[1].x[0] = 0;
          g->dual[1].x[1] = 1.0;
          g->dual[1].x[2] = 0.0;

          g->dual[2].x[0] = 0;
          g->dual[2].x[1] = 0.0;
          g->dual[2].x[2] = 1.0 / dZ_dz;        
          
          g->tang[0].x[0] = 1.0;
          g->tang[0].x[1] = 0.0;
          g->tang[0].x[2] = 0.0;

          g->tang[1].x[0] = 0.0;
          g->tang[1].x[1] = 1.0;
          g->tang[1].x[2] = 0.0;        
          
          g->tang[2].x[0] = 0;
          g->tang[2].x[1] = 0.0;
          g->tang[2].x[2] = 1.0 * dZ_dz;
          
          if (inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z)
            g->Jc = 0; // assumes asymptotics of psi ~ r^2 as r -> 0
          else
            g->Jc = 1/fout2[DPSI_R_I] * dZ_dz;

          g->B.x[0] = 0.0; // no radial component
          g->B.x[1] = 0.0;
          g->B.x[2] = fout2[DPSI_R_I]; // diff(psi,r,2)

          // \nabla X \hat{b} = 0 on axis because field is purely in Z direction
          g->curlbhat.x[0] = 0.0;
          g->curlbhat.x[1] = 0.0;
          g->curlbhat.x[2] = 0.0;
        }
        else {
          double fout[3]; // first derivative of psi is needed
          evcub->eval_cubic_wgrad(0.0, rz, fout, evcub->ctx);
        
          // e^1
          g->dual[0].x[0] = fout[DPSI_R_I]; // dpsi/dr
          g->dual[0].x[1] = 0.0; // no toroidal component
          g->dual[0].x[2] = fout[DPSI_Z_I]; // dspi/dz
        
          if (inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z) {
            // for sqrt(psi) as radial coordinate e^1 = grad(psi)/2*sqrt(psi)
            g->dual[0].x[0] = g->dual[0].x[0]/(2*floor_sqrt(fout[0]));
            g->dual[0].x[2] = g->dual[0].x[2]/(2*floor_sqrt(fout[0]));
          }

          // e^2 is just e^phi
          g->dual[1].x[0] = 0;
          g->dual[1].x[1] = 1.0/(rz[0]*rz[0]);
          g->dual[1].x[2] = 0.0;

          // e^3 is just sigma_3
          g->dual[2].x[0] = 0;
          g->dual[2].x[1] = 0.0;
          g->dual[2].x[2] = 1.0 / dZ_dz;

          // e_1 points along the radial direction
          g->tang[0].x[0] = 1/g->dual[0].x[0];
          g->tang[0].x[1] = 0.0;
          g->tang[0].x[2] = 0.0;

          // e_2
          g->tang[1].x[0] = 0;
          g->tang[1].x[1] = 1.0;
          g->tang[1].x[2] = 0.0;

          // e_3
          g->tang[2].x[0] = -fout[DPSI_Z_I]/fout[DPSI_R_I];
          g->tang[2].x[1] = 0.0;
          g->tang[2].x[2] = 1.0 * dZ_dz;
          
          if (inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z)
            g->Jc = 2*floor_sqrt(fout[PSI_I])*rz[0]/fout[DPSI_R_I] * dZ_dz;
          else
            g->Jc = rz[0]/fout[DPSI_R_I] * dZ_dz;
          
          g->B.x[0] = -fout[DPSI_Z_I]/rz[0];
          g->B.x[1] = 0.0;
          g->B.x[2] = fout[DPSI_R_I]/rz[0];

          // \nabla X \hat{b}
          curlbhat_func(rz[0], rz[1], rzp[2], evcub, &g->curlbhat);
        }
      }
    }
  }
  
  cleanup:

  if (true != status) {
    gkyl_mirror_grid_gen_release(geo);
    geo = 0;
    fprintf(stderr, "gkyl_mirror_grid_gen_inew failed to generate a grid\n");
  }
  
  gkyl_dg_basis_ops_evalf_release(evcub);
  
  return geo;
}

struct gkyl_mirror_grid_gen *
gkyl_mirror_grid_gen_int_inew(const struct gkyl_mirror_grid_gen_inp *inp)
{
  struct gkyl_mirror_grid_gen *geo = gkyl_malloc(sizeof *geo);
  geo->gg_x = gkyl_malloc(sizeof *geo->gg_x);

  geo->gg_x->fl_coord = inp->fl_coord;
  geo->gg_x->include_axis = inp->include_axis;

  int nr = inp->nrcells, nz = inp->nzcells;
  int cells[] = { nr, nz };
  double lower[2] = { inp->R[0], inp->Z[0] };
  double upper[2] = { inp->R[1], inp->Z[1] };

  struct gkyl_rect_grid gridRZ;
  gkyl_rect_grid_init(&gridRZ, 2, lower, upper, cells);

  struct gkyl_basis_ops_evalf *evcub =
    gkyl_dg_basis_ops_evalf_new(&gridRZ, inp->psiRZ);

  do {
    const char *fname = inp->psi_cubic_fname ? inp->psi_cubic_fname : "psi_cubic.gkyl";
    if (inp->write_psi_cubic)
      gkyl_dg_basis_ops_evalf_write_cubic(evcub, fname);
  } while (0);

  // Construct grid in RZ plane
  
  enum { NPSI, NAL, NZ };
  
  geo->nodes_rza = gkyl_array_new(GKYL_DOUBLE, 3, inp->nrange.volume);
  geo->nodes_psi = gkyl_array_new(GKYL_DOUBLE, 1, inp->nrange.volume);
  geo->nodes_geom = gkyl_array_new(GKYL_USER, sizeof(struct gkyl_mirror_grid_gen_geom), inp->nrange.volume);

  // Take special care for setting dpsi
  double psi_lo = inp->comp_grid->lower[0];
  double psi_up = inp->comp_grid->upper[0];
  bool inc_axis = inp->include_axis;
  // Adjust if we are using sqrt(psi) as radial coordinate
  double psic_lo = psi_lo, psic_up = psi_up;
  if (inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z) {
    // When we include axis, psi_lo is ignored
    psic_lo = inc_axis ? 0.0 : sqrt(psi_lo);
    psic_up = sqrt(psi_up);
  }
  double dpsi = (psic_up-psic_lo)/(gkyl_range_shape(&inp->global,0));

  double dz = inp->comp_grid->dx[NZ];
  double dalpha = inp->comp_grid->dx[NAL];

  double z_lo = inp->comp_grid->lower[NZ] + (inp->local.lower[NZ] - inp->global.lower[NZ])*inp->comp_grid->dx[NZ];
  double alpha_lo = inp->comp_grid->lower[NAL] + (inp->local.lower[NAL] - inp->global.lower[NAL])*inp->comp_grid->dx[NAL];
  psi_lo = inp->comp_grid->lower[NPSI] + (inp->local.lower[NPSI] - inp->global.lower[NPSI])*inp->comp_grid->dx[NPSI];
  psi_up = inp->comp_grid->lower[NPSI] + (inp->local.upper[NPSI] - inp->global.lower[NPSI])*inp->comp_grid->dx[NPSI];

  // Adjust if we are using sqrt(psi) as radial coordinate
  psic_lo = psi_lo;
  psic_up = psi_up;
  if (inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z) {
    // When we include axis, psi_lo is ignored
    psic_lo = sqrt(psi_lo);
    psic_up = sqrt(psi_up);
  }

  double dels[2] = {1.0/sqrt(3), 1.0-1.0/sqrt(3) };
  z_lo = z_lo + dels[1]*dz/2.0;
  psic_lo = psic_lo + dels[1]*dpsi/2.0;
  alpha_lo = alpha_lo + dels[1]*dalpha/2.0;

  double rlow = lower[0], rup = upper[0];
  double rmin = rlow + 1e-8*(rup-rlow);

  struct psirz_ctx pctx = { .evcub = evcub };

  // Compute node locations
  bool status = true;
  for (int iz=inp->nrange.lower[NZ]; iz<=inp->nrange.upper[NZ]; ++iz) {
    double zcurr = calc_running_coord(z_lo, iz-inp->nrange.lower[NZ], dz);

    double Zcurr;
    inp->position_map->maps[2](0.0, &zcurr,  &Zcurr,  inp->position_map->ctxs[2]);
    zcurr = Zcurr; // update zcurr to be the mapped value

    double psi_min[1], psi_max[1];
    evcub->eval_cubic(0.0, (double[2]) { rmin, zcurr }, psi_min, evcub->ctx);
    evcub->eval_cubic(0.0, (double[2]) { rup, zcurr }, psi_max, evcub->ctx);

    for (int ipsi=inp->nrange.lower[NPSI]; ipsi<=inp->nrange.upper[NPSI]; ++ipsi) {

      double psic_curr = calc_running_coord(psic_lo, ipsi-inp->nrange.lower[NPSI], dpsi);
      pctx.Z = zcurr;

      // We continue to do root-finding for psi and not sqrt(psi)
      double psi_curr = psic_curr;
      if (inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z)
        psi_curr = psic_curr*psic_curr;
      
      pctx.psi =  psi_curr;
      
      struct gkyl_qr_res root = gkyl_ridders(psirz, &pctx, rmin, rup,
        psi_min[0]-psi_curr, psi_max[0]-psi_curr,
        100, 1e-10);

      if (root.status) {
        status = false;
        goto cleanup;
      }
      
      for (int ia=inp->nrange.lower[NAL]; ia<=inp->nrange.upper[NAL]; ++ia){
        double alpha_curr = calc_running_coord(alpha_lo, ia-inp->nrange.lower[NAL], dalpha);
        int idx[3] = { ipsi, ia, iz };
        double *rz = gkyl_array_fetch(geo->nodes_rza, gkyl_range_idx(&inp->nrange, idx));
        rz[0] = root.res; rz[1] = zcurr; rz[2] = alpha_curr;
        double *psi_coord = gkyl_array_fetch(geo->nodes_psi, gkyl_range_idx(&inp->nrange, idx));
        psi_coord[0] = psi_curr; // store psi coordinate
      }
    }
  }

  enum { PSI_I, DPSI_R_I, DPSI_Z_I };
  
  // Compute geometry at nodes
  for (int iz=inp->nrange.lower[NZ]; iz<=inp->nrange.upper[NZ]; ++iz) {
    
    for (int ipsi=inp->nrange.lower[NPSI]; ipsi<=inp->nrange.upper[NPSI]; ++ipsi) {

      for (int ia=inp->nrange.lower[NAL]; ia<=inp->nrange.upper[NAL]; ++ia){
        int idx[3] = { ipsi, ia, iz };
        long loc = gkyl_range_idx(&inp->nrange, idx);
        
        const double *rzp = gkyl_array_cfetch(geo->nodes_rza, loc);
        double rz[2] = { rzp[0], rzp[1] };

        double dZ_dz = gkyl_position_map_slope(inp->position_map, 2, rzp[1],
            dz, iz, &inp->nrange); // derivative of the position map in Z direction
        
        struct gkyl_mirror_grid_gen_geom *g = gkyl_array_fetch(geo->nodes_geom, loc);
        
        double fout[3]; // first derivative of psi is needed
        evcub->eval_cubic_wgrad(0.0, rz, fout, evcub->ctx);
        
        // e^1
        g->dual[0].x[0] = fout[DPSI_R_I]; // dpsi/dr
        g->dual[0].x[1] = 0.0; // no toroidal component
        g->dual[0].x[2] = fout[DPSI_Z_I]; // dspi/dz
        
        if (inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z) {
          // for sqrt(psi) as radial coordinate e^1 = grad(psi)/2*sqrt(psi)
          g->dual[0].x[0] = g->dual[0].x[0]/(2*floor_sqrt(fout[0]));
          g->dual[0].x[2] = g->dual[0].x[2]/(2*floor_sqrt(fout[0]));
        }

        // e^2 is just e^phi
        g->dual[1].x[0] = 0;
        g->dual[1].x[1] = 1.0/(rz[0]*rz[0]);
        g->dual[1].x[2] = 0.0;

        // e^3 is just sigma_3
        g->dual[2].x[0] = 0;
        g->dual[2].x[1] = 0.0;
        g->dual[2].x[2] = 1.0 / dZ_dz;  

        // e_1 points along the radial direction
        g->tang[0].x[0] = 1/g->dual[0].x[0];
        g->tang[0].x[1] = 0.0;
        g->tang[0].x[2] = 0.0;

        // e_2
        g->tang[1].x[0] = 0;
        g->tang[1].x[1] = 1.0;
        g->tang[1].x[2] = 0.0;

        // e_3
        g->tang[2].x[0] = -fout[DPSI_Z_I]/fout[DPSI_R_I];
        g->tang[2].x[1] = 0.0;
        g->tang[2].x[2] = 1.0 * dZ_dz;
        
        if (inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z)
          g->Jc = 2*floor_sqrt(fout[PSI_I])*rz[0]/fout[DPSI_R_I] * dZ_dz;
        else
          g->Jc = rz[0]/fout[DPSI_R_I] * dZ_dz;
        
        g->B.x[0] = -fout[DPSI_Z_I]/rz[0];
        g->B.x[1] = 0.0;
        g->B.x[2] = fout[DPSI_R_I]/rz[0];

        // \nabla X \hat{b}
        curlbhat_func(rz[0], rz[1], rzp[2], evcub, &g->curlbhat);
      }
    }
  }
  
  cleanup:

  if (true != status) {
    gkyl_mirror_grid_gen_release(geo);
    geo = 0;
    fprintf(stderr, "gkyl_mirror_grid_gen_inew failed to generate a grid\n");
  }
  
  gkyl_dg_basis_ops_evalf_release(evcub);
  
  return geo;
}

struct gkyl_mirror_grid_gen *
gkyl_mirror_grid_gen_surf_inew(const struct gkyl_mirror_grid_gen_inp *inp)
{
  struct gkyl_mirror_grid_gen *geo = gkyl_malloc(sizeof *geo);
  geo->gg_x = gkyl_malloc(sizeof *geo->gg_x);

  geo->gg_x->fl_coord = inp->fl_coord;
  geo->gg_x->include_axis = inp->include_axis;

  int nr = inp->nrcells, nz = inp->nzcells;
  int cells[] = { nr, nz };
  double lower[2] = { inp->R[0], inp->Z[0] };
  double upper[2] = { inp->R[1], inp->Z[1] };

  struct gkyl_rect_grid gridRZ;
  gkyl_rect_grid_init(&gridRZ, 2, lower, upper, cells);

  struct gkyl_basis_ops_evalf *evcub =
    gkyl_dg_basis_ops_evalf_new(&gridRZ, inp->psiRZ);

  do {
    const char *fname = inp->psi_cubic_fname ? inp->psi_cubic_fname : "psi_cubic.gkyl";
    if (inp->write_psi_cubic)
      gkyl_dg_basis_ops_evalf_write_cubic(evcub, fname);
  } while (0);

  // Construct grid in RZ plane
  
  enum { NPSI, NAL, NZ };
  
  geo->nodes_rza = gkyl_array_new(GKYL_DOUBLE, 3, inp->nrange.volume);
  geo->nodes_psi = gkyl_array_new(GKYL_DOUBLE, 1, inp->nrange.volume);
  geo->nodes_geom = gkyl_array_new(GKYL_USER, sizeof(struct gkyl_mirror_grid_gen_geom), inp->nrange.volume);

  // Take special care for setting dpsi
  double psi_lo = inp->comp_grid->lower[0];
  double psi_up = inp->comp_grid->upper[0];
  bool inc_axis = inp->include_axis;
  // Adjust if we are using sqrt(psi) as radial coordinate
  double psic_lo = psi_lo, psic_up = psi_up;
  if (inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z) {
    // When we include axis, psi_lo is ignored
    psic_lo = inc_axis ? 0.0 : sqrt(psi_lo);
    psic_up = sqrt(psi_up);
  }
  double dpsi = (psic_up-psic_lo)/(gkyl_range_shape(&inp->global,0));

  double dz = inp->comp_grid->dx[NZ];
  double dalpha = inp->comp_grid->dx[NAL];

  double z_lo = inp->comp_grid->lower[NZ] + (inp->local.lower[NZ] - inp->global.lower[NZ])*inp->comp_grid->dx[NZ];
  double alpha_lo = inp->comp_grid->lower[NAL] + (inp->local.lower[NAL] - inp->global.lower[NAL])*inp->comp_grid->dx[NAL];
  psi_lo = inp->comp_grid->lower[NPSI] + (inp->local.lower[NPSI] - inp->global.lower[NPSI])*inp->comp_grid->dx[NPSI];
  psi_up = inp->comp_grid->lower[NPSI] + (inp->local.upper[NPSI] - inp->global.lower[NPSI])*inp->comp_grid->dx[NPSI];

  // Adjust if we are using sqrt(psi) as radial coordinate
  psic_lo = psi_lo;
  psic_up = psi_up;
  if (inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z) {
    // When we include axis, psi_lo is ignored
    psic_lo = inp->dir==0 && inc_axis ? 0.0 : sqrt(psi_lo);
    psic_up = sqrt(psi_up);
  }

  double dels[2] = {1.0/sqrt(3), 1.0-1.0/sqrt(3) };
  z_lo += inp->dir == 2 ? 0.0 : dels[1]*dz/2.0;
  psic_lo += inp->dir == 0 ? 0.0 : dels[1]*dpsi/2.0;
  alpha_lo += inp->dir == 1 ? 0. : dels[1]*dalpha/2.0;

  double rlow = lower[0], rup = upper[0];
  double rmin = rlow + 1e-8*(rup-rlow);

  struct psirz_ctx pctx = { .evcub = evcub };

  // Compute node locations
  bool status = true;
  for (int iz=inp->nrange.lower[NZ]; iz<=inp->nrange.upper[NZ]; ++iz) {
    double zcurr = calc_running_coord(z_lo, iz-inp->nrange.lower[NZ], dz);
    zcurr = inp->dir==2 ? z_lo + iz*dz: calc_running_coord(z_lo, iz-inp->nrange.lower[NZ], dz);

    double Zcurr;
    inp->position_map->maps[2](0.0, &zcurr,  &Zcurr,  inp->position_map->ctxs[2]);
    zcurr = Zcurr; // update zcurr to be the mapped value

    double psi_min[1], psi_max[1];
    evcub->eval_cubic(0.0, (double[2]) { rmin, zcurr }, psi_min, evcub->ctx);
    evcub->eval_cubic(0.0, (double[2]) { rup, zcurr }, psi_max, evcub->ctx);

    for (int ipsi=inp->nrange.lower[NPSI]; ipsi<=inp->nrange.upper[NPSI]; ++ipsi) {

      if (inc_axis && (ipsi == 0) && inp->dir==0) {
        for (int ia=inp->nrange.lower[NAL]; ia<=inp->nrange.upper[NAL]; ++ia){
          double alpha_curr = calc_running_coord(alpha_lo, ia-inp->nrange.lower[NAL], dalpha);
          int idx[3] = { ipsi, ia, iz };
          double *rz = gkyl_array_fetch(geo->nodes_rza, gkyl_range_idx(&inp->nrange, idx));
          rz[0] = 0.0; rz[1] = zcurr; rz[2] = alpha_curr;
        }
      }
      else {
        double psic_curr = inp->dir == 0 ? psic_lo + ipsi*dpsi : calc_running_coord(psic_lo, ipsi-inp->nrange.lower[NPSI], dpsi) ;
        pctx.Z = zcurr;

        // We continue to do root-finding for psi and not sqrt(psi)
        double psi_curr = psic_curr;
        if (inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z)
          psi_curr = psic_curr*psic_curr;
        
        pctx.psi =  psi_curr;
        
        struct gkyl_qr_res root = gkyl_ridders(psirz, &pctx, rmin, rup,
          psi_min[0]-psi_curr, psi_max[0]-psi_curr,
          100, 1e-10);

        if (root.status) {
          status = false;
          goto cleanup;
        }
        
        for (int ia=inp->nrange.lower[NAL]; ia<=inp->nrange.upper[NAL]; ++ia){
          double alpha_curr = inp->dir==1 ? alpha_lo + ia*dalpha : calc_running_coord(alpha_lo, ia-inp->nrange.lower[NAL], dalpha);
          int idx[3] = { ipsi, ia, iz };
          double *rz = gkyl_array_fetch(geo->nodes_rza, gkyl_range_idx(&inp->nrange, idx));
          rz[0] = root.res; rz[1] = zcurr; rz[2] = alpha_curr;
          double *psi_coord = gkyl_array_fetch(geo->nodes_psi, gkyl_range_idx(&inp->nrange, idx));
          psi_coord[0] = psi_curr; // store psi coordinate
        }
      }
    }
  }

  enum { PSI_I, DPSI_R_I, DPSI_Z_I };
  
  // Compute geometry at nodes
  for (int iz=inp->nrange.lower[NZ]; iz<=inp->nrange.upper[NZ]; ++iz) {
    
    for (int ipsi=inp->nrange.lower[NPSI]; ipsi<=inp->nrange.upper[NPSI]; ++ipsi) {

      for (int ia=inp->nrange.lower[NAL]; ia<=inp->nrange.upper[NAL]; ++ia){
        int idx[3] = { ipsi, ia, iz };
        long loc = gkyl_range_idx(&inp->nrange, idx);
        
        const double *rzp = gkyl_array_cfetch(geo->nodes_rza, loc);
        double rz[2] = { rzp[0], rzp[1] };

        double dZ_dz = gkyl_position_map_slope(inp->position_map, 2, rzp[1],
            dz, iz, &inp->nrange); // derivative of the position map in Z direction
        
        struct gkyl_mirror_grid_gen_geom *g = gkyl_array_fetch(geo->nodes_geom, loc);
        
        if (inc_axis && (ipsi == 0) && inp->dir==0) {
          double fout2[4]; // second derivative of psi is needed
          evcub->eval_cubic_wgrad2(0.0, rz, fout2, evcub->ctx);

          // On-axis the coordinate system breaks down. Below we choose
          // some reasonable defaults for the tnagent and
          // duals. However, the Jacobians and magnetic field are
          // correct and computed using the estimated asymptotic
          // behavior of psi as r -> 0.
          
          g->dual[0].x[0] = 1.0;
          g->dual[0].x[1] = 0.0;
          g->dual[0].x[2] = 0.0;

          g->dual[1].x[0] = 0;
          g->dual[1].x[1] = 1.0;
          g->dual[1].x[2] = 0.0;

          g->dual[2].x[0] = 0;
          g->dual[2].x[1] = 0.0;
          g->dual[2].x[2] = 1.0 / dZ_dz;

          g->tang[0].x[0] = 1.0;
          g->tang[0].x[1] = 0.0;
          g->tang[0].x[2] = 0.0;

          g->tang[1].x[0] = 0.0;
          g->tang[1].x[1] = 1.0;
          g->tang[1].x[2] = 0.0;        
          
          g->tang[2].x[0] = 0;
          g->tang[2].x[1] = 0.0;
          g->tang[2].x[2] = 1.0 * dZ_dz;
          
          if (inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z)
            g->Jc = 0; // assumes asymptotics of psi ~ r^2 as r -> 0
          else
            g->Jc = 1/fout2[DPSI_R_I] * dZ_dz;

          g->B.x[0] = 0.0; // no radial component
          g->B.x[1] = 0.0;
          g->B.x[2] = fout2[DPSI_R_I]; // diff(psi,r,2)

          // \nabla X \hat{b} = 0 on axis because field is purely in Z direction
          g->curlbhat.x[0] = 0.0;
          g->curlbhat.x[1] = 0.0;
          g->curlbhat.x[2] = 0.0;
        }
        else {
          double fout[3]; // first derivative of psi is needed
          evcub->eval_cubic_wgrad(0.0, rz, fout, evcub->ctx);
        
          // e^1
          g->dual[0].x[0] = fout[DPSI_R_I]; // dpsi/dr
          g->dual[0].x[1] = 0.0; // no toroidal component
          g->dual[0].x[2] = fout[DPSI_Z_I]; // dspi/dz
        
          if (inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z) {
            // for sqrt(psi) as radial coordinate e^1 = grad(psi)/2*sqrt(psi)
            g->dual[0].x[0] = g->dual[0].x[0]/(2*floor_sqrt(fout[0]));
            g->dual[0].x[2] = g->dual[0].x[2]/(2*floor_sqrt(fout[0]));
          }

          // e^2 is just e^phi
          g->dual[1].x[0] = 0;
          g->dual[1].x[1] = 1.0/(rz[0]*rz[0]);
          g->dual[1].x[2] = 0.0;

          // e^3 is just sigma_3
          g->dual[2].x[0] = 0;
          g->dual[2].x[1] = 0.0;
          g->dual[2].x[2] = 1.0 / dZ_dz;

          // e_1 points along the radial direction
          g->tang[0].x[0] = 1/g->dual[0].x[0];
          g->tang[0].x[1] = 0.0;
          g->tang[0].x[2] = 0.0;

          // e_2
          g->tang[1].x[0] = 0;
          g->tang[1].x[1] = 1.0;
          g->tang[1].x[2] = 0.0;

          // e_3
          g->tang[2].x[0] = -fout[DPSI_Z_I]/fout[DPSI_R_I];
          g->tang[2].x[1] = 0.0;
          g->tang[2].x[2] = 1.0 * dZ_dz;
          
          if (inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z)
            g->Jc = 2*floor_sqrt(fout[PSI_I])*rz[0]/fout[DPSI_R_I] * dZ_dz;
          else
            g->Jc = rz[0]/fout[DPSI_R_I] * dZ_dz;
          
          g->B.x[0] = -fout[DPSI_Z_I]/rz[0];
          g->B.x[1] = 0.0;
          g->B.x[2] = fout[DPSI_R_I]/rz[0];

          // \nabla X \hat{b}
          curlbhat_func(rz[0], rz[1], rzp[2], evcub, &g->curlbhat);
        }
      }
    }
  }
  
  cleanup:

  if (true != status) {
    gkyl_mirror_grid_gen_release(geo);
    geo = 0;
    fprintf(stderr, "gkyl_mirror_grid_gen_inew failed to generate a grid\n");
  }
  
  gkyl_dg_basis_ops_evalf_release(evcub);
  
  return geo;
}


bool
gkyl_mirror_grid_gen_is_include_axis(const struct gkyl_mirror_grid_gen *geom)
{
  return geom->gg_x->include_axis;
}

enum gkyl_mirror_grid_gen_field_line_coord
  gkyl_mirror_grid_gen_fl_coord(const struct gkyl_mirror_grid_gen *geom)
{
  return geom->gg_x->fl_coord;
}
  
void
gkyl_mirror_grid_gen_release(struct gkyl_mirror_grid_gen *geom)
{
  gkyl_array_release(geom->nodes_rza);
  gkyl_array_release(geom->nodes_psi);
  gkyl_array_release(geom->nodes_geom);
  gkyl_free(geom->gg_x);
  gkyl_free(geom);
}