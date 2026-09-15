#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_math.h>
#include <gkyl_gk_geometry_node_coords_priv.h>
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

  Br = -dpsidZ/r_curr;
  Bz = dpsidR/r_curr;
  bmag = sqrt(Br*Br+Bz*Bz);

  dBrdR  = -d2psidRdZ/r_curr - Br/r_curr;
  dBrdZ  = -d2psidZ2/r_curr;
  dBzdR  = d2psidR2/r_curr - Bz/r_curr;
  dBzdZ  = d2psidRdZ/r_curr;

  dBdR = 1/bmag*(Br*dBrdR + Bz*dBzdR);
  dBdZ = 1/bmag*(Br*dBrdZ + Bz*dBzdZ);

  // Get the polar components (contravariant, upperscript components on tangent basis)
  curlbhat->x[0] =  0.0; // R component ^1
  curlbhat->x[1] = (dBrdZ-dBzdR)/(bmag*r_curr) + (dBdR*Bz-dBdZ*Br)/(bmag*bmag*r_curr); // Phi component ^2
  curlbhat->x[2] = 0.0;
}


// kind = -1: corners, -2: volume quadrature, 0..2: surface quadrature.
// All three paths use the same coordinates and chain rule.
static struct gkyl_mirror_grid_gen *
mggen_new(const struct gkyl_mirror_grid_gen_inp *inp, int kind)
{
  struct gkyl_mirror_grid_gen *geo = gkyl_malloc(sizeof *geo);
  geo->gg_x = gkyl_malloc(sizeof *geo->gg_x);
  geo->gg_x->fl_coord = inp->fl_coord;
  geo->gg_x->include_axis = inp->include_axis;
  geo->nodes_rza = gkyl_array_new(GKYL_DOUBLE, 3, inp->nrange.volume);
  geo->nodes_psi = gkyl_array_new(GKYL_DOUBLE, 1, inp->nrange.volume);
  geo->nodes_geom = gkyl_array_new(GKYL_USER, sizeof(struct gkyl_mirror_grid_gen_geom), inp->nrange.volume);

  struct gkyl_rect_grid gridRZ;
  gkyl_rect_grid_init(&gridRZ, 2, (double[]) {inp->R[0], inp->Z[0]},
    (double[]) {inp->R[1], inp->Z[1]}, (int[]) {inp->nrcells, inp->nzcells});
  struct gkyl_basis_ops_evalf *evcub = gkyl_dg_basis_ops_evalf_new(&gridRZ, inp->psiRZ);
  if (inp->write_psi_cubic)
    gkyl_dg_basis_ops_evalf_write_cubic(evcub, inp->psi_cubic_fname ? inp->psi_cubic_fname : "psi_cubic.gkyl");

  // The input radial bounds are psi bounds; this option samples uniformly
  // in sqrt(psi), before applying any user radial position map.
  struct gkyl_rect_grid comp_grid = *inp->comp_grid;
  bool sqrt_psi = inp->fl_coord == GKYL_GEOMETRY_MIRROR_GRID_GEN_SQRT_PSI_CART_Z;
  if (sqrt_psi) {
    comp_grid.lower[0] = inp->include_axis ? 0.0 : sqrt(comp_grid.lower[0]);
    comp_grid.upper[0] = sqrt(comp_grid.upper[0]);
    comp_grid.dx[0] = (comp_grid.upper[0]-comp_grid.lower[0])/comp_grid.cells[0];
  }
  bool quad[3];
  int num_nodes[3];
  for (int dim=0; dim<3; ++dim) {
    quad[dim] = kind==-2 || (kind>=0 && kind!=dim);
    num_nodes[dim] = (gkyl_range_shape(&inp->nrange, dim)-(quad[dim] ? 0 : 1))
      /gkyl_range_shape(&inp->local, dim);
  }

  const struct gkyl_position_map *pmap = inp->position_map;
  double rmin = inp->R[0] + 1e-8*(inp->R[1]-inp->R[0]), rmax = inp->R[1];
  struct psirz_ctx pctx = {.evcub = evcub};
  bool status = true;
  for (int iz=inp->nrange.lower[2]; iz<=inp->nrange.upper[2]; ++iz) {
    double zcomp = gk_geometry_node_coord(&comp_grid, &inp->local, &inp->global,
      &inp->nrange, 2, iz, num_nodes[2], quad[2]);
    double zcurr;
    pmap->maps[2](0.0, &zcomp, &zcurr, pmap->ctxs[2]);
    // Use one computational cell as the nominal difference step, with
    // second-order inward differences near the global endpoints.
    double dZ_dz = gkyl_position_map_slope(pmap, 2, zcomp, comp_grid.dx[2],
      comp_grid.lower[2], comp_grid.upper[2]);
    double psi_min, psi_max;
    evcub->eval_cubic(0.0, (double[]) {rmin, zcurr}, &psi_min, evcub->ctx);
    evcub->eval_cubic(0.0, (double[]) {rmax, zcurr}, &psi_max, evcub->ctx);
    pctx.Z = zcurr;

    for (int ipsi=inp->nrange.lower[0]; ipsi<=inp->nrange.upper[0]; ++ipsi) {
      double psic = gk_geometry_node_coord(&comp_grid, &inp->local, &inp->global,
        &inp->nrange, 0, ipsi, num_nodes[0], quad[0]);
      double psi_map;
      pmap->maps[0](0.0, &psic, &psi_map, pmap->ctxs[0]);
      double dPsi_dpsi = gkyl_position_map_slope(pmap, 0, psic, comp_grid.dx[0],
        comp_grid.lower[0], comp_grid.upper[0]);
      double psi_curr = sqrt_psi ? psi_map*psi_map : psi_map;
      bool on_axis = inp->include_axis && !quad[0] && ipsi==inp->nrange.lower[0]
        && inp->local.lower[0]==inp->global.lower[0];
      double radius = 0.0;
      if (on_axis)
        psi_curr = 0.0;
      else {
        pctx.psi = psi_curr;
        struct gkyl_qr_res root = gkyl_ridders(psirz, &pctx, rmin, rmax,
          psi_min-psi_curr, psi_max-psi_curr, 100, 1e-10);
        if (root.status) { status = false; goto cleanup; }
        radius = root.res;
      }

      for (int ia=inp->nrange.lower[1]; ia<=inp->nrange.upper[1]; ++ia) {
        double alpha_comp = gk_geometry_node_coord(&comp_grid, &inp->local, &inp->global,
          &inp->nrange, 1, ia, num_nodes[1], quad[1]);
        double alpha_curr;
        pmap->maps[1](0.0, &alpha_comp, &alpha_curr, pmap->ctxs[1]);
        double dAlpha_dalpha = gkyl_position_map_slope(pmap, 1, alpha_comp, comp_grid.dx[1],
          comp_grid.lower[1], comp_grid.upper[1]);
        int idx[] = {ipsi, ia, iz};
        long loc = gkyl_range_idx(&inp->nrange, idx);
        double *rzp = gkyl_array_fetch(geo->nodes_rza, loc);
        rzp[0] = radius; rzp[1] = zcurr; rzp[2] = alpha_curr;
        double *psi_coord = gkyl_array_fetch(geo->nodes_psi, loc);
        psi_coord[0] = psi_curr;
        struct gkyl_mirror_grid_gen_geom *geom = gkyl_array_fetch(geo->nodes_geom, loc);
        *geom = (struct gkyl_mirror_grid_gen_geom) {0};
        double rz[] = {radius, zcurr};
        if (on_axis) {
          double fout[4];
          evcub->eval_cubic_wgrad2(0.0, rz, fout, evcub->ctx);
          // The radial/poloidal basis is singular at the axis. Keep the
          // existing defaults and the asymptotic psi ~ R^2 Jacobian.
          geom->dual[0].x[0] = 1.0/dPsi_dpsi;
          geom->dual[1].x[1] = 1.0/dAlpha_dalpha;
          geom->dual[2].x[2] = 1.0/dZ_dz;
          geom->tang[0].x[0] = dPsi_dpsi;
          geom->tang[1].x[1] = dAlpha_dalpha;
          geom->tang[2].x[2] = dZ_dz;
          geom->Jc = sqrt_psi ? 0.0 : dPsi_dpsi*dAlpha_dalpha*dZ_dz/fout[1];
          geom->B.x[2] = fout[1];
        }
        else {
          double fout[3];
          evcub->eval_cubic_wgrad(0.0, rz, fout, evcub->ctx);
          double radial_scale = dPsi_dpsi*(sqrt_psi ? 2.0*floor_sqrt(fout[0]) : 1.0);
          geom->dual[0].x[0] = fout[1]/radial_scale;
          geom->dual[0].x[2] = fout[2]/radial_scale;
          geom->dual[1].x[1] = 1.0/(radius*radius*dAlpha_dalpha);
          geom->dual[2].x[2] = 1.0/dZ_dz;
          geom->tang[0].x[0] = radial_scale/fout[1];
          geom->tang[1].x[1] = dAlpha_dalpha;
          // R(psi, M(z)): both components of the field-line tangent
          // carry dM/dz. This preserves e_3 dot grad(psi) = 0.
          geom->tang[2].x[0] = -fout[2]/fout[1]*dZ_dz;
          geom->tang[2].x[2] = dZ_dz;
          geom->Jc = radial_scale*radius/fout[1]*dAlpha_dalpha*dZ_dz;
          geom->B.x[0] = -fout[2]/radius;
          geom->B.x[2] = fout[1]/radius;
          curlbhat_func(radius, zcurr, alpha_curr, evcub, &geom->curlbhat);
        }
      }
    }
  }

cleanup:
  gkyl_dg_basis_ops_evalf_release(evcub);
  if (!status) {
    gkyl_mirror_grid_gen_release(geo);
    fprintf(stderr, "gkyl_mirror_grid_gen failed to generate a grid\n");
    return 0;
  }
  return geo;
}

struct gkyl_mirror_grid_gen *
gkyl_mirror_grid_gen_inew(const struct gkyl_mirror_grid_gen_inp *inp)
{
  return mggen_new(inp, -1);
}

struct gkyl_mirror_grid_gen *
gkyl_mirror_grid_gen_int_inew(const struct gkyl_mirror_grid_gen_inp *inp)
{
  return mggen_new(inp, -2);
}

struct gkyl_mirror_grid_gen *
gkyl_mirror_grid_gen_surf_inew(const struct gkyl_mirror_grid_gen_inp *inp)
{
  return mggen_new(inp, inp->dir);
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