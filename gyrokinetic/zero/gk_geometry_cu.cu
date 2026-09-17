/* -*- c++ -*- */
extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_math.h>
#include <gkyl_util.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_priv.h>
#include <assert.h>
}

__global__ static void
gk_geometry_set_corn_cu_kernel(struct gk_geometry *gk_geom,
  struct gkyl_array *mc2p, struct gkyl_array *mc2nu_pos, struct gkyl_array *bmag, struct gkyl_array *bmag_inv,
  struct gkyl_array *mc2p_nodal, struct gkyl_array *mc2nu_pos_nodal, struct gkyl_array *bmag_nodal,
  struct gkyl_array *mc2p_deflated, struct gkyl_array *mc2nu_pos_deflated
 )
{
  // Expansions.
  gk_geom->geo_corn.mc2p = mc2p;
  gk_geom->geo_corn.mc2nu_pos = mc2nu_pos;
  gk_geom->geo_corn.bmag = bmag;
  gk_geom->geo_corn.bmag_inv = bmag_inv;
  // Nodal.
  gk_geom->geo_corn.mc2p_nodal = mc2p_nodal;
  gk_geom->geo_corn.mc2nu_pos_nodal = mc2nu_pos_nodal;
  gk_geom->geo_corn.bmag_nodal = bmag_nodal;
  // Deflated.
  gk_geom->geo_corn.mc2p_deflated = mc2p_deflated;
  gk_geom->geo_corn.mc2nu_pos_deflated = mc2nu_pos_deflated;
}

__global__ static void
gk_geometry_set_int_cu_kernel(struct gk_geometry *gk_geom,
  struct gkyl_array *mc2p, struct gkyl_array *bmag, struct gkyl_array *g_ij, struct gkyl_array *g_ij_neut,
  struct gkyl_array *dxdz, struct gkyl_array *dzdx, struct gkyl_array *dualmag, struct gkyl_array *normals,
  struct gkyl_array *jacobgeo, struct gkyl_array *jacobgeo_ghost, struct gkyl_array *jacobgeo_inv, struct gkyl_array *gij,
  struct gkyl_array *gij_neut, struct gkyl_array *b_i, struct gkyl_array *bcart, struct gkyl_array *cmag,
  struct gkyl_array *jacobtot, struct gkyl_array *jacobtot_inv,
  struct gkyl_array *gxxj, struct gkyl_array *gxyj, struct gkyl_array *gyyj, struct gkyl_array *gxzj,
  struct gkyl_array *eps2, struct gkyl_array *rtg33inv, struct gkyl_array *dualcurlbhatoverB, struct gkyl_array *bioverJB,
  struct gkyl_array *B3, struct gkyl_array *dualcurlbhat, struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *mc2p_nodal,
  struct gkyl_array *bmag_nodal, struct gkyl_array *ddtheta_nodal, struct gkyl_array *ddpsi_nodal, struct gkyl_array *curlbhat_nodal,
  struct gkyl_array *dualcurlbhat_nodal, struct gkyl_array *jacobgeo_nodal, struct gkyl_array *g_ij_nodal,
  struct gkyl_array *g_ij_neut_nodal, struct gkyl_array *dxdz_nodal, struct gkyl_array *dzdx_nodal,
  struct gkyl_array *dualmag_nodal, struct gkyl_array *normals_nodal, struct gkyl_array *gij_neut_nodal,
  struct gkyl_array *b_i_nodal, struct gkyl_array *b_i_nodal_fd, struct gkyl_array *bcart_nodal,
  struct gkyl_array *B3_nodal, struct gkyl_array *dualcurlbhatoverB_nodal, struct gkyl_array *rtg33inv_nodal,
  struct gkyl_array *bioverJB_nodal, struct gkyl_array *qprofile
 )
{
  // Expansions.
  gk_geom->geo_int.mc2p = mc2p;
  gk_geom->geo_int.bmag = bmag;
  gk_geom->geo_int.g_ij = g_ij;
  gk_geom->geo_int.g_ij_neut = g_ij_neut;
  gk_geom->geo_int.dxdz = dxdz;
  gk_geom->geo_int.dzdx = dzdx;
  gk_geom->geo_int.dualmag = dualmag;
  gk_geom->geo_int.normals = normals;
  gk_geom->geo_int.jacobgeo = jacobgeo;
  gk_geom->geo_int.jacobgeo_ghost = jacobgeo_ghost;
  gk_geom->geo_int.jacobgeo_inv = jacobgeo_inv;
  gk_geom->geo_int.gij = gij;
  gk_geom->geo_int.gij_neut = gij_neut;
  gk_geom->geo_int.b_i = b_i;
  gk_geom->geo_int.bcart = bcart;
  gk_geom->geo_int.cmag = cmag;
  gk_geom->geo_int.jacobtot = jacobtot;
  gk_geom->geo_int.jacobtot_inv = jacobtot_inv;
  gk_geom->geo_int.gxxj = gxxj;
  gk_geom->geo_int.gxyj = gxyj;
  gk_geom->geo_int.gyyj = gyyj;
  gk_geom->geo_int.gxzj = gxzj;
  gk_geom->geo_int.eps2 = eps2;
  gk_geom->geo_int.rtg33inv = rtg33inv;
  gk_geom->geo_int.dualcurlbhatoverB = dualcurlbhatoverB;
  gk_geom->geo_int.dualcurlbhat = dualcurlbhat;
  gk_geom->geo_int.bioverJB = bioverJB;
  gk_geom->geo_int.B3 = B3;
  gk_geom->geo_int.qprofile = qprofile;
  // Nodal.
  gk_geom->geo_int.mc2p_nodal_fd = mc2p_nodal_fd;
  gk_geom->geo_int.mc2p_nodal = mc2p_nodal;
  gk_geom->geo_int.bmag_nodal = bmag_nodal;
  gk_geom->geo_int.ddtheta_nodal = ddtheta_nodal;
  gk_geom->geo_int.ddpsi_nodal = ddpsi_nodal;
  gk_geom->geo_int.curlbhat_nodal = curlbhat_nodal;
  gk_geom->geo_int.dualcurlbhat_nodal = dualcurlbhat_nodal;
  gk_geom->geo_int.jacobgeo_nodal = jacobgeo_nodal;
  gk_geom->geo_int.g_ij_nodal = g_ij_nodal;
  gk_geom->geo_int.g_ij_neut_nodal = g_ij_neut_nodal;
  gk_geom->geo_int.dxdz_nodal = dxdz_nodal;
  gk_geom->geo_int.dzdx_nodal = dzdx_nodal;
  gk_geom->geo_int.dualmag_nodal = dualmag_nodal;
  gk_geom->geo_int.normals_nodal = normals_nodal;
  gk_geom->geo_int.gij_neut_nodal = gij_neut_nodal;
  gk_geom->geo_int.b_i_nodal = b_i_nodal;
  gk_geom->geo_int.b_i_nodal_fd = b_i_nodal_fd;
  gk_geom->geo_int.bcart_nodal = bcart_nodal;
  gk_geom->geo_int.B3_nodal = B3_nodal;
  gk_geom->geo_int.dualcurlbhatoverB_nodal = dualcurlbhatoverB_nodal;
  gk_geom->geo_int.rtg33inv_nodal = rtg33inv_nodal;
  gk_geom->geo_int.bioverJB_nodal = bioverJB_nodal;
}

__global__ static void
gk_geometry_set_surf_cu_kernel(struct gk_geometry *gk_geom, int dir,
  struct gkyl_array *bmag, struct gkyl_array *jacobgeo, struct gkyl_array *jacobgeo_ratio,
  struct gkyl_array *b_i, struct gkyl_array *cmag, struct gkyl_array *jacobtot_inv,
  struct gkyl_array *B3, struct gkyl_array *normcurlbhat, struct gkyl_array *normals, struct gkyl_array *lenr,
  struct gkyl_array *bimpactangle, struct gkyl_array *deltats,
  struct gkyl_array *mc2p_nodal_fd, struct gkyl_array *mc2p_nodal, struct gkyl_array *bmag_nodal,
  struct gkyl_array *curlbhat_nodal, struct gkyl_array *normcurlbhat_nodal, struct gkyl_array *ddtheta_nodal,
  struct gkyl_array *ddpsi_nodal,
  struct gkyl_array *jacobgeo_nodal, struct gkyl_array *b_i_nodal, struct gkyl_array *b_i_nodal_fd,
  struct gkyl_array *cmag_nodal, struct gkyl_array *jacobtot_inv_nodal, struct gkyl_array *g_ij_nodal,
  struct gkyl_array *dxdz_nodal, struct gkyl_array *dzdx_nodal, struct gkyl_array *normals_nodal,
  struct gkyl_array *dualmag_nodal, struct gkyl_array *bcart_nodal, struct gkyl_array *B3_nodal,
  struct gkyl_array *lenr_nodal, struct gkyl_array *bimpactangle_nodal, struct gkyl_array *deltats_nodal
 )
{
  // Expansions.
  gk_geom->geo_surf[dir].bmag = bmag;
  gk_geom->geo_surf[dir].jacobgeo = jacobgeo;
  gk_geom->geo_surf[dir].jacobgeo_ratio = jacobgeo_ratio;
  gk_geom->geo_surf[dir].b_i = b_i;
  gk_geom->geo_surf[dir].cmag = cmag;
  gk_geom->geo_surf[dir].jacobtot_inv = jacobtot_inv;
  gk_geom->geo_surf[dir].B3 = B3;
  gk_geom->geo_surf[dir].normcurlbhat = normcurlbhat;
  gk_geom->geo_surf[dir].normals = normals;
  gk_geom->geo_surf[dir].lenr = lenr;
  gk_geom->geo_surf[dir].bimpactangle = bimpactangle;
  gk_geom->geo_surf[dir].deltats = deltats;
  // Nodal.
  gk_geom->geo_surf[dir].mc2p_nodal_fd = mc2p_nodal_fd;
  gk_geom->geo_surf[dir].mc2p_nodal = mc2p_nodal;
  gk_geom->geo_surf[dir].bmag_nodal = bmag_nodal;
  gk_geom->geo_surf[dir].curlbhat_nodal = curlbhat_nodal;
  gk_geom->geo_surf[dir].normcurlbhat_nodal = normcurlbhat_nodal;
  gk_geom->geo_surf[dir].ddtheta_nodal = ddtheta_nodal;
  gk_geom->geo_surf[dir].ddpsi_nodal = ddpsi_nodal;
  gk_geom->geo_surf[dir].jacobgeo_nodal = jacobgeo_nodal;
  gk_geom->geo_surf[dir].b_i_nodal = b_i_nodal;
  gk_geom->geo_surf[dir].b_i_nodal_fd = b_i_nodal_fd;
  gk_geom->geo_surf[dir].cmag_nodal = cmag_nodal;
  gk_geom->geo_surf[dir].jacobtot_inv_nodal = jacobtot_inv_nodal;
  gk_geom->geo_surf[dir].g_ij_nodal = g_ij_nodal;
  gk_geom->geo_surf[dir].dxdz_nodal = dxdz_nodal;
  gk_geom->geo_surf[dir].dzdx_nodal = dzdx_nodal;
  gk_geom->geo_surf[dir].normals_nodal = normals_nodal;
  gk_geom->geo_surf[dir].dualmag_nodal = dualmag_nodal;
  gk_geom->geo_surf[dir].bcart_nodal = bcart_nodal;
  gk_geom->geo_surf[dir].B3_nodal = B3_nodal;
  gk_geom->geo_surf[dir].lenr_nodal = lenr_nodal;
  gk_geom->geo_surf[dir].bimpactangle_nodal = bimpactangle_nodal;
  gk_geom->geo_surf[dir].deltats_nodal = deltats_nodal;
}

// Host-side wrapper for set_corn_cu_kernel
void
gkyl_geometry_set_corn_cu(struct gk_geometry *gk_geom, struct gk_geom_corn *geo_corn)
{
  gk_geometry_set_corn_cu_kernel<<<1,1>>>(gk_geom, 
    geo_corn->mc2p->on_dev, geo_corn->mc2nu_pos->on_dev, geo_corn->bmag->on_dev, geo_corn->bmag_inv->on_dev,
    geo_corn->mc2p_nodal->on_dev, geo_corn->mc2nu_pos_nodal->on_dev, geo_corn->bmag_nodal->on_dev,
    geo_corn->mc2p_deflated->on_dev, geo_corn->mc2nu_pos_deflated->on_dev);
}

// Host-side wrapper for set_int_cu_kernel
void
gkyl_geometry_set_int_cu(struct gk_geometry *gk_geom, struct gk_geom_int *geo_int)
{
  gk_geometry_set_int_cu_kernel<<<1,1>>>(gk_geom,
    geo_int->mc2p->on_dev, geo_int->bmag->on_dev, geo_int->g_ij->on_dev, geo_int->g_ij_neut->on_dev,
    geo_int->dxdz->on_dev, geo_int->dzdx->on_dev, geo_int->dualmag->on_dev, geo_int->normals->on_dev,
    geo_int->jacobgeo->on_dev, geo_int->jacobgeo_ghost->on_dev, geo_int->jacobgeo_inv->on_dev, geo_int->gij->on_dev,
    geo_int->gij_neut->on_dev, geo_int->b_i->on_dev, geo_int->bcart->on_dev, geo_int->cmag->on_dev,
    geo_int->jacobtot->on_dev, geo_int->jacobtot_inv->on_dev,
    geo_int->gxxj->on_dev, geo_int->gxyj->on_dev, geo_int->gyyj->on_dev, geo_int->gxzj->on_dev,
    geo_int->eps2->on_dev, geo_int->rtg33inv->on_dev, geo_int->dualcurlbhatoverB->on_dev, geo_int->bioverJB->on_dev,
    geo_int->B3->on_dev, geo_int->dualcurlbhat->on_dev,
    geo_int->mc2p_nodal_fd->on_dev, geo_int->mc2p_nodal->on_dev, geo_int->bmag_nodal->on_dev,
    geo_int->ddtheta_nodal->on_dev, geo_int->ddpsi_nodal->on_dev, geo_int->curlbhat_nodal->on_dev, geo_int->dualcurlbhat_nodal->on_dev,
    geo_int->jacobgeo_nodal->on_dev, geo_int->g_ij_nodal->on_dev, geo_int->g_ij_neut_nodal->on_dev,
    geo_int->dxdz_nodal->on_dev, geo_int->dzdx_nodal->on_dev, geo_int->dualmag_nodal->on_dev,
    geo_int->normals_nodal->on_dev, geo_int->gij_neut_nodal->on_dev, geo_int->b_i_nodal->on_dev,
    geo_int->b_i_nodal_fd->on_dev, geo_int->bcart_nodal->on_dev, geo_int->B3_nodal->on_dev,
    geo_int->dualcurlbhatoverB_nodal->on_dev, geo_int->rtg33inv_nodal->on_dev, geo_int->bioverJB_nodal->on_dev,
    geo_int->qprofile->on_dev);
}

// Host-side wrapper for set_surf_cu_kernel
void
gkyl_geometry_set_surf_cu(struct gk_geometry *gk_geom, struct gk_geom_surf *geo_surf, int dir)
{
  gk_geometry_set_surf_cu_kernel<<<1,1>>>(gk_geom, dir,
    geo_surf->bmag->on_dev, geo_surf->jacobgeo->on_dev, geo_surf->jacobgeo_ratio->on_dev,
    geo_surf->b_i->on_dev, geo_surf->cmag->on_dev, geo_surf->jacobtot_inv->on_dev, geo_surf->B3->on_dev,
    geo_surf->normcurlbhat->on_dev, geo_surf->normals->on_dev, geo_surf->lenr->on_dev, geo_surf->bimpactangle->on_dev, geo_surf->deltats->on_dev,
    geo_surf->mc2p_nodal_fd->on_dev, geo_surf->mc2p_nodal->on_dev, geo_surf->bmag_nodal->on_dev,
    geo_surf->curlbhat_nodal->on_dev, geo_surf->normcurlbhat_nodal->on_dev, geo_surf->ddtheta_nodal->on_dev, geo_surf->ddpsi_nodal->on_dev,
    geo_surf->jacobgeo_nodal->on_dev, geo_surf->b_i_nodal->on_dev, geo_surf->b_i_nodal_fd->on_dev,
    geo_surf->cmag_nodal->on_dev, geo_surf->jacobtot_inv_nodal->on_dev, geo_surf->g_ij_nodal->on_dev,
    geo_surf->dxdz_nodal->on_dev, geo_surf->dzdx_nodal->on_dev, geo_surf->normals_nodal->on_dev,
    geo_surf->dualmag_nodal->on_dev, geo_surf->bcart_nodal->on_dev, geo_surf->B3_nodal->on_dev,
    geo_surf->lenr_nodal->on_dev, geo_surf->bimpactangle->on_dev, geo_surf->deltats_nodal->on_dev);
}

struct gk_geom_corn*
gk_geometry_corn_cu_dev_alloc(struct gk_geom_corn up_corn_host)
{
  struct gk_geom_corn *up_corn_dev = (struct gk_geom_corn*) gkyl_malloc(sizeof(struct gk_geom_corn));
  // Expansions.
  up_corn_dev->mc2p = gkyl_array_cu_dev_new(up_corn_host.mc2p->type,
    up_corn_host.mc2p->ncomp, up_corn_host.mc2p->size);
  up_corn_dev->mc2nu_pos = gkyl_array_cu_dev_new(up_corn_host.mc2nu_pos->type,
    up_corn_host.mc2nu_pos->ncomp, up_corn_host.mc2nu_pos->size);
  up_corn_dev->bmag = gkyl_array_cu_dev_new(up_corn_host.bmag->type,
    up_corn_host.bmag->ncomp, up_corn_host.bmag->size);
  up_corn_dev->bmag_inv = gkyl_array_cu_dev_new(up_corn_host.bmag_inv->type,
    up_corn_host.bmag_inv->ncomp, up_corn_host.bmag_inv->size);
  up_corn_dev->mc2p_deflated = gkyl_array_cu_dev_new(up_corn_host.mc2p_deflated->type,
    up_corn_host.mc2p_deflated->ncomp, up_corn_host.mc2p_deflated->size);
  up_corn_dev->mc2nu_pos_deflated = gkyl_array_cu_dev_new(up_corn_host.mc2nu_pos_deflated->type,
    up_corn_host.mc2nu_pos_deflated->ncomp, up_corn_host.mc2nu_pos_deflated->size);
  // Nodal.
  up_corn_dev->mc2p_nodal = gkyl_array_cu_dev_new(up_corn_host.mc2p_nodal->type,
    up_corn_host.mc2p_nodal->ncomp, up_corn_host.mc2p_nodal->size);
  up_corn_dev->mc2nu_pos_nodal = gkyl_array_cu_dev_new(up_corn_host.mc2nu_pos_nodal->type,
    up_corn_host.mc2nu_pos_nodal->ncomp, up_corn_host.mc2nu_pos_nodal->size);
  up_corn_dev->bmag_nodal = gkyl_array_cu_dev_new(up_corn_host.bmag_nodal->type,
    up_corn_host.bmag_nodal->ncomp, up_corn_host.bmag_nodal->size);
  return up_corn_dev;
}

struct gk_geom_int*
gk_geometry_int_cu_dev_alloc(struct gk_geom_int up_int_host)
{
  struct gk_geom_int *up_int_dev = (struct gk_geom_int*) gkyl_malloc(sizeof(struct gk_geom_int));
  // Expansions.
  up_int_dev->mc2p = gkyl_array_cu_dev_new(up_int_host.mc2p->type,
    up_int_host.mc2p->ncomp, up_int_host.mc2p->size);
  up_int_dev->bmag = gkyl_array_cu_dev_new(up_int_host.bmag->type,
    up_int_host.bmag->ncomp, up_int_host.bmag->size);
  up_int_dev->g_ij = gkyl_array_cu_dev_new(up_int_host.g_ij->type,
    up_int_host.g_ij->ncomp, up_int_host.g_ij->size);
  up_int_dev->g_ij_neut = gkyl_array_cu_dev_new(up_int_host.g_ij_neut->type,
    up_int_host.g_ij_neut->ncomp, up_int_host.g_ij_neut->size);
  up_int_dev->dxdz = gkyl_array_cu_dev_new(up_int_host.dxdz->type,
    up_int_host.dxdz->ncomp, up_int_host.dxdz->size);
  up_int_dev->dzdx = gkyl_array_cu_dev_new(up_int_host.dzdx->type,
    up_int_host.dzdx->ncomp, up_int_host.dzdx->size);
  up_int_dev->dualmag = gkyl_array_cu_dev_new(up_int_host.dualmag->type,
    up_int_host.dualmag->ncomp, up_int_host.dualmag->size);
  up_int_dev->normals = gkyl_array_cu_dev_new(up_int_host.normals->type,
    up_int_host.normals->ncomp, up_int_host.normals->size);
  up_int_dev->jacobgeo = gkyl_array_cu_dev_new(up_int_host.jacobgeo->type,
    up_int_host.jacobgeo->ncomp, up_int_host.jacobgeo->size);
  up_int_dev->jacobgeo_ghost = gkyl_array_cu_dev_new(up_int_host.jacobgeo_ghost->type,
    up_int_host.jacobgeo_ghost->ncomp, up_int_host.jacobgeo_ghost->size);
  up_int_dev->jacobgeo_inv = gkyl_array_cu_dev_new(up_int_host.jacobgeo_inv->type,
    up_int_host.jacobgeo_inv->ncomp, up_int_host.jacobgeo_inv->size);
  up_int_dev->gij = gkyl_array_cu_dev_new(up_int_host.gij->type,
    up_int_host.gij->ncomp, up_int_host.gij->size);
  up_int_dev->gij_neut = gkyl_array_cu_dev_new(up_int_host.gij_neut->type,
    up_int_host.gij_neut->ncomp, up_int_host.gij_neut->size);
  up_int_dev->b_i = gkyl_array_cu_dev_new(up_int_host.b_i->type,
    up_int_host.b_i->ncomp, up_int_host.b_i->size);
  up_int_dev->bcart = gkyl_array_cu_dev_new(up_int_host.bcart->type,
    up_int_host.bcart->ncomp, up_int_host.bcart->size);
  up_int_dev->cmag = gkyl_array_cu_dev_new(up_int_host.cmag->type, 
    up_int_host.cmag->ncomp, up_int_host.cmag->size);
  up_int_dev->jacobtot = gkyl_array_cu_dev_new(up_int_host.jacobtot->type,
    up_int_host.jacobtot->ncomp, up_int_host.jacobtot->size);
  up_int_dev->jacobtot_inv = gkyl_array_cu_dev_new(up_int_host.jacobtot_inv->type,
    up_int_host.jacobtot_inv->ncomp, up_int_host.jacobtot_inv->size);
  up_int_dev->gxxj = gkyl_array_cu_dev_new(up_int_host.gxxj->type, 
    up_int_host.gxxj->ncomp, up_int_host.gxxj->size);
  up_int_dev->gxyj = gkyl_array_cu_dev_new(up_int_host.gxyj->type, 
    up_int_host.gxyj->ncomp, up_int_host.gxyj->size);
  up_int_dev->gyyj = gkyl_array_cu_dev_new(up_int_host.gyyj->type, 
    up_int_host.gyyj->ncomp, up_int_host.gyyj->size);
  up_int_dev->gxzj = gkyl_array_cu_dev_new(up_int_host.gxzj->type, 
    up_int_host.gxzj->ncomp, up_int_host.gxzj->size);
  up_int_dev->eps2 = gkyl_array_cu_dev_new(up_int_host.eps2->type, 
    up_int_host.eps2->ncomp, up_int_host.eps2->size);
  up_int_dev->rtg33inv = gkyl_array_cu_dev_new(up_int_host.rtg33inv->type, 
    up_int_host.rtg33inv->ncomp, up_int_host.rtg33inv->size);
  up_int_dev->dualcurlbhatoverB = gkyl_array_cu_dev_new(up_int_host.dualcurlbhatoverB->type,
    up_int_host.dualcurlbhatoverB->ncomp, up_int_host.dualcurlbhatoverB->size);
  up_int_dev->bioverJB = gkyl_array_cu_dev_new(up_int_host.bioverJB->type, 
    up_int_host.bioverJB->ncomp, up_int_host.bioverJB->size);
  up_int_dev->B3 = gkyl_array_cu_dev_new(up_int_host.B3->type, 
    up_int_host.B3->ncomp, up_int_host.B3->size);
  up_int_dev->dualcurlbhat = gkyl_array_cu_dev_new(up_int_host.dualcurlbhat->type,
    up_int_host.dualcurlbhat->ncomp, up_int_host.dualcurlbhat->size);
  up_int_dev->qprofile = gkyl_array_cu_dev_new(up_int_host.qprofile->type,
    up_int_host.qprofile->ncomp, up_int_host.qprofile->size);
  // Nodal.
  up_int_dev->mc2p_nodal_fd = gkyl_array_cu_dev_new(up_int_host.mc2p_nodal_fd->type,
    up_int_host.mc2p_nodal_fd->ncomp, up_int_host.mc2p_nodal_fd->size);
  up_int_dev->mc2p_nodal = gkyl_array_cu_dev_new(up_int_host.mc2p_nodal->type,
    up_int_host.mc2p_nodal->ncomp, up_int_host.mc2p_nodal->size);
  up_int_dev->bmag_nodal = gkyl_array_cu_dev_new(up_int_host.bmag_nodal->type,
    up_int_host.bmag_nodal->ncomp, up_int_host.bmag_nodal->size);
  up_int_dev->ddtheta_nodal = gkyl_array_cu_dev_new(up_int_host.ddtheta_nodal->type,
    up_int_host.ddtheta_nodal->ncomp, up_int_host.ddtheta_nodal->size);
  up_int_dev->ddpsi_nodal = gkyl_array_cu_dev_new(up_int_host.ddpsi_nodal->type,
    up_int_host.ddpsi_nodal->ncomp, up_int_host.ddpsi_nodal->size);
  up_int_dev->curlbhat_nodal = gkyl_array_cu_dev_new(up_int_host.curlbhat_nodal->type,
    up_int_host.curlbhat_nodal->ncomp, up_int_host.curlbhat_nodal->size);
  up_int_dev->dualcurlbhat_nodal = gkyl_array_cu_dev_new(up_int_host.dualcurlbhat_nodal->type,
    up_int_host.dualcurlbhat_nodal->ncomp, up_int_host.dualcurlbhat_nodal->size);
  up_int_dev->jacobgeo_nodal = gkyl_array_cu_dev_new(up_int_host.jacobgeo_nodal->type,
    up_int_host.jacobgeo_nodal->ncomp, up_int_host.jacobgeo_nodal->size);
  up_int_dev->g_ij_nodal = gkyl_array_cu_dev_new(up_int_host.g_ij_nodal->type,
    up_int_host.g_ij_nodal->ncomp, up_int_host.g_ij_nodal->size);
  up_int_dev->g_ij_neut_nodal = gkyl_array_cu_dev_new(up_int_host.g_ij_neut_nodal->type,
    up_int_host.g_ij_neut_nodal->ncomp, up_int_host.g_ij_neut_nodal->size);
  up_int_dev->dxdz_nodal = gkyl_array_cu_dev_new(up_int_host.dxdz_nodal->type,
    up_int_host.dxdz_nodal->ncomp, up_int_host.dxdz_nodal->size);
  up_int_dev->dzdx_nodal = gkyl_array_cu_dev_new(up_int_host.dzdx_nodal->type,
    up_int_host.dzdx_nodal->ncomp, up_int_host.dzdx_nodal->size);
  up_int_dev->dualmag_nodal = gkyl_array_cu_dev_new(up_int_host.dualmag_nodal->type,
    up_int_host.dualmag_nodal->ncomp, up_int_host.dualmag_nodal->size);
  up_int_dev->normals_nodal = gkyl_array_cu_dev_new(up_int_host.normals_nodal->type,
    up_int_host.normals_nodal->ncomp, up_int_host.normals_nodal->size);
  up_int_dev->gij_neut_nodal = gkyl_array_cu_dev_new(up_int_host.gij_neut_nodal->type,
    up_int_host.gij_neut_nodal->ncomp, up_int_host.gij_neut_nodal->size);
  up_int_dev->b_i_nodal = gkyl_array_cu_dev_new(up_int_host.b_i_nodal->type,
    up_int_host.b_i_nodal->ncomp, up_int_host.b_i_nodal->size);
  up_int_dev->b_i_nodal_fd = gkyl_array_cu_dev_new(up_int_host.b_i_nodal_fd->type,
    up_int_host.b_i_nodal_fd->ncomp, up_int_host.b_i_nodal_fd->size);
  up_int_dev->bcart_nodal = gkyl_array_cu_dev_new(up_int_host.bcart_nodal->type,
    up_int_host.bcart_nodal->ncomp, up_int_host.bcart_nodal->size);
  up_int_dev->B3_nodal = gkyl_array_cu_dev_new(up_int_host.B3_nodal->type,
    up_int_host.B3_nodal->ncomp, up_int_host.B3_nodal->size);
  up_int_dev->dualcurlbhatoverB_nodal = gkyl_array_cu_dev_new(up_int_host.dualcurlbhatoverB_nodal->type,
    up_int_host.dualcurlbhatoverB_nodal->ncomp, up_int_host.dualcurlbhatoverB_nodal->size);
  up_int_dev->rtg33inv_nodal = gkyl_array_cu_dev_new(up_int_host.rtg33inv_nodal->type,
    up_int_host.rtg33inv_nodal->ncomp, up_int_host.rtg33inv_nodal->size);
  up_int_dev->bioverJB_nodal = gkyl_array_cu_dev_new(up_int_host.bioverJB_nodal->type,
    up_int_host.bioverJB_nodal->ncomp, up_int_host.bioverJB_nodal->size);
  return up_int_dev;
}

struct gk_geom_surf*
gk_geometry_surf_cu_dev_alloc(struct gk_geom_surf up_surf_host)
{
  struct gk_geom_surf *up_surf_dev = (struct gk_geom_surf*) gkyl_malloc(sizeof(struct gk_geom_surf));
  // Expansions.
  up_surf_dev->bmag = gkyl_array_cu_dev_new(up_surf_host.bmag->type,
    up_surf_host.bmag->ncomp, up_surf_host.bmag->size);
  up_surf_dev->jacobgeo = gkyl_array_cu_dev_new(up_surf_host.jacobgeo->type,
    up_surf_host.jacobgeo->ncomp, up_surf_host.jacobgeo->size);
  up_surf_dev->jacobgeo_ratio = gkyl_array_cu_dev_new(up_surf_host.jacobgeo_ratio->type,
    up_surf_host.jacobgeo_ratio->ncomp, up_surf_host.jacobgeo_ratio->size);
  up_surf_dev->b_i = gkyl_array_cu_dev_new(up_surf_host.b_i->type,
    up_surf_host.b_i->ncomp, up_surf_host.b_i->size);
  up_surf_dev->cmag = gkyl_array_cu_dev_new(up_surf_host.cmag->type,
    up_surf_host.cmag->ncomp, up_surf_host.cmag->size);
  up_surf_dev->jacobtot_inv = gkyl_array_cu_dev_new(up_surf_host.jacobtot_inv->type,
    up_surf_host.jacobtot_inv->ncomp, up_surf_host.jacobtot_inv->size);
  up_surf_dev->B3 = gkyl_array_cu_dev_new(up_surf_host.B3->type,
    up_surf_host.B3->ncomp, up_surf_host.B3->size);
  up_surf_dev->normcurlbhat = gkyl_array_cu_dev_new(up_surf_host.normcurlbhat->type,
    up_surf_host.normcurlbhat->ncomp, up_surf_host.normcurlbhat->size);
  up_surf_dev->normals = gkyl_array_cu_dev_new(up_surf_host.normals->type,
    up_surf_host.normals->ncomp, up_surf_host.normals->size);
  up_surf_dev->lenr = gkyl_array_cu_dev_new(up_surf_host.lenr->type,
    up_surf_host.lenr->ncomp, up_surf_host.lenr->size);
  up_surf_dev->bimpactangle = gkyl_array_cu_dev_new(up_surf_host.bimpactangle->type,
    up_surf_host.bimpactangle->ncomp, up_surf_host.bimpactangle->size);
  up_surf_dev->deltats = gkyl_array_cu_dev_new(up_surf_host.deltats->type,
    up_surf_host.deltats->ncomp, up_surf_host.deltats->size);
  // Nodal.
  up_surf_dev->mc2p_nodal_fd = gkyl_array_cu_dev_new(up_surf_host.mc2p_nodal_fd->type,
    up_surf_host.mc2p_nodal_fd->ncomp, up_surf_host.mc2p_nodal_fd->size);
  up_surf_dev->mc2p_nodal = gkyl_array_cu_dev_new(up_surf_host.mc2p_nodal->type,
    up_surf_host.mc2p_nodal->ncomp, up_surf_host.mc2p_nodal->size);
  up_surf_dev->bmag_nodal = gkyl_array_cu_dev_new(up_surf_host.bmag_nodal->type,
    up_surf_host.bmag_nodal->ncomp, up_surf_host.bmag_nodal->size);
  up_surf_dev->curlbhat_nodal = gkyl_array_cu_dev_new(up_surf_host.curlbhat_nodal->type,
    up_surf_host.curlbhat_nodal->ncomp, up_surf_host.curlbhat_nodal->size);
  up_surf_dev->normcurlbhat_nodal = gkyl_array_cu_dev_new(up_surf_host.normcurlbhat_nodal->type,
    up_surf_host.normcurlbhat_nodal->ncomp, up_surf_host.normcurlbhat_nodal->size);
  up_surf_dev->ddtheta_nodal = gkyl_array_cu_dev_new(up_surf_host.ddtheta_nodal->type,
    up_surf_host.ddtheta_nodal->ncomp, up_surf_host.ddtheta_nodal->size);
  up_surf_dev->ddpsi_nodal = gkyl_array_cu_dev_new(up_surf_host.ddpsi_nodal->type,
    up_surf_host.ddpsi_nodal->ncomp, up_surf_host.ddpsi_nodal->size);
  up_surf_dev->jacobgeo_nodal = gkyl_array_cu_dev_new(up_surf_host.jacobgeo_nodal->type,
    up_surf_host.jacobgeo_nodal->ncomp, up_surf_host.jacobgeo_nodal->size);
  up_surf_dev->b_i_nodal = gkyl_array_cu_dev_new(up_surf_host.b_i_nodal->type,
    up_surf_host.b_i_nodal->ncomp, up_surf_host.b_i_nodal->size);
  up_surf_dev->b_i_nodal_fd = gkyl_array_cu_dev_new(up_surf_host.b_i_nodal_fd->type,
    up_surf_host.b_i_nodal_fd->ncomp, up_surf_host.b_i_nodal_fd->size);
  up_surf_dev->cmag_nodal = gkyl_array_cu_dev_new(up_surf_host.cmag_nodal->type,
    up_surf_host.cmag_nodal->ncomp, up_surf_host.cmag_nodal->size);
  up_surf_dev->jacobtot_inv_nodal = gkyl_array_cu_dev_new(up_surf_host.jacobtot_inv_nodal->type,
    up_surf_host.jacobtot_inv_nodal->ncomp, up_surf_host.jacobtot_inv_nodal->size);
  up_surf_dev->g_ij_nodal = gkyl_array_cu_dev_new(up_surf_host.g_ij_nodal->type,
    up_surf_host.g_ij_nodal->ncomp, up_surf_host.g_ij_nodal->size);
  up_surf_dev->dxdz_nodal = gkyl_array_cu_dev_new(up_surf_host.dxdz_nodal->type,
    up_surf_host.dxdz_nodal->ncomp, up_surf_host.dxdz_nodal->size);
  up_surf_dev->dzdx_nodal = gkyl_array_cu_dev_new(up_surf_host.dzdx_nodal->type,
    up_surf_host.dzdx_nodal->ncomp, up_surf_host.dzdx_nodal->size);
  up_surf_dev->normals_nodal = gkyl_array_cu_dev_new(up_surf_host.normals_nodal->type,
    up_surf_host.normals_nodal->ncomp, up_surf_host.normals_nodal->size);
  up_surf_dev->dualmag_nodal = gkyl_array_cu_dev_new(up_surf_host.dualmag_nodal->type,
    up_surf_host.dualmag_nodal->ncomp, up_surf_host.dualmag_nodal->size);
  up_surf_dev->bcart_nodal = gkyl_array_cu_dev_new(up_surf_host.bcart_nodal->type,
    up_surf_host.bcart_nodal->ncomp, up_surf_host.bcart_nodal->size);
  up_surf_dev->B3_nodal = gkyl_array_cu_dev_new(up_surf_host.B3_nodal->type,
    up_surf_host.B3_nodal->ncomp, up_surf_host.B3_nodal->size);
  up_surf_dev->lenr_nodal = gkyl_array_cu_dev_new(up_surf_host.lenr_nodal->type,
    up_surf_host.lenr_nodal->ncomp, up_surf_host.lenr_nodal->size);
  up_surf_dev->bimpactangle_nodal = gkyl_array_cu_dev_new(up_surf_host.bimpactangle_nodal->type,
    up_surf_host.bimpactangle_nodal->ncomp, up_surf_host.bimpactangle_nodal->size);
  up_surf_dev->deltats_nodal = gkyl_array_cu_dev_new(up_surf_host.deltats_nodal->type,
    up_surf_host.deltats_nodal->ncomp, up_surf_host.deltats_nodal->size);
  return up_surf_dev;
}

// CPU interface to create and track a GPU object
struct gk_geometry* 
gkyl_gk_geometry_cu_dev_new(struct gk_geometry* geo_host, struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gk_geometry *up = (struct gk_geometry*) gkyl_malloc(sizeof(struct gk_geometry));

  up->geometry_id = geometry_inp->geometry_id;
  up->basis = geometry_inp->basis;
  up->local = geometry_inp->local;
  up->local_ext = geometry_inp->local_ext;
  up->global = geometry_inp->global;
  up->global_ext = geometry_inp->global_ext;
  up->grid = geometry_inp->grid;
  gk_geometry_set_nodal_ranges(up) ;
  if (up->grid.ndim > 1) {
    gkyl_cart_modal_serendip(&up->surf_basis, up->grid.ndim-1, up->basis.poly_order);
    up->num_surf_basis = up->surf_basis.num_basis;
  }
  else {
    up->num_surf_basis = 1;
  }
  up->geqdsk_sign_convention = geo_host->geqdsk_sign_convention;
  up->half_domain = geo_host->half_domain;
  up->has_LCFS = geo_host->has_LCFS;
  if (up->has_LCFS) {
    up->x_LCFS = geo_host->x_LCFS;
    // Check that the split happens within the domain.
    assert((up->grid.lower[0] <= up->x_LCFS) && (up->x_LCFS <= up->grid.upper[0]));
    // Check that the split happens at a cell boundary;
    double needint = (up->x_LCFS - up->grid.lower[0])/up->grid.dx[0];
    double rem_floor = fabs(needint-floor(needint));
    double rem_ceil = fabs(needint-ceil(needint));
    if (rem_floor < 1.0e-12) {
      up->idx_LCFS_lo = (int) floor(needint);
    }
    else if (rem_ceil < 1.0e-12) {
      up->idx_LCFS_lo = (int) ceil(needint);
    }
    else {
      fprintf(stderr, "x_LCFS = %.9e must be at a cell boundary.\n", up->x_LCFS);
      assert(false);
    }
  }

  // Function pointers to twistshift function.
  up->parallel_lower_bc_shift_func = geometry_inp->parallel_lower_bc_shift_func;
  up->parallel_upper_bc_shift_func = geometry_inp->parallel_upper_bc_shift_func;
  up->parallel_lower_bc_shift_ctx  = geometry_inp->parallel_lower_bc_shift_ctx ;
  up->parallel_upper_bc_shift_ctx  = geometry_inp->parallel_upper_bc_shift_ctx ;

  struct gk_geom_corn *geo_corn_dev = gk_geometry_corn_cu_dev_alloc(geo_host->geo_corn);
  struct gk_geom_int *geo_int_dev = gk_geometry_int_cu_dev_alloc(geo_host->geo_int);
  struct gk_geom_surf *geo_surf_dev[up->grid.ndim];
  for (int dir=0; dir<up->grid.ndim; ++dir)
    geo_surf_dev[dir] = gk_geometry_surf_cu_dev_alloc(geo_host->geo_surf[dir]);

  // Expansions.
  gkyl_array_copy(geo_corn_dev->mc2p, geo_host->geo_corn.mc2p);
  gkyl_array_copy(geo_corn_dev->mc2p_deflated, geo_host->geo_corn.mc2p_deflated);
  gkyl_array_copy(geo_corn_dev->mc2nu_pos, geo_host->geo_corn.mc2nu_pos);
  gkyl_array_copy(geo_corn_dev->mc2nu_pos_deflated, geo_host->geo_corn.mc2nu_pos_deflated);
  gkyl_array_copy(geo_corn_dev->bmag, geo_host->geo_corn.bmag);
  gkyl_array_copy(geo_corn_dev->bmag_inv, geo_host->geo_corn.bmag_inv);
  // Nodal.
  gkyl_array_copy(geo_corn_dev->mc2p_nodal, geo_host->geo_corn.mc2p_nodal);
  gkyl_array_copy(geo_corn_dev->mc2nu_pos_nodal, geo_host->geo_corn.mc2nu_pos_nodal);
  gkyl_array_copy(geo_corn_dev->bmag_nodal, geo_host->geo_corn.bmag_nodal);

  // Expansions.
  gkyl_array_copy(geo_int_dev->mc2p, geo_host->geo_int.mc2p);
  gkyl_array_copy(geo_int_dev->bmag, geo_host->geo_int.bmag);
  gkyl_array_copy(geo_int_dev->g_ij, geo_host->geo_int.g_ij);
  gkyl_array_copy(geo_int_dev->g_ij_neut, geo_host->geo_int.g_ij_neut);
  gkyl_array_copy(geo_int_dev->dxdz, geo_host->geo_int.dxdz);
  gkyl_array_copy(geo_int_dev->dzdx, geo_host->geo_int.dzdx);
  gkyl_array_copy(geo_int_dev->dualmag, geo_host->geo_int.dualmag);
  gkyl_array_copy(geo_int_dev->normals, geo_host->geo_int.normals);
  gkyl_array_copy(geo_int_dev->jacobgeo, geo_host->geo_int.jacobgeo);
  gkyl_array_copy(geo_int_dev->jacobgeo_ghost, geo_host->geo_int.jacobgeo_ghost);
  gkyl_array_copy(geo_int_dev->jacobgeo_inv, geo_host->geo_int.jacobgeo_inv);
  gkyl_array_copy(geo_int_dev->gij, geo_host->geo_int.gij);
  gkyl_array_copy(geo_int_dev->gij_neut, geo_host->geo_int.gij_neut);
  gkyl_array_copy(geo_int_dev->b_i, geo_host->geo_int.b_i);
  gkyl_array_copy(geo_int_dev->bcart, geo_host->geo_int.bcart);
  gkyl_array_copy(geo_int_dev->cmag, geo_host->geo_int.cmag);
  gkyl_array_copy(geo_int_dev->jacobtot, geo_host->geo_int.jacobtot);
  gkyl_array_copy(geo_int_dev->jacobtot_inv, geo_host->geo_int.jacobtot_inv);
  gkyl_array_copy(geo_int_dev->gxxj, geo_host->geo_int.gxxj);
  gkyl_array_copy(geo_int_dev->gxyj, geo_host->geo_int.gxyj);
  gkyl_array_copy(geo_int_dev->gyyj, geo_host->geo_int.gyyj);
  gkyl_array_copy(geo_int_dev->gxzj, geo_host->geo_int.gxzj);
  gkyl_array_copy(geo_int_dev->eps2, geo_host->geo_int.eps2);
  gkyl_array_copy(geo_int_dev->rtg33inv, geo_host->geo_int.rtg33inv);
  gkyl_array_copy(geo_int_dev->dualcurlbhatoverB, geo_host->geo_int.dualcurlbhatoverB);
  gkyl_array_copy(geo_int_dev->bioverJB, geo_host->geo_int.bioverJB);
  gkyl_array_copy(geo_int_dev->B3, geo_host->geo_int.B3);
  gkyl_array_copy(geo_int_dev->dualcurlbhat, geo_host->geo_int.dualcurlbhat);
  gkyl_array_copy(geo_int_dev->qprofile, geo_host->geo_int.qprofile);
  // Nodal.
  gkyl_array_copy(geo_int_dev->mc2p_nodal_fd, geo_host->geo_int.mc2p_nodal_fd); 
  gkyl_array_copy(geo_int_dev->mc2p_nodal, geo_host->geo_int.mc2p_nodal); 
  gkyl_array_copy(geo_int_dev->bmag_nodal, geo_host->geo_int.bmag_nodal); 
  gkyl_array_copy(geo_int_dev->ddtheta_nodal, geo_host->geo_int.ddtheta_nodal); 
  gkyl_array_copy(geo_int_dev->ddpsi_nodal, geo_host->geo_int.ddpsi_nodal); 
  gkyl_array_copy(geo_int_dev->curlbhat_nodal, geo_host->geo_int.curlbhat_nodal); 
  gkyl_array_copy(geo_int_dev->dualcurlbhat_nodal, geo_host->geo_int.dualcurlbhat_nodal); 
  gkyl_array_copy(geo_int_dev->jacobgeo_nodal, geo_host->geo_int.jacobgeo_nodal); 
  gkyl_array_copy(geo_int_dev->g_ij_nodal, geo_host->geo_int.g_ij_nodal); 
  gkyl_array_copy(geo_int_dev->g_ij_neut_nodal, geo_host->geo_int.g_ij_neut_nodal); 
  gkyl_array_copy(geo_int_dev->dxdz_nodal, geo_host->geo_int.dxdz_nodal); 
  gkyl_array_copy(geo_int_dev->dzdx_nodal, geo_host->geo_int.dzdx_nodal); 
  gkyl_array_copy(geo_int_dev->dualmag_nodal, geo_host->geo_int.dualmag_nodal); 
  gkyl_array_copy(geo_int_dev->normals_nodal, geo_host->geo_int.normals_nodal); 
  gkyl_array_copy(geo_int_dev->gij_neut_nodal, geo_host->geo_int.gij_neut_nodal); 
  gkyl_array_copy(geo_int_dev->b_i_nodal, geo_host->geo_int.b_i_nodal); 
  gkyl_array_copy(geo_int_dev->b_i_nodal_fd, geo_host->geo_int.b_i_nodal_fd); 
  gkyl_array_copy(geo_int_dev->bcart_nodal, geo_host->geo_int.bcart_nodal); 
  gkyl_array_copy(geo_int_dev->B3_nodal, geo_host->geo_int.B3_nodal); 
  gkyl_array_copy(geo_int_dev->dualcurlbhatoverB_nodal, geo_host->geo_int.dualcurlbhatoverB_nodal); 
  gkyl_array_copy(geo_int_dev->rtg33inv_nodal, geo_host->geo_int.rtg33inv_nodal); 
  gkyl_array_copy(geo_int_dev->bioverJB_nodal, geo_host->geo_int.bioverJB_nodal); 

  for (int dir=0; dir<up->grid.ndim; ++dir) {
    // Expansions.
    gkyl_array_copy(geo_surf_dev[dir]->bmag, geo_host->geo_surf[dir].bmag);
    gkyl_array_copy(geo_surf_dev[dir]->jacobgeo, geo_host->geo_surf[dir].jacobgeo);
    gkyl_array_copy(geo_surf_dev[dir]->jacobgeo_ratio, geo_host->geo_surf[dir].jacobgeo_ratio);
    gkyl_array_copy(geo_surf_dev[dir]->b_i, geo_host->geo_surf[dir].b_i);
    gkyl_array_copy(geo_surf_dev[dir]->cmag, geo_host->geo_surf[dir].cmag);
    gkyl_array_copy(geo_surf_dev[dir]->jacobtot_inv, geo_host->geo_surf[dir].jacobtot_inv);
    gkyl_array_copy(geo_surf_dev[dir]->B3, geo_host->geo_surf[dir].B3);
    gkyl_array_copy(geo_surf_dev[dir]->normcurlbhat, geo_host->geo_surf[dir].normcurlbhat);
    gkyl_array_copy(geo_surf_dev[dir]->normals, geo_host->geo_surf[dir].normals);
    gkyl_array_copy(geo_surf_dev[dir]->lenr, geo_host->geo_surf[dir].lenr);
    gkyl_array_copy(geo_surf_dev[dir]->bimpactangle, geo_host->geo_surf[dir].bimpactangle);
    gkyl_array_copy(geo_surf_dev[dir]->deltats, geo_host->geo_surf[dir].deltats);
    // Nodal.
    gkyl_array_copy(geo_surf_dev[dir]->mc2p_nodal_fd, geo_host->geo_surf[dir].mc2p_nodal_fd);
    gkyl_array_copy(geo_surf_dev[dir]->mc2p_nodal, geo_host->geo_surf[dir].mc2p_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->bmag_nodal, geo_host->geo_surf[dir].bmag_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->curlbhat_nodal, geo_host->geo_surf[dir].curlbhat_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->normcurlbhat_nodal, geo_host->geo_surf[dir].normcurlbhat_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->ddtheta_nodal, geo_host->geo_surf[dir].ddtheta_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->ddpsi_nodal, geo_host->geo_surf[dir].ddpsi_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->jacobgeo_nodal, geo_host->geo_surf[dir].jacobgeo_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->b_i_nodal, geo_host->geo_surf[dir].b_i_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->b_i_nodal_fd, geo_host->geo_surf[dir].b_i_nodal_fd);
    gkyl_array_copy(geo_surf_dev[dir]->cmag_nodal, geo_host->geo_surf[dir].cmag_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->jacobtot_inv_nodal, geo_host->geo_surf[dir].jacobtot_inv_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->g_ij_nodal, geo_host->geo_surf[dir].g_ij_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->dxdz_nodal, geo_host->geo_surf[dir].dxdz_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->dzdx_nodal, geo_host->geo_surf[dir].dzdx_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->normals_nodal, geo_host->geo_surf[dir].normals_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->dualmag_nodal, geo_host->geo_surf[dir].dualmag_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->bcart_nodal, geo_host->geo_surf[dir].bcart_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->B3_nodal, geo_host->geo_surf[dir].B3_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->lenr_nodal, geo_host->geo_surf[dir].lenr_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->bimpactangle_nodal, geo_host->geo_surf[dir].bimpactangle_nodal);
    gkyl_array_copy(geo_surf_dev[dir]->deltats_nodal, geo_host->geo_surf[dir].deltats_nodal);
  }

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);

  // Initialize the device geometry object
  struct gk_geometry *up_cu = (struct gk_geometry*) gkyl_cu_malloc(sizeof(struct gk_geometry));
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gk_geometry), GKYL_CU_MEMCPY_H2D);
  gkyl_geometry_set_corn_cu(up_cu, geo_corn_dev);
  gkyl_geometry_set_int_cu(up_cu, geo_int_dev);
  for (int dir=0; dir<up->grid.ndim; ++dir)
    gkyl_geometry_set_surf_cu(up_cu, geo_surf_dev[dir], dir);

  up->on_dev = up_cu;

  // Geometry object should store host pointer.
  // Expansions.
  up->geo_corn.mc2p = geo_corn_dev->mc2p;
  up->geo_corn.mc2nu_pos = geo_corn_dev->mc2nu_pos;
  up->geo_corn.bmag = geo_corn_dev->bmag;
  up->geo_corn.bmag_inv = geo_corn_dev->bmag_inv;
  up->geo_corn.mc2p_deflated = geo_corn_dev->mc2p_deflated;
  up->geo_corn.mc2nu_pos_deflated = geo_corn_dev->mc2nu_pos_deflated;
  // Nodal.
  up->geo_corn.mc2p_nodal = geo_corn_dev->mc2p_nodal;
  up->geo_corn.mc2nu_pos_nodal = geo_corn_dev->mc2nu_pos_nodal;
  up->geo_corn.bmag_nodal = geo_corn_dev->bmag_nodal;
  gkyl_free(geo_corn_dev);

  // Expansions.
  up->geo_int.mc2p = geo_int_dev->mc2p;
  up->geo_int.bmag = geo_int_dev->bmag;
  up->geo_int.g_ij = geo_int_dev->g_ij;
  up->geo_int.g_ij_neut = geo_int_dev->g_ij_neut;
  up->geo_int.dxdz = geo_int_dev->dxdz;
  up->geo_int.dzdx = geo_int_dev->dzdx;
  up->geo_int.dualmag = geo_int_dev->dualmag;
  up->geo_int.normals = geo_int_dev->normals;
  up->geo_int.jacobgeo = geo_int_dev->jacobgeo;
  up->geo_int.jacobgeo_ghost = geo_int_dev->jacobgeo_ghost;
  up->geo_int.jacobgeo_inv = geo_int_dev->jacobgeo_inv;
  up->geo_int.gij = geo_int_dev->gij;
  up->geo_int.gij_neut = geo_int_dev->gij_neut;
  up->geo_int.b_i = geo_int_dev->b_i;
  up->geo_int.bcart = geo_int_dev->bcart;
  up->geo_int.cmag = geo_int_dev->cmag;
  up->geo_int.jacobtot = geo_int_dev->jacobtot;
  up->geo_int.jacobtot_inv = geo_int_dev->jacobtot_inv;
  up->geo_int.gxxj = geo_int_dev->gxxj;
  up->geo_int.gxyj = geo_int_dev->gxyj;
  up->geo_int.gyyj = geo_int_dev->gyyj;
  up->geo_int.gxzj = geo_int_dev->gxzj;
  up->geo_int.eps2 = geo_int_dev->eps2;
  up->geo_int.rtg33inv = geo_int_dev->rtg33inv;
  up->geo_int.dualcurlbhatoverB = geo_int_dev->dualcurlbhatoverB;
  up->geo_int.bioverJB = geo_int_dev->bioverJB;
  up->geo_int.B3 = geo_int_dev->B3;
  up->geo_int.dualcurlbhat = geo_int_dev->dualcurlbhat;
  up->geo_int.qprofile = geo_int_dev->qprofile;
  // Nodal.
  up->geo_int.mc2p_nodal_fd = geo_int_dev->mc2p_nodal_fd;
  up->geo_int.mc2p_nodal = geo_int_dev->mc2p_nodal;
  up->geo_int.bmag_nodal = geo_int_dev->bmag_nodal;
  up->geo_int.ddtheta_nodal = geo_int_dev->ddtheta_nodal;
  up->geo_int.ddpsi_nodal = geo_int_dev->ddpsi_nodal;
  up->geo_int.curlbhat_nodal = geo_int_dev->curlbhat_nodal;
  up->geo_int.dualcurlbhat_nodal = geo_int_dev->dualcurlbhat_nodal;
  up->geo_int.jacobgeo_nodal = geo_int_dev->jacobgeo_nodal;
  up->geo_int.g_ij_nodal = geo_int_dev->g_ij_nodal;
  up->geo_int.g_ij_neut_nodal = geo_int_dev->g_ij_neut_nodal;
  up->geo_int.dxdz_nodal = geo_int_dev->dxdz_nodal;
  up->geo_int.dzdx_nodal = geo_int_dev->dzdx_nodal;
  up->geo_int.dualmag_nodal = geo_int_dev->dualmag_nodal;
  up->geo_int.normals_nodal = geo_int_dev->normals_nodal;
  up->geo_int.gij_neut_nodal = geo_int_dev->gij_neut_nodal;
  up->geo_int.b_i_nodal = geo_int_dev->b_i_nodal;
  up->geo_int.b_i_nodal_fd = geo_int_dev->b_i_nodal_fd;
  up->geo_int.bcart_nodal = geo_int_dev->bcart_nodal;
  up->geo_int.B3_nodal = geo_int_dev->B3_nodal;
  up->geo_int.dualcurlbhatoverB_nodal = geo_int_dev->dualcurlbhatoverB_nodal;
  up->geo_int.rtg33inv_nodal = geo_int_dev->rtg33inv_nodal;
  up->geo_int.bioverJB_nodal = geo_int_dev->bioverJB_nodal;
  gkyl_free(geo_int_dev);

  for (int dir=0; dir<up->grid.ndim; ++dir) {
    // Expansions.
    up->geo_surf[dir].bmag = geo_surf_dev[dir]->bmag;
    up->geo_surf[dir].jacobgeo = geo_surf_dev[dir]->jacobgeo;
    up->geo_surf[dir].jacobgeo_ratio = geo_surf_dev[dir]->jacobgeo_ratio;
    up->geo_surf[dir].b_i = geo_surf_dev[dir]->b_i;
    up->geo_surf[dir].cmag = geo_surf_dev[dir]->cmag;
    up->geo_surf[dir].jacobtot_inv = geo_surf_dev[dir]->jacobtot_inv;
    up->geo_surf[dir].B3 = geo_surf_dev[dir]->B3;
    up->geo_surf[dir].normcurlbhat = geo_surf_dev[dir]->normcurlbhat;
    up->geo_surf[dir].normals = geo_surf_dev[dir]->normals;
    up->geo_surf[dir].lenr = geo_surf_dev[dir]->lenr;
    up->geo_surf[dir].bimpactangle = geo_surf_dev[dir]->bimpactangle;
    up->geo_surf[dir].deltats = geo_surf_dev[dir]->deltats;
    // Nodal.
    up->geo_surf[dir].mc2p_nodal_fd = geo_surf_dev[dir]->mc2p_nodal_fd;
    up->geo_surf[dir].mc2p_nodal = geo_surf_dev[dir]->mc2p_nodal;
    up->geo_surf[dir].bmag_nodal = geo_surf_dev[dir]->bmag_nodal;
    up->geo_surf[dir].curlbhat_nodal = geo_surf_dev[dir]->curlbhat_nodal;
    up->geo_surf[dir].normcurlbhat_nodal = geo_surf_dev[dir]->normcurlbhat_nodal;
    up->geo_surf[dir].ddtheta_nodal = geo_surf_dev[dir]->ddtheta_nodal;
    up->geo_surf[dir].ddpsi_nodal = geo_surf_dev[dir]->ddpsi_nodal;
    up->geo_surf[dir].jacobgeo_nodal = geo_surf_dev[dir]->jacobgeo_nodal;
    up->geo_surf[dir].b_i_nodal = geo_surf_dev[dir]->b_i_nodal;
    up->geo_surf[dir].b_i_nodal_fd = geo_surf_dev[dir]->b_i_nodal_fd;
    up->geo_surf[dir].cmag_nodal = geo_surf_dev[dir]->cmag_nodal;
    up->geo_surf[dir].jacobtot_inv_nodal = geo_surf_dev[dir]->jacobtot_inv_nodal;
    up->geo_surf[dir].g_ij_nodal = geo_surf_dev[dir]->g_ij_nodal;
    up->geo_surf[dir].dxdz_nodal = geo_surf_dev[dir]->dxdz_nodal;
    up->geo_surf[dir].dzdx_nodal = geo_surf_dev[dir]->dzdx_nodal;
    up->geo_surf[dir].normals_nodal = geo_surf_dev[dir]->normals_nodal;
    up->geo_surf[dir].dualmag_nodal = geo_surf_dev[dir]->dualmag_nodal;
    up->geo_surf[dir].bcart_nodal = geo_surf_dev[dir]->bcart_nodal;
    up->geo_surf[dir].B3_nodal = geo_surf_dev[dir]->B3_nodal;
    up->geo_surf[dir].lenr_nodal = geo_surf_dev[dir]->lenr_nodal;
    up->geo_surf[dir].bimpactangle_nodal = geo_surf_dev[dir]->bimpactangle_nodal;
    up->geo_surf[dir].deltats_nodal = geo_surf_dev[dir]->deltats_nodal;
    gkyl_free(geo_surf_dev[dir]);
  }

  // Copy metadata.
  up->io_meta_basic = gkyl_msgpack_map_elem_clone(geo_host->io_meta_basic_len, geo_host->io_meta_basic);
  up->io_meta_basic_len = geo_host->io_meta_basic_len;
 
  return up;
}

