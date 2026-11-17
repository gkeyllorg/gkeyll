#pragma once
// Per-surface metadata for the generated nodal surface-evaluation tables
// (gkyl_vlasov_flux_surf_mod2nod_tables_*.h) and their per-node kernels. One instance per
// surface (vst_<tag>_ph_v<d>_meta / vst_<tag>_ph_x<d>_meta) is emitted next
// to the tables. Loop bounds and array extents are compiled into the node
// functions as literals; this record is for the host-side launch logic
// (block sizing, shared-memory sizing, output offsets).
struct gkyl_vlasov_surf_meta {
  int nb; // phase-basis coefficients (f/hamil length; stage-1 loop bound)
  int na; // distinct outer (chi) shapes: Cm columns, G-column length
  int nmv; // distinct inner (psi) shapes: shared V columns
  int no; // outer nodes (Cm rows; i extent)
  int ni; // inner nodes (V rows; j extent; flat node n = i*ni + j)
  int ns; // sparse-Hamiltonian support size (vel_sparse_idx length)
  int fhat_off; // base offset of this surface's block in Fhat_nodal/alpha_quad
};
