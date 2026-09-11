-- Tests skipped by the Lua regression suite.
-- Remove an entry manually to re-enable the test.
-- gpu: Lua tests whose GPU variant timed out (CPU variant still runs).
return {
   tests = {
      "rt_gk_sheath_1x2v_p1",
      "rt_gk_sheath_1x2v_p1_cons",
      "rt_gk_sheath_2x2v_p1",
      "rt_gk_sheath_3x2v_p1",
      "rt_gk_sheath_bgk_1x2v_p1",
      "rt_gk_sheath_nonuniformx_1x2v_p1",
      "rt_gk_sheath_nonuniformx_2x2v_p1",
      "rt_gk_sheath_nonuniformx_3x2v_p1",
      "rt_gk_sheath_nonuniformv_1x2v_p1",
      "rt_gk_sheath_nonuniformv_2x2v_p1",
      "rt_gk_sheath_nonuniformv_3x2v_p1",
      "rt_gk_bgk_periodic_sodshock_1x2v_p1",
      "rt_gk_bgk_relax_1x2v_p1",
      "rt_gk_bgk_cross_relax_1x2v_p1",
      "rt_gk_lbo_cross_relax_1x2v_p1",
      "rt_gk_ar_react_nonuniformv_1x2v_p1",
      "rt_gk_rad_1x2v_p1",
      "rt_gk_rad_nonuniformv_1x2v_p1",
      "rt_gk_li_react_nonuniformv_3x2v_p1",
      "rt_gk_li_react_3x2v_p1",
      "rt_gk_lbo_relax_bimaxwellian_nonuniformv_3x2v_p1",
      "rt_gk_bgk_periodic_sodshock_1x2v_p1",  -- crashes compareFiles (C assertion in gkyl_aligned_free_)
   },
   gpu = {
   },
}
