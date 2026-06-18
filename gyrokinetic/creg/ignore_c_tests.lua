-- Tests skipped by the C regression suite.
-- Remove an entry manually to re-enable the test.
-- gpu: C tests whose GPU variant timed out (CPU variant still runs).
return {
   tests = {
      "rt_gk_mirror_boltz_elc_poa_1x2v_p1",
      "rt_gk_sheath_nonuniformx_3x2v_p1",
      "rt_gk_multib_step_sol_2x2v_p1",
      "rt_gk_wham_3x2v_p1",
      "rt_gk_multib_step_nonuniform_2x2v_p1",
      "rt_gk_wham_boltz_elc_poa_1x2v_p1",
      "rt_gk_multib_step_2x2v_p1",
      "rt_gk_wham_nonuniformx_1x2v_p1",
      "rt_gk_wham_1x2v_p1_static_field",
      "rt_gk_wham_nonuniformx_3x2v_p1",
      "rt_gk_wham_2xIC_3x2v_p1",
      "rt_gk_lbo_relax_bimaxwellian_nonuniformv_3x2v_p1",
      "rt_gk_leaky_bag_boltz_sheath_1x2v_p1_nux_edge",
      "rt_gk_sheath_3x2v_p1",
      "rt_gk_ltx_1x2v_p1",
      "rt_gk_multib_asdex_2x2v_p1",
      "rt_gk_bgk_periodic_sodshock_1x2v_p1",    -- crashes compareFiles (C assertion in gkyl_aligned_free_)
      "rt_gk_bgk_im_periodic_sodshock_1x2v_p1", -- crashes compareFiles (C assertion in gkyl_aligned_free_)
      "rt_gk_wham_1xIC_2x2v_p1",                -- crashes compareFiles (C assertion in gkyl_aligned_free_)
   },
   gpu = {
   },
}
