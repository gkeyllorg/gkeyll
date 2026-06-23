-- Tests skipped by the C regression suite.
-- Remove an entry manually to re-enable the test.
-- gpu: C tests whose GPU variant timed out (CPU variant still runs).
return {
   tests = {
      "rt_can_pb_bgk_surf_cylindrical_sodshock_im_1x3v_p2",
      "rt_can_pb_bgk_surf_annulus_sodshock_im_2x2v_p2",
      "rt_can_pb_bgk_surf_sphere_sodshock_im_2x2v_p1",
      "rt_vlasov_neut_bgk_sodshock_1x3v_p2",
      "rt_can_pb_bgk_surf_sphere_sodshock_im_1x2v_p2",
      "rt_can_pb_ex_bgk_surf_flat_sq_ic",
      "rt_can_pb_bgk_surf_sphere_sodshock_im_2x2v_p2",
      "rt_escreen_sr",
      "rt_hyper_vlasov_tm",
      "rt_vlasov_kerntm",
      "rt_vlasov_moments",
      "rt_can_pb_bgk_surf_cylindrical_sodshock_im_3x3v_p1",
      "rt_can_pb_bgk_surf_toroidal_sodshock_im_3x3v_p1",
      "rt_vlasov_lbo_cross_1x2v_p2",
      "rt_vp_sheath_Aext_1x2v_p2",
      "rt_vlasov_neut_lbo_sodshock_1x3v_p2",
      "rt_can_pb_neut_bgk_sodshock_im_1x2v_p2",
      "rt_vlasov_sheath_bgk_1x1v_p2",
      "rt_can_pb_bgk_surf_annulus_sodshock_im_2x2v_p1",
      "rt_vlasov_weibel_2x2v_p2",
      "rt_vlasov_es_shock_lbo_1x3v",
      "rt_can_pb_bgk_surf_toroidal_sodshock_im_2x3v_p2",
      "rt_dg_euler_kh_2d",
      "rt_can_pb_bgk_surf_annulus_sodshock_im_1x2v_p2",
      "rt_vp_sheath_Bext_1x2v_p2",
   },
   gpu = {
   },
}
