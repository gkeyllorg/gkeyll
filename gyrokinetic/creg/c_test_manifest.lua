-- Tests skipped by the C regression suite.
-- Remove an entry manually to re-enable the test.
-- gpu: C tests whose GPU variant timed out (CPU variant still runs).
return {
   ignore = {
      tests = {
         "rt_gk_ltx_iwl_num_miller_3x2v_p1",
         "rt_gk_multib_asdex_solonly_3x2v_p1",
         "rt_gk_multib_nstx_solonly_3x2v_p1",
         "rt_gk_multib_tcv_x21_3x2v_p1",
         "rt_gk_sheath_3x2v_p1_cons",
      },
      gpu = {
         "rt_gk_ltx_iwl_num_miller_3x2v_p1",
         "rt_gk_multib_asdex_solonly_3x2v_p1",
         "rt_gk_multib_nstx_solonly_3x2v_p1",
         "rt_gk_multib_sheath_1x2v_p1",
         "rt_gk_multib_step_sol_1x2v_p1",
         "rt_gk_multib_tcv_x21_3x2v_p1",
         "rt_gk_sheath_3x2v_p1_cons",
      },
   },
   moat = {
      "rt_gk_ion_sound_1x2v_p1", "rt_gk_lbo_relax_1x2v_p1",
      "rt_gk_sheath_bgk_1x2v_p1",
   },
   parallel = {
      { name = "rt_gk_cbc_3x2v_p1", cuts = { 1, 1, 4 } },
      { name = "rt_gk_d3d_iwl_2x2v_p1", cuts = { 1, 4 } },
      { name = "rt_gk_ion_sound_adiabatic_elc_1x2v_p1", cuts = { 4 } },
      { name = "rt_gk_lapd_cyl_3x2v_p1", cuts = { 1, 1, 4 } },
      { name = "rt_gk_ltx_1x2v_p1", cuts = { 4 } },
      { name = "rt_gk_mirror_boltz_elc_1x2v_p1", cuts = { 4 } },
      { name = "rt_gk_multib_step_sol_2x2v_p1", cuts = { 1, 3 } },
      { name = "rt_gk_sheath_2x2v_p1", cuts = { 1, 4 } },
      { name = "rt_gk_tcv_iwl_adapt_source_3x2v_p1", cuts = { 1, 1, 4 } },
      { name = "rt_gk_wham_boltz_elc_poa_1x2v_p1", cuts = { 4 } },
   },
}
