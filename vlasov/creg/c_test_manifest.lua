-- Tests skipped by the C regression suite.
-- Remove an entry manually to re-enable the test.
-- gpu: C tests whose GPU variant timed out (CPU variant still runs).
return {
   ignore = {
      tests = {
         "rt_triad_bgk_spherical_blast_1x3v_p2",
         "rt_triad_bgk_spherical_blast_1x3v_p1",
         "rt_escreen_sr",
         "rt_hyper_vlasov_tm",
         "rt_vlasov_kerntm",
         "rt_vlasov_moments",
         "rt_diffusion_1x",
      },
      gpu = {
         "rt_dg_5m_mom_beach_p2",
         "rt_vlasov_sr_freestream",
         "rt_dg_diffusion_gen_3x",
         "rt_escreen_sr",
         "rt_dg_5m_mom_beach_p3",
         "rt_hyper_vlasov_tm",
         "rt_dg_diffusion_gen_2x",
         "rt_diffusion_1x",
         "rt_vlasov_kerntm",
         "rt_diffusion_2x",
      },
   },
   moat = {
      "rt_vlasov_landau_damping_1x1v", "rt_vlasov_twostream_1x1v",
      "rt_vlasov_es_shock_1x1v", "rt_vlasov_bgk_relax_1x1v",
   },
   parallel = {
      { name = "rt_vlasov_twostream_1x1v", cuts = { 4 } },
      { name = "rt_vlasov_weibel_2x2v", cuts = { 2, 2 } },
   },
}
