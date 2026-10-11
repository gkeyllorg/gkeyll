-- Tests skipped by the C regression suite.
-- Remove an entry manually to re-enable the test.
-- gpu: C tests whose GPU variant timed out (CPU variant still runs).
return {
   ignore = {
      tests = {
      },
      gpu = {
         "rt_dg_diffusion_gen_3x",
         "rt_dg_diffusion_gen_2x",
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
