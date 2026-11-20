-- Tests skipped by the C regression suite.
-- Remove an entry manually to re-enable the test.
-- gpu: C tests whose GPU variant does not run (CPU variant still runs). The general
-- diffusion tests are listed because the general-stencil hyper DG advance is not
-- implemented on the device (the updater asserts there), not because of a timeout.
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
      "rt_pkpm_landau_damping_1x", "rt_pkpm_em_advect_1x",
      "rt_pkpm_neut_lbo_sodshock_1x", "rt_pkpm_travel_pulse_lbo_1x",
   },
   parallel = {
      { name = "rt_vlasov_twostream_1x1v", cuts = { 4 } },
      { name = "rt_vlasov_weibel_2x2v", cuts = { 2, 2 } },
      { name = "rt_pkpm_travel_pulse_lbo_2x", cuts = { 2, 2 } },
   },
}
