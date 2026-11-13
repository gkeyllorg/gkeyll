-- Tests skipped by the Lua regression suite.
-- Remove an entry manually to re-enable the test.
-- gpu: Lua tests whose GPU variant does not run (CPU variant still runs). The general
-- diffusion tests are listed because the general-stencil hyper DG advance is not
-- implemented on the device (the updater asserts there), not because of a timeout.
return {
   ignore = {
      tests = {
      },
      gpu = {
         "rt_dg_diffusion_gen_2x_ser_p1",
         "rt_dg_diffusion_gen_2x_ser_p2",
         "rt_dg_diffusion_gen_2x_ser_p3",
         "rt_dg_diffusion_gen_2x_tensor_p2",
         "rt_dg_diffusion_gen_2x_tensor_p3",
         "rt_dg_diffusion_gen_3x_ser_p1",
         "rt_dg_diffusion_gen_3x_ser_p2",
      },
   },
   moat = {
      "rt_vlasov_landau_damping_1x1v_ser_p2", "rt_vlasov_twostream_1x1v_ser_p2",
      "rt_vlasov_es_shock_1x1v_ser_p2", "rt_vlasov_bgk_relax_1x1v_ser_p2",
   },
   parallel = {},
}
