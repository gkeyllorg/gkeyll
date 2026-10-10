-- Tests skipped by the Lua regression suite.
-- Remove an entry manually to re-enable the test.
-- gpu: Lua tests whose GPU variant timed out (CPU variant still runs).
return {
   ignore = {
      tests = {
      },
      gpu = {
      },
   },
   moat = {
      "rt_vlasov_landau_damping_1x1v_ser_p2", "rt_vlasov_twostream_1x1v_ser_p2",
      "rt_vlasov_es_shock_1x1v_ser_p2", "rt_vlasov_bgk_relax_1x1v_ser_p2",
   },
   parallel = {},
}
