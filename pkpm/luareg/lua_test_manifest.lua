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
      "rt_pkpm_landau_damping_1x", "rt_pkpm_em_advect_1x",
      "rt_pkpm_neut_lbo_sodshock_1x", "rt_pkpm_travel_pulse_lbo_1x",
   },
   parallel = {},
}
