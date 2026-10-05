-- Tests skipped by the Lua regression suite.
-- Remove an entry manually to re-enable the test.
-- gpu: Lua tests whose GPU variant timed out (CPU variant still runs).
return {
   ignore = {
      tests = {
         "rt_pkpm_alf_wave_1x_p1",
         "rt_pkpm_2d_travel_pulse_p1",
         "rt_pkpm_alf_wave_1x_p2",
         "rt_pkpm_ot_p1",
         "rt_pkpm_alf_wave_explicit_1x_p1",
         "rt_pkpm_alf_wave_explicit_1x_p2",
         "rt_pkpm_periodic_es_shock_p2",
         "rt_pkpm_es_shock_p2",
         "rt_pkpm_periodic_es_shock_p1",
      },
      gpu = {
      },
   },
   moat = {
      "rt_pkpm_landau_damping_p1", "rt_pkpm_em_advect_p1",
      "rt_pkpm_neut_sodshock_p1", "rt_pkpm_travel_pulse_p1",
   },
   parallel = {},
}
