-- Tests skipped by the C regression suite.
-- Remove an entry manually to re-enable the test.
-- gpu: C tests whose GPU variant timed out (CPU variant still runs).
return {
   ignore = {
      tests = {
         "rt_pkpm_alf_wave_p1",
         "rt_pkpm_alf_wave_explicit_1x_p1",
         "rt_pkpm_alf_wave_1x_p1",
         "rt_pkpm_alf_wave_explicit_1x_p2",
         "rt_pkpm_periodic_es_shock_p2",
         "rt_pkpm_ot_p1",
         "rt_pkpm_alf_wave_1x_p2",
         "rt_pkpm_es_shock_p2",
      },
      gpu = {
         "rt_pkpm_alf_wave_explicit_1x_p2",
         "rt_pkpm_alf_wave_1x_p2",
         "rt_pkpm_alf_wave_explicit_1x_p1",
         "rt_pkpm_alf_wave_1x_p1",
         "rt_pkpm_travel_pulse_p2",
         "rt_pkpm_alf_soliton_1x_p2",
      },
   },
   moat = {
      "rt_pkpm_landau_damping_1x", "rt_pkpm_em_advect_1x",
      "rt_pkpm_neut_lbo_sodshock_1x", "rt_pkpm_travel_pulse_lbo_1x",
   },
   parallel = {
      { name = "rt_pkpm_alf_soliton_1x_p2", cuts = { 4 } },
      { name = "rt_pkpm_travel_pulse_lbo_2x", cuts = { 2, 2 } },
   },
}
