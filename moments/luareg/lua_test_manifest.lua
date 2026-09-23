-- Tests skipped by the Lua regression suite.
-- Remove an entry manually to re-enable the test.
-- gpu: Lua tests whose GPU variant timed out (CPU variant still runs).
return {
   ignore = {
      tests = {
         "rt_gr_ultra_rel_bhl_static",
         "rt_gr_bhl_spinning_neutronstar",
         "rt_gr_mhd_bhl_static_tetrad",
         "rt_gr_mhd_bhl_static",
         "rt_gr_mhd_bhl_spinning",
         "rt_gr_ultra_rel_bhl_static_tetrad",
         "rt_gr_bhl_static",
         "rt_gr_multifluid_bhl_spinning",
         "rt_gr_mhd_brio_wu_tetrad",
         "rt_gr_multifluid_bhl_spinning_tetrad",
         "rt_gr_multifluid_bhl_static_tetrad",
         "rt_gr_ultra_rel_bhl_spinning",
         "rt_gr_mhd_strong_blast_tetrad",
         "rt_euler_rgfm_shock_bubble",
         "rt_10m_ot_grad_closure",
         "rt_gr_mhd_brio_wu",
         "rt_gr_bhl_static_neutronstar",
         "rt_gr_ultra_rel_bhl_spinning_tetrad",
         "rt_gr_mhd_bhl_spinning_tetrad",
         "rt_gr_multifluid_bhl_static",
         "rt_gr_bhl_spinning",
         "rt_gr_mhd_blackhole_collapse",
         "rt_10m_ot",
         "rt_10m_burch_grad_closure",
         "rt_gr_mhd_strong_blast",
         "rt_gr_kh_2d",
         "rt_10m_ot_nn_closure",
         "rt_10m_burch",
         "rt_gr_wald_magnetosphere_spinning_neutronstar",
         "rt_iso_euler_mixture_shock_bubble",
         "rt_gr_wald_magnetosphere_static_neutronstar",
      },
      gpu = {
      },
   },
   moat = {
      "rt_euler_sodshock", "rt_mhd_brio_wu", "rt_10m_riem", "rt_5m_gem",
      "rt_maxwell_plane_wave_1d", "rt_gr_einstein_plane_shock",
   },
   parallel = {},
}
