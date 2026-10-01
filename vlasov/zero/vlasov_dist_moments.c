#include <gkyl_vlasov_dist_moments.h>
#include <gkyl_vlasov_dist_moments_priv.h>

static dist_moments_new_t dist_moments_new[] ={
  [GKYL_DIST_TYPE_LTE] = gkyl_vlasov_lte_moments_new,
};

struct gkyl_vlasov_dist_moments* gkyl_vlasov_dist_moments_inew(const struct gkyl_vlasov_dist_moments_inp *inp) 
{
  return dist_moments_new[inp->dist_id](inp);
}


void 
gkyl_vlasov_dist_density_moment_advance(struct gkyl_vlasov_dist_moments *dist_moms, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local, 
  const struct gkyl_array *fin, struct gkyl_array *density_out)
{
  dist_moms->density_moment(dist_moms, phase_local, conf_local, fin, density_out);
}

void 
gkyl_vlasov_dist_moments_advance(struct gkyl_vlasov_dist_moments *dist_moms, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local, 
  const struct gkyl_array *fin, struct gkyl_array *moms_out)
{
  dist_moms->moments(dist_moms, phase_local, conf_local, fin, moms_out);
}

void 
gkyl_vlasov_dist_moments_release(struct gkyl_vlasov_dist_moments *dist_moms)
{
  dist_moms->release(dist_moms);
}
