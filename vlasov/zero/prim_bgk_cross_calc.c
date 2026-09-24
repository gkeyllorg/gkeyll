#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_prim_bgk_cross_calc.h>

void
gkyl_prim_bgk_cross_calc_advance(struct gkyl_basis basis,
  int vdim_phys, const struct gkyl_array* m0sdeltas,
  double massself, const struct gkyl_array* primsself,
  double massother, const struct gkyl_array* primsother,
  const struct gkyl_range *range, struct gkyl_array* crossprims)
{
  unsigned num_basis = basis.num_basis;
  assert(num_basis <= 20); // MF 2022/11/20: Hardcoded to a max of 3x p2 ser.
  unsigned ndim = basis.ndim;
  unsigned poly_order = basis.poly_order;
  mul_op_t mul_op = choose_ser_mul_kern(ndim, poly_order);

  unsigned udim = primsself->ncomp/num_basis-1;
  unsigned u_num_basis = udim*num_basis;
  double massDiff = 0.5*(massself-massother)/vdim_phys;
  double massSum = massself+massother;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *m0sdeltas_d = gkyl_array_cfetch(m0sdeltas, loc);
    const double *primsself_d = gkyl_array_cfetch(primsself, loc);
    const double *primsother_d = gkyl_array_cfetch(primsother, loc);

    double *crossprims_d = gkyl_array_fetch(crossprims, loc);

    // Pointers to each primitive moment (for simplicity).
    const double *uself = primsself_d;
    const double *vtsqself = &primsself_d[u_num_basis];
    const double *uother = primsother_d;
    const double *vtsqother = &primsother_d[u_num_basis];
    double *ucross = crossprims_d;
    double *vtsqcross = &crossprims_d[u_num_basis];

    // u_sr = u_si - u_ri:
    array_set1(u_num_basis, ucross, 1., uself);
    array_acc1(u_num_basis, ucross,-1., uother);

    // v_tsr^2 = (u_si-u_ri)^2 = (u_si-u_ri) . (u_si-u_ri)
    double vbuf[20*3]; // MF 2022/11/20: Hardcoded to 3x p2 ser (3 components).
    for (int k=0; k<num_basis; k++) vtsqcross[k] = 0.;
    for (int d=0; d<udim; d++) {
      mul_op(ucross+d*num_basis, ucross+d*num_basis, vbuf);
      for (int k=0; k<num_basis; k++) vtsqcross[k] += vbuf[k];
    }

    // v_tsr^2 = m_s*v_ts^2-m_r*v_tr^2+((m_s-m_r)/(2*vdim_phys))*(u_si-u_ri)^2
    for (int k=0; k<num_basis; k++)
      vtsqcross[k] = massDiff*vtsqcross[k]+massself*vtsqself[k]-massother*vtsqother[k];;

    // vbuf = -0.5*delta_s*(beta+1)*(u_si - u_ri) = u_sri - u_si:
    for (int d=0; d<udim; d++) {
      mul_op(m0sdeltas_d, ucross+d*num_basis, vbuf+d*num_basis);
      for (int k=d*num_basis; k<(d+1)*num_basis; k++)
        vbuf[k] *= -0.5;
    }
    // v_tsr^2 = -(delta_s*(beta+1)/(m_s+m_r))*(m_s*v_ts^2-m_r*v_tr^2+((m_s-m_r)/(2*vdim_phys))*(u_si-u_ri)^2)
    mul_op(m0sdeltas_d, vtsqcross, vtsqcross);
    for (int k=0; k<num_basis; k++)
      vtsqcross[k] = -vtsqcross[k]/massSum;

    // u_sr = u_si-0.5*delta_s*(beta+1)*(u_si - u_ri) = u_si + vbuf:
    for (int k=0; k<u_num_basis; k++)
      ucross[k] = uself[k]+vbuf[k];

    // v_tsr^2 = -(u_sri-u_ri).(u_sri-u_si)/vdim_phys
    //   -(delta_s*(beta+1)/(m_s+m_r))*(m_s*v_ts^2-m_r*v_tr^2+((m_s-m_r)/(2*vdim_phys))*(u_si-u_ri)^2)
    // Need the dot product (u_sri-u_ri).(u_sri-u_si):
    for (int d=0; d<udim; d++) {
      double compbuf[20]; // MF 2022/11/20: Hardcoded to 3x p2 ser.
      for (int k=0; k<num_basis; k++)
        compbuf[k] = ucross[d*num_basis+k]-uother[d*num_basis+k];
      mul_op(compbuf, vbuf+d*num_basis, compbuf);
      for (int k=0; k<num_basis; k++)
        vtsqcross[k] += -compbuf[k]/vdim_phys;
    }

    // v_tsr^2 = v_ts^2-(u_sri-u_ri).(u_sri-u_si)/vdim_phys
    //   -(delta_s*(beta+1)/(m_s+m_r))*(m_s*v_ts^2-m_r*v_tr^2+((m_s-m_r)/(2*vdim_phys))*(u_si-u_ri)^2)
    for (int k=0; k<num_basis; k++)
      vtsqcross[k] += vtsqself[k];
  }
}
