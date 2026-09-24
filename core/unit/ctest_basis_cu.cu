/* -*- c -*- */

#include <stdio.h>

extern "C" {
#include <gkyl_basis.h>
#include <gkyl_util.h>
#include <gkyl_alloc.h>
  int dev_cu_ser_2d(struct gkyl_basis *basis);
}

__global__
void
ker_dev_cu_ser_2d(struct gkyl_basis *basis, int *nfail)
{
  *nfail = 0;

  GKYL_CU_CHECK( basis->ndim == 2, nfail);
  GKYL_CU_CHECK( basis->poly_order == 2, nfail);
  GKYL_CU_CHECK( basis->num_basis == 8, nfail);
  GKYL_CU_CHECK( basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY, nfail);

  double z[128], b[128];

  z[0] = 0.0; z[1] = 0.0;
  basis->eval(z, b);

  double fin[128], fout[128];
  for (int i=0; i<basis->num_basis; ++i) {
    fin[i] = 1.0;
    fout[i] = 0.0;
  }

  basis->flip_odd_sign(0, fin, fout);
  GKYL_CU_CHECK( fin[0] == fout[0], nfail);
  GKYL_CU_CHECK( -fin[1] == fout[1], nfail);
  GKYL_CU_CHECK( fin[2] == fout[2], nfail);  
  GKYL_CU_CHECK( -fin[3] == fout[3], nfail);
  GKYL_CU_CHECK( fin[4] == fout[4], nfail);
  GKYL_CU_CHECK( fin[5] == fout[5], nfail);
  GKYL_CU_CHECK( fin[6] == fout[6], nfail);
  GKYL_CU_CHECK( -fin[7] == fout[7], nfail);

  basis->flip_odd_sign(1, fin, fout);
  GKYL_CU_CHECK( fin[0] == fout[0], nfail);
  GKYL_CU_CHECK( fin[1] == fout[1], nfail);
  GKYL_CU_CHECK( -fin[2] == fout[2], nfail);  
  GKYL_CU_CHECK( -fin[3] == fout[3], nfail);
  GKYL_CU_CHECK( fin[4] == fout[4], nfail);
  GKYL_CU_CHECK( fin[5] == fout[5], nfail);
  GKYL_CU_CHECK( -fin[6] == fout[6], nfail);
  GKYL_CU_CHECK( fin[7] == fout[7], nfail);

  basis->flip_even_sign(0, fin, fout);
  GKYL_CU_CHECK( -fin[0] == fout[0], nfail);
  GKYL_CU_CHECK( fin[1] == fout[1], nfail);
  GKYL_CU_CHECK( -fin[2] == fout[2], nfail);  
  GKYL_CU_CHECK( fin[3] == fout[3], nfail);
  GKYL_CU_CHECK( -fin[4] == fout[4], nfail);
  GKYL_CU_CHECK( -fin[5] == fout[5], nfail);
  GKYL_CU_CHECK( -fin[6] == fout[6], nfail);
  GKYL_CU_CHECK( fin[7] == fout[7], nfail);

  basis->flip_even_sign(1, fin, fout);
  GKYL_CU_CHECK( -fin[0] == fout[0], nfail);
  GKYL_CU_CHECK( -fin[1] == fout[1], nfail);
  GKYL_CU_CHECK( fin[2] == fout[2], nfail);  
  GKYL_CU_CHECK( fin[3] == fout[3], nfail);
  GKYL_CU_CHECK( -fin[4] == fout[4], nfail);
  GKYL_CU_CHECK( -fin[5] == fout[5], nfail);
  GKYL_CU_CHECK( fin[6] == fout[6], nfail);
  GKYL_CU_CHECK( -fin[7] == fout[7], nfail);  

  double nodes[128*128];
  basis->node_list(nodes);

  GKYL_CU_CHECK( nodes[0] == -1, nfail);
  GKYL_CU_CHECK( nodes[1] == -1, nfail);

  GKYL_CU_CHECK( nodes[2] == 0, nfail);
  GKYL_CU_CHECK( nodes[3] == -1, nfail);

  GKYL_CU_CHECK( nodes[4] == 1, nfail);
  GKYL_CU_CHECK( nodes[5] == -1, nfail);
}

int
dev_cu_ser_2d(struct gkyl_basis *basis)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));
  ker_dev_cu_ser_2d<<<1,1>>>(basis, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;
}
