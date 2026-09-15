#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_hyper_dg.h>

static struct gkyl_array*
mkarr1(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (use_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

//#define mkarr1(use_gpu, nc, size) (fprintf(stderr, "mkarr1: %d\n", __LINE__), mkarr1_(use_gpu, nc, size))

int hyper_dg_kernel_test(const gkyl_hyper_dg *slvr);

void
test_vlasov_3x3v_p1_(bool use_gpu)
{
  // initialize grid and ranges
  int cdim = 3, vdim = 3;
  int pdim = cdim+vdim;

  int cells[] = {8, 8, 8, 8, 8, 8};
  int ghost[] = {1, 1, 1, 0, 0, 0};
  double lower[] = {0., 0., 0., -1., -1., -1.};
  double upper[] = {1., 1., 1., 1., 1., 1.};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  // initialize basis
  int poly_order = 1;
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in velocity space). */
    gkyl_cart_modal_hybrid(&basis, cdim, vdim);
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // initialize eqn
  struct gkyl_dg_eqn *eqn;
  enum gkyl_field_id field_id = GKYL_FIELD_NULL;
  enum gkyl_model_id model_id = GKYL_MODEL_DEFAULT;
  eqn = gkyl_dg_vlasov_new(&confBasis, &basis, &confRange, &phaseRange, model_id, field_id, use_gpu);

  // initialize hyper_dg slvr 
  // FIELD_NULL so only configuration space update, no velocity space update
  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2};
  int zero_flux_flags[2*GKYL_MAX_DIM] = {0, 0, 0, 0, 0, 0};
  int num_up_dirs = cdim; 

  gkyl_hyper_dg *slvr;
  slvr = gkyl_hyper_dg_new(&phaseGrid, &basis, eqn, num_up_dirs, up_dirs, zero_flux_flags, 1, use_gpu);

  // initialize arrays
  struct gkyl_array *fin, *rhs, *cflrate;
  struct gkyl_array *fin_h, *rhs_h;
  
  fin = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  rhs = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  cflrate = mkarr1(use_gpu, 1, phaseRange_ext.volume);

  // set initial condition
  int nf = phaseRange_ext.volume*basis.num_basis;
  double *fin_d;
  if (use_gpu) {
    fin_h = mkarr1(false, basis.num_basis, phaseRange_ext.volume);
    fin_d = fin_h->data;
  } else {
    fin_d = fin->data;
  }
  for(int i=0; i< nf; i++) {
    fin_d[i] = (double)(2*i+11 % nf) / nf  * ((i%2 == 0) ? 1 : -1);
  }
  if (use_gpu) gkyl_array_copy(fin, fin_h);

  // run hyper_dg_advance
  int nrep = 10;
  for(int n=0; n<nrep; n++) {
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_vlasov_set_auxfields(eqn,
      (struct gkyl_dg_vlasov_auxfields){.field = 0, .cot_vec = 0, 
      .alpha_surf = 0, .sgn_alpha_surf = 0, .const_sgn_alpha = 0 }); // must set EM fields to use

    gkyl_hyper_dg_advance(slvr, &phaseRange, fin, cflrate, rhs);
  }

  // get linear index of first non-ghost cell
  // 1-indexed for interfacing with G2 Lua layer
  int idx[] = {1, 1, 1, 1, 1, 1};
  int linl = gkyl_range_idx(&phaseRange, idx);

  rhs_h = mkarr1(false, basis.num_basis, phaseRange_ext.volume);
  gkyl_array_copy(rhs_h, rhs);

  // check that ghost cells are empty
  double val = 0;
  double *rhs_d;
  int i = 0;
  while (val==0) {
    rhs_d = gkyl_array_fetch(rhs_h, i);
    val = rhs_d[0];
    if(val==0) i++;
  }
  TEST_CHECK(i == linl);

  // check data in first non-ghost cell
  rhs_d = gkyl_array_fetch(rhs_h, linl);

  // printf("\nfirst cell rhs\n");
  //for(int i=0; i<rhs->ncomp; i++) printf("%.16e\n", rhs_d[i]);
  //for(int i=0; i<rhs->ncomp; i++) printf("  TEST_CHECK( gkyl_compare_double(rhs_d[%d], %.16e, 1e-12) );\n", i, rhs_d[i]);
  TEST_CHECK( gkyl_compare_double(rhs_d[0], 3.8735549798502387e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[1], 1.6959259934263994e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[2], -4.8984121582866758e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[3], 5.0094536289262006e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[4], -1.2449694827016953e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[5], 3.5633095838952239e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[6], 3.4307333299092484e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[7], 2.5317237348814032e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[8], -2.5305142336381611e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[9], 1.6976562666517527e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[10], -1.3188217206179452e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[11], 1.1566279578878277e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[12], -1.0487240328165591e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[13], 9.4320746037960195e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[14], -9.4533250999594411e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[15], 1.0693888393608047e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[16], 3.3731060719665238e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[17], 1.4232063291190894e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[18], -1.1238200859276988e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[19], 8.4493461115432069e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[20], 3.7168983240080964e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[21], -3.4422939181393599e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[22], -3.5083447363027062e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[23], 2.1375305784636044e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[24], -2.0825739029971018e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[25], 1.7488681040300975e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[26], -2.0179782511044245e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[27], 2.0435157947684868e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[28], -1.6636355898951987e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[29], 1.3442391639742803e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[30], -1.0853876848206312e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[31], 1.0965697471325935e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[32], -2.0951764524448709e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[33], 1.9681787690876064e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[34], -1.6327928661319586e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[35], 1.5449923691405468e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[36], -1.0957323983061325e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[37], 9.7521904885660291e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[38], -1.4225790936968075e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[39], 1.0938543097641011e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[40], -1.0386797997392399e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[41], -3.4880098226412404e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[42], -3.2719711091781413e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[43], 3.4275232843717916e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[44], -2.2836372706768060e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[45], 2.0046113232551104e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[46], -1.6286898011836502e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[47], 3.6396222957973507e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[48], -1.8214304335963199e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[49], 1.9421272746405052e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[50], -2.1630418974003465e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[51], 2.3930071574221582e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[52], -2.2853209071137993e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[53], 2.3025619211231863e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[54], -1.5458366066464473e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[55], 7.1515312063012813e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[56], -7.0629012979087076e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[57], 3.5412075288314718e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[58], -3.2236018558034843e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[59], 2.6160256284435640e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[60], -2.3954488886068940e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[61], 2.4102710211474815e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[62], -1.2413660608970357e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[63], 3.1404312005369288e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[64], 3.8478463062984543e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[65], 1.6383893567995464e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[66], -9.8254789277861470e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[67], 1.0874361679351932e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[68], 3.4860754298984871e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[69], -4.0144054180961835e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[70], -2.2701650346193428e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[71], 2.0111297913693406e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[72], -2.0080280134615986e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[73], 1.4867323173528218e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[74], -4.7783457601632309e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[75], 4.8533857084924703e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[76], -1.7159060407537222e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[77], 6.5386083037954048e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[78], -6.4865099825069255e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[79], 9.2319002778701043e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[80], -2.9241883912606315e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[81], 2.7051634774399659e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[82], -2.6999541275348928e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[83], 1.6379088171963765e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[84], -2.3768994612383693e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[85], 2.3844039231662258e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[86], -1.3397291027561050e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[87], 1.1915163460450563e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[88], -1.1526371087392411e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[89], 8.9360187378766600e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[90], -3.5968321533258347e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[91], 3.5646551271490665e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[92], -1.9084162474935276e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[93], 2.0133047360266591e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[94], -1.5776559149735849e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[95], 3.3655193755196017e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[96], 3.4765338749530885e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[97], 1.4286349564885294e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[98], -1.0112338482133872e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[99], 1.1150930863744215e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[100], 3.4518779637733150e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[101], -3.8558779637733149e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[102], -2.1677973959573716e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[103], 2.0592000772753732e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[104], -1.8526642147186362e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[105], 1.4725144649176752e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[106], -4.8521922798635471e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[107], 4.8411612448282195e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[108], -1.6957738962444672e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[109], 5.0229451664849600e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[110], -6.5552330699909449e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[111], 9.5738749391218370e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[112], -3.0587723580006806e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[113], 2.7061391432184497e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[114], -2.7135368204699201e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[115], 1.6658891463003339e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[116], -2.5232806031213709e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[117], 2.3727686083107582e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[118], -1.5348828550137744e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[119], 1.2053183736637598e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[120], -1.1410362029896483e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[121], 8.9018513779220996e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[122], -3.5401022042449036e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[123], 3.3559361556510758e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[124], -1.9028085031190408e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[125], 1.9995127996263243e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[126], -1.5892669582706860e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[127], 3.3599151830280540e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[128], 3.4338880304476271e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[129], 1.4269014807977600e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[130], -1.1635169628913534e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[131], 9.6639648064887425e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[132], 3.4738880304476267e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[133], -3.8338880304476275e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[134], -2.1630741319338554e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[135], 2.0676527431156593e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[136], -1.8565090575750190e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[137], 1.4728594501959885e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[138], -4.7295692886806586e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[139], 4.9222976237120477e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[140], -1.6896409618195747e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[141], 6.5060001864290484e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[142], -5.0282663545078714e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[143], 9.7937742723787069e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[144], -3.0588791152729172e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[145], 2.7152337758829059e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[146], -2.7054978742680060e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[147], 1.6664574133151095e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[148], -2.3674068408778592e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[149], 2.5247172067018912e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[150], -1.5248151675977429e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[151], 1.2031229322141296e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[152], -9.9197930392987512e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[153], 1.0353108517654821e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[154], -3.4020199323497820e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[155], 3.3660199786841687e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[156], -2.0518750723189335e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[157], 1.8543965617470551e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[158], -1.5910122910463183e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[159], 3.3616635264207474e+01, 1e-12) );

/*   // get linear index of some other cell */
/*   // 1-indexed for interfacing with G2 Lua layer */
  int idx2[] = {6, 3, 5, 8, 2, 1};
  int linl2 = gkyl_range_idx(&phaseRange, idx2);
  rhs_d = gkyl_array_fetch(rhs_h, linl2);

  // printf("second cell rhs\n");
  //for(int i=0; i<rhs->ncomp; i++) printf("%.16e\n", rhs_d[i]);
  //for(int i=0; i<rhs->ncomp; i++) printf("  TEST_CHECK( gkyl_compare_double(rhs_d[%d], %.16e, 1e-12) );\n", i, rhs_d[i]);
  TEST_CHECK( gkyl_compare_double(rhs_d[0], 7.0283701215299565e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[1], 4.7846177739417136e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[2], -3.2440580569657463e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[3], 4.7834435172957100e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[4], -3.3964806386044817e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[5], 3.6656129906922845e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[6], 1.0600153622120099e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[7], 9.8197929504110462e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[8], -1.1346864925832745e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[9], 9.9519188298106030e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[10], -3.8427607132713874e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[11], 4.6685208259857504e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[12], -6.1166112907734686e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[13], 5.0806785012409009e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[14], -3.4404297692371081e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[15], 6.1730531849577105e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[16], 1.3601589806692580e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[17], 5.5606774714677393e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[18], -4.6853462430480434e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[19], 4.8960263481450980e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[20], 9.8840947991088379e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[21], -1.0715759504421176e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[22], -1.5649149976122467e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[23], 9.0646874012320907e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[24], -1.0551869062512243e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[25], 9.9534967043384000e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[26], -8.9669280348657580e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[27], 1.0859866842168806e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[28], -9.5787783660800088e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[29], 5.9003194017868232e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[30], -4.2858735788656148e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[31], 6.1644570050895027e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[32], -9.3693889519560514e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[33], 1.0459267060273295e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[34], -9.6279354665355200e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[35], 5.1838536940564353e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[36], -4.6930351157227243e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[37], 5.7893976582206406e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[38], -5.5600505574220719e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[39], 4.2501170277772061e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[40], -5.7370813214632818e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[41], -7.5952097854402589e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[42], -1.5112521038242042e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[43], 1.5166281201815158e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[44], -9.3845022817631502e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[45], 9.5950426777163528e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[46], -9.5028250006771700e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[47], 1.4343455956661325e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[48], -8.8669420429927314e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[49], 1.1596206998529810e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[50], -9.9656228044263742e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[51], 8.9575579031755979e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[52], -1.0405867197686671e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[53], 1.1273492436504714e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[54], -5.2423958891702576e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[55], 4.5661639299989652e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[56], -6.0855838051600962e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[57], 1.4781355074540281e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[58], -1.2684443586352376e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[59], 1.2637828515811871e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[60], -8.1030760766532168e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[61], 9.6461812198928470e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[62], -7.3547047664692542e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[63], 1.4322009894672004e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[64], 7.1935746220124452e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[65], 4.7913444708426276e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[66], -3.7409857712756967e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[67], 6.1153234258921039e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[68], 1.0731283386012138e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[69], -8.8591657399897272e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[70], -8.6398599759795090e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[71], 9.2095110512420121e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[72], -9.5968996273635057e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[73], 5.5514941660377154e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[74], -3.3078263954184003e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[75], 4.8436138128923325e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[76], -4.7834501948352298e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[77], 4.1396081713442648e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[78], -5.6626818497470516e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[79], 3.4747011836898047e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[80], -1.4264972931737645e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[81], 9.0887141852893933e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[82], -1.0611788447974614e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[83], 9.9679463000264775e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[84], -8.8872500641148591e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[85], 1.0423038256781317e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[86], -8.2151469282290449e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[87], 4.2729166191391300e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[88], -4.6960818962577818e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[89], 5.2657335222913041e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[90], -1.3954254583181825e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[91], 1.3993631687859335e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[92], -8.4764645663605506e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[93], 1.0850803211201104e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[94], -9.0949727090642483e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[95], 1.4720160200709708e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[96], 1.1058159072558467e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[97], 5.5661072828562872e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[98], -3.8038399453442359e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[99], 6.1387574319713288e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[100], 1.0811599960760687e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[101], -6.8515999607606970e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[102], -7.7181652109862142e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[103], 9.1271728627705045e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[104], -8.7568400872410479e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[105], 5.5961302495084126e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[106], -3.3578250907422031e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[107], 4.8466142788859123e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[108], -4.7844638330728891e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[109], 3.2949002915137193e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[110], -5.6653312461354503e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[111], 3.4666695262149414e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[112], -1.5147327839312067e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[113], 9.0008888762480964e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[114], -1.0574963170297912e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[115], 9.8817565354966675e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[116], -9.7729562651682940e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[117], 1.0461811224142504e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[118], -9.0234403715862186e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[119], 4.2278649030128754e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[120], -4.6503108968744243e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[121], 5.2665396986558520e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[122], -1.5562513629902506e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[123], 1.4807822157984455e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[124], -8.5681017412682138e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[125], 1.0895865064874867e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[126], -9.1407519709951146e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[127], 1.4811801000515698e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[128], 1.0631700627503820e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[129], 5.5643743484313944e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[130], -4.6450398387617000e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[131], 5.3011431100062829e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[132], 1.1031700627503864e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[133], -6.6317006275038415e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[134], -7.6750483997030528e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[135], 9.1740170558386822e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[136], -8.7606839925974327e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[137], 5.5964757760526041e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[138], -3.3071697856301448e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[139], 4.8931199852680621e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[140], -4.7783303573821193e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[141], 4.1321225722476250e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[142], -4.8237168583476468e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[143], 3.4886594595406240e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[144], -1.5147434200350182e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[145], 9.0483740986404442e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[146], -1.0528531614336353e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[147], 9.8823238650114433e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[148], -8.9281643904511640e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[149], 1.1302678080039013e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[150], -9.0133717466701867e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[151], 4.2256689202973689e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[152], -3.8123353440751522e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[153], 6.1005831288686224e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[154], -1.4773905314027536e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[155], 1.4817905584783423e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[156], -9.4060856304734799e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[157], 1.0061831652021979e+02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[158], -9.1424963662707484e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[159], 1.4813548947674269e+02, 1e-12) );

  // clean up
  gkyl_array_release(fin);
  gkyl_array_release(rhs);
  gkyl_array_release(rhs_h);
  gkyl_array_release(cflrate);

  gkyl_hyper_dg_release(slvr);
  gkyl_dg_eqn_release(eqn);

  if (use_gpu) {
    gkyl_array_release(fin_h);
  }
}

void
test_hyper3x_dg_vlasov_3x3v_p1_ho()
{
  test_vlasov_3x3v_p1_(false);
}

TEST_LIST = {
  { "test_hyper3x_dg_vlasov_3x3v_p1_ho", test_hyper3x_dg_vlasov_3x3v_p1_ho },
  { NULL, NULL },
};

