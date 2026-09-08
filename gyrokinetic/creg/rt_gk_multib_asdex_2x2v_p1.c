#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_efit.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_mpi_comm.h>
#include <gkyl_null_comm.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>

struct gk_asdex_ctx {
  int cdim, vdim; // Dimensionality.

  double charge_elc; // Electron charge.
  double charge_ion; // Ion charge.
  double mass_elc; // Electron mass.
  double mass_ion; // Ion mass.
  double Te; // Electron temperature
  double Ti; // Ion temperature
  double c_s; // sound speed
  double nu_elc; // electron collision frequency
  double nu_ion; // ion collision frequency
  double B0; // Magnetic field.
  double n0; // Density.

  // Source parameters
  double psi_src; // Location.
  double lambda_src; // Width.
  double ndot_src; // Particle source rate.
  double Te_src; // Electron source temperature.
  double Ti_src; // Ion source temperature.

  // Domain parameters.            
  char geqdsk_file[128]; // File with equilibrium.
  double psi_axis; // Psi at the magnetic axis.
  double psi_sep; // Psi at the separatrix.
  double psi_min_core, psi_max_sol, psi_min_pf; // Psi extents.
  double Lx_core; // Box size in x
  // Z location of the X-point on psi=psi_min.
  double z_xpt_psi_sep_lo, z_xpt_psi_sep_up;

  // Grid.
  int Npsi_sol; // Number of cells in psi in the SOL.
  int Npsi_pf; // Number of cells in psi in the private flux.
  int Npsi_core; // Number of cells in psi in the core.
  int Ntheta_divertor; // Number of cells in theta in the divertor.
  int Ntheta_sol; // Number of cells in theta in the (upper) SOL (and core).
  int Nvpar, Nmu; // Number of cells in vpar,mu.
  int cells_v[2]; // Number of cells in all directions.
  int num_blocks; // Number of blocks.

  // Physical velocity space limits
  double vpar_max_elc; // Parallel velocity extents for electrons.
  double mu_max_elc; // Maximum magnetic moment for electrons.
  double vpar_max_ion; // Parallel velocity extents for ions.
  double mu_max_ion; // Maximum magnetic moment for ions.

  // Computational velocity space limits
  double vpar_min_elc_c, vpar_max_elc_c;
  double mu_min_elc_c, mu_max_elc_c;
  double vpar_min_ion_c, vpar_max_ion_c;
  double mu_min_ion_c, mu_max_ion_c;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

void divertor_plate_func_out(double s, double* RZ)
{
//  // Straight plate from (1.5966,-1.1421) to (1.6888,-0.8781).
//  RZ[0] = 1.5966 + (1.6888 - 1.5966)*s;
//  RZ[1] = -1.1421 + (-0.8781 - (-1.1421))*s;
  // Cubic spline approximation to the plate.
  const int npts = 100;
//  double t[] = {
//    0.        , 0.01010101, 0.02020202, 0.03030303, 0.04040404, 0.05050505,
//    0.06060606, 0.07070707, 0.08080808, 0.09090909, 0.1010101 , 0.11111111,
//    0.12121212, 0.13131313, 0.14141414, 0.15151515, 0.16161616, 0.17171717,
//    0.18181818, 0.19191919, 0.2020202 , 0.21212121, 0.22222222, 0.23232323,
//    0.24242424, 0.25252525, 0.26262626, 0.27272727, 0.28282828, 0.29292929,
//    0.3030303 , 0.31313131, 0.32323232, 0.33333333, 0.34343434, 0.35353535,
//    0.36363636, 0.37373737, 0.38383838, 0.39393939, 0.4040404 , 0.41414141,
//    0.42424242, 0.43434343, 0.44444444, 0.45454545, 0.46464646, 0.47474747,
//    0.48484848, 0.49494949, 0.50505051, 0.51515152, 0.52525253, 0.53535354,
//    0.54545455, 0.55555556, 0.56565657, 0.57575758, 0.58585859, 0.5959596 ,
//    0.60606061, 0.61616162, 0.62626263, 0.63636364, 0.64646465, 0.65656566,
//    0.66666667, 0.67676768, 0.68686869, 0.6969697 , 0.70707071, 0.71717172,
//    0.72727273, 0.73737374, 0.74747475, 0.75757576, 0.76767677, 0.77777778,
//    0.78787879, 0.7979798 , 0.80808081, 0.81818182, 0.82828283, 0.83838384,
//    0.84848485, 0.85858586, 0.86868687, 0.87878788, 0.88888889, 0.8989899 ,
//    0.90909091, 0.91919192, 0.92929293, 0.93939394, 0.94949495, 0.95959596,
//    0.96969697, 0.97979798, 0.98989899, 1.        ,
//  };
//  double R[] = {
//    1.58211701, 1.58373013, 1.58522876, 1.58662408, 1.58792727, 1.58914951,
//    1.59030198, 1.59139587, 1.59244236, 1.59345262, 1.59443785, 1.59540921,
//    1.59637789, 1.59735404, 1.59833796, 1.59932451, 1.60030849, 1.60128472,
//    1.60224816, 1.60320156, 1.60415862, 1.60513384, 1.6061417 , 1.60719629,
//    1.608296  , 1.6094194 , 1.61054379, 1.61164649, 1.61270478, 1.61369639,
//    1.61461937, 1.61550099, 1.6163708 , 1.61725834, 1.61819317, 1.6192048 ,
//    1.62030643, 1.6214671 , 1.62264849, 1.62381228, 1.62492016, 1.62595178,
//    1.62693809, 1.62791909, 1.62893481, 1.63002528, 1.63123012, 1.63256489,
//    1.63400806, 1.63553491, 1.63712072, 1.63874077, 1.64037203, 1.642006  ,
//    1.64364166, 1.64527805, 1.64691421, 1.64854896, 1.65017582, 1.65178259,
//    1.65335681, 1.65488603, 1.65635783, 1.65777956, 1.65920642, 1.66070072,
//    1.66232475, 1.66412714, 1.66609031, 1.66817673, 1.67034882, 1.67256901,
//    1.67480116, 1.67704013, 1.67931122, 1.681641  , 1.68405605, 1.68658294,
//    1.68924066, 1.69201572, 1.69488593, 1.69782911, 1.70082309, 1.70384595,
//    1.70688438, 1.70993515, 1.7129956 , 1.71606307, 1.71913489, 1.72220864,
//    1.72528362, 1.72836005, 1.73143813, 1.73451804, 1.7376    , 1.7406842 ,
//    1.74377084, 1.74686013, 1.74995225, 1.75304743,
//  };
//  double Z[] = {
//    -1.1797302 , -1.17584427, -1.17192351, -1.16797136, -1.16399128, -1.15998669,
//    -1.15596104, -1.15191777, -1.14786031, -1.14379211, -1.13971661, -1.13563724,
//    -1.13155745, -1.12748036, -1.12340608, -1.11933308, -1.11525979, -1.11118464,
//    -1.10710612, -1.10302482, -1.09894425, -1.09486817, -1.09080032, -1.08674434,
//    -1.08269988, -1.0786615 , -1.07462347, -1.07058004, -1.06652547, -1.06245411,
//    -1.05836552, -1.05426666, -1.05016509, -1.04606836, -1.04198404, -1.03791967,
//    -1.03387851, -1.02985215, -1.02583029, -1.02180258, -1.01775872, -1.01369345,
//    -1.00961605, -1.00553836, -1.00147223, -0.99742949, -0.99342188, -0.98945418,
//    -0.9855205 , -0.981614  , -0.97772786, -0.97385526, -0.96998977, -0.96612844,
//    -0.96227006, -0.95841347, -0.95455747, -0.95070083, -0.94684033, -0.94297071,
//    -0.93908658, -0.93518255, -0.93125327, -0.92730239, -0.92335531, -0.91944071,
//    -0.91558724, -0.91181723, -0.90812254, -0.90448578, -0.90088956, -0.8973165 ,
//    -0.89375003, -0.8901916 , -0.88666037, -0.88317622, -0.87975903, -0.87642869,
//    -0.87319966, -0.87006321, -0.86700441, -0.8640083 , -0.86105996, -0.85814459,
//    -0.85525204, -0.85237755, -0.84951673, -0.84666513, -0.84381835, -0.84097233,
//    -0.83812611, -0.83528023, -0.83243524, -0.82959169, -0.82675013, -0.82391111,
//    -0.82107517, -0.81824287, -0.81541475, -0.81259136,
//  };
  // Extended plates:
  double t[] = {
    0.        , 0.01010101, 0.02020202, 0.03030303, 0.04040404, 0.05050505,
    0.06060606, 0.07070707, 0.08080808, 0.09090909, 0.1010101 , 0.11111111,
    0.12121212, 0.13131313, 0.14141414, 0.15151515, 0.16161616, 0.17171717,
    0.18181818, 0.19191919, 0.2020202 , 0.21212121, 0.22222222, 0.23232323,
    0.24242424, 0.25252525, 0.26262626, 0.27272727, 0.28282828, 0.29292929,
    0.3030303 , 0.31313131, 0.32323232, 0.33333333, 0.34343434, 0.35353535,
    0.36363636, 0.37373737, 0.38383838, 0.39393939, 0.4040404 , 0.41414141,
    0.42424242, 0.43434343, 0.44444444, 0.45454545, 0.46464646, 0.47474747,
    0.48484848, 0.49494949, 0.50505051, 0.51515152, 0.52525253, 0.53535354,
    0.54545455, 0.55555556, 0.56565657, 0.57575758, 0.58585859, 0.5959596 ,
    0.60606061, 0.61616162, 0.62626263, 0.63636364, 0.64646465, 0.65656566,
    0.66666667, 0.67676768, 0.68686869, 0.6969697 , 0.70707071, 0.71717172,
    0.72727273, 0.73737374, 0.74747475, 0.75757576, 0.76767677, 0.77777778,
    0.78787879, 0.7979798 , 0.80808081, 0.81818182, 0.82828283, 0.83838384,
    0.84848485, 0.85858586, 0.86868687, 0.87878788, 0.88888889, 0.8989899 ,
    0.90909091, 0.91919192, 0.92929293, 0.93939394, 0.94949495, 0.95959596,
    0.96969697, 0.97979798, 0.98989899, 1.        ,
  };
  double R[] = {
    1.57598284, 1.57749872, 1.57719396, 1.57796683, 1.57991018, 1.58038241,
    1.5809759 , 1.58399201, 1.58425793, 1.58426253, 1.58499896, 1.58614281,
    1.58753427, 1.58902876, 1.5904817 , 1.59179306, 1.59302565, 1.59427949,
    1.59565459, 1.59723624, 1.59893569, 1.60055062, 1.60187675, 1.60274568,
    1.60335474, 1.60411485, 1.60542574, 1.60728998, 1.60926429, 1.61088179,
    1.61186578, 1.61252968, 1.61330099, 1.61451591, 1.61611502, 1.61793178,
    1.61980027, 1.62157904, 1.62315876, 1.62443226, 1.6252949 , 1.62583334,
    1.62645555, 1.62760027, 1.62957002, 1.63208393, 1.63470503, 1.63700839,
    1.63888443, 1.64056222, 1.6422873 , 1.64422197, 1.6463157 , 1.64848453,
    1.65066913, 1.65286514, 1.65507565, 1.65730495, 1.65958062, 1.66195139,
    1.66446661, 1.66714803, 1.6699525 , 1.6728276 , 1.67572691, 1.67866399,
    1.68168688, 1.68484367, 1.68816064, 1.69163081, 1.69524436, 1.69899094,
    1.70285788, 1.70683192, 1.71089965, 1.7150371 , 1.71920156, 1.72334843,
    1.72743844, 1.73147165, 1.73546566, 1.7394382 , 1.7434051 , 1.74737814,
    1.7513686 , 1.75538652, 1.75942234, 1.76345072, 1.76744589, 1.77138913,
    1.77529078, 1.7791686 , 1.78304031, 1.78692   , 1.79081551, 1.79473402,
    1.79868275, 1.80266888, 1.80669962, 1.81078217,
  };
  double Z[] = {
    -1.23058861, -1.22520737, -1.21958678, -1.21420037, -1.2089721 , -1.20333331,
    -1.19784428, -1.19319607, -1.18774139, -1.18215984, -1.17672523, -1.171336  ,
    -1.1659685 , -1.16060623, -1.15523274, -1.14983804, -1.14443584, -1.13904528,
    -1.1336855 , -1.12837252, -1.12308551, -1.11777961, -1.11240953, -1.10693817,
    -1.10141171, -1.09592495, -1.09057007, -1.08534602, -1.08014713, -1.07486217,
    -1.069425  , -1.06390939, -1.05841614, -1.05302462, -1.0477214 , -1.04246792,
    -1.03722574, -1.03196133, -1.02664761, -1.02125793, -1.01576636, -1.01019997,
    -1.00467492, -0.99931589, -0.99420628, -0.98925267, -0.98431435, -0.97925461,
    -0.97404232, -0.96875963, -0.96349423, -0.95830444, -0.95317355, -0.94807301,
    -0.94298211, -0.93789756, -0.9328184 , -0.92774445, -0.92269006, -0.91768278,
    -0.91275056, -0.90790634, -0.90312784, -0.8983877 , -0.89366215, -0.88896287,
    -0.88432198, -0.87977162, -0.8753329 , -0.87101007, -0.86680593, -0.86272171,
    -0.85875197, -0.85488954, -0.85112709, -0.84744285, -0.84378957, -0.84011741,
    -0.8363825 , -0.83258511, -0.82874522, -0.82488291, -0.8210157 , -0.81715554,
    -0.81331368, -0.80950011, -0.80570491, -0.80190214, -0.79806541, -0.79417556,
    -0.79024315, -0.78628637, -0.78232333, -0.77836847, -0.77442978, -0.77051462,
    -0.76663036, -0.76278437, -0.75898399, -0.7552366 ,
  };
  // Find indices in t that bound s.
  int idx_tlo, idx_tup;
  if (s < 1e-8) {
    idx_tlo = 0;
    idx_tup = 0;
  }
  else if (fabs(s-1.0) < 1e-8) {
    idx_tlo = npts-1;
    idx_tup = npts-1;
  }
  else {
    for (int i=0; i<npts-1; i++) {
      if (t[i] <= s && s < t[i+1]) {
        idx_tlo = i;
        idx_tup = i+1;
        break;
      }
    }
  }
  // Interpolate the value of R and Z.
  double Dt = t[idx_tup]-t[idx_tlo];
  if (idx_tlo == idx_tup) {
    RZ[0] = R[idx_tlo];
    RZ[1] = Z[idx_tlo];
  }
  else {
    RZ[0] = ((s-t[idx_tlo])/Dt)*R[idx_tup] + ((t[idx_tup]-s)/Dt)*R[idx_tlo];
    RZ[1] = ((s-t[idx_tlo])/Dt)*Z[idx_tup] + ((t[idx_tup]-s)/Dt)*Z[idx_tlo];
  }
}

void divertor_plate_func_in(double s, double* RZ)
{
//  // Straight plate from (1.2686,-1.0520) to (1.1886,-0.7294).
//  RZ[0] = 1.2686 + (1.1886 - 1.2686)*s;
//  RZ[1] = -1.0520 + (-0.7294 - (-1.0520))*s;
  // Cubic spline approximation to the plate.
  const int npts = 100;
//  double t[] = {
//    0.        , 0.01010101, 0.02020202, 0.03030303, 0.04040404, 0.05050505,
//    0.06060606, 0.07070707, 0.08080808, 0.09090909, 0.1010101 , 0.11111111,
//    0.12121212, 0.13131313, 0.14141414, 0.15151515, 0.16161616, 0.17171717,
//    0.18181818, 0.19191919, 0.2020202 , 0.21212121, 0.22222222, 0.23232323,
//    0.24242424, 0.25252525, 0.26262626, 0.27272727, 0.28282828, 0.29292929,
//    0.3030303 , 0.31313131, 0.32323232, 0.33333333, 0.34343434, 0.35353535,
//    0.36363636, 0.37373737, 0.38383838, 0.39393939, 0.4040404 , 0.41414141,
//    0.42424242, 0.43434343, 0.44444444, 0.45454545, 0.46464646, 0.47474747,
//    0.48484848, 0.49494949, 0.50505051, 0.51515152, 0.52525253, 0.53535354,
//    0.54545455, 0.55555556, 0.56565657, 0.57575758, 0.58585859, 0.5959596 ,
//    0.60606061, 0.61616162, 0.62626263, 0.63636364, 0.64646465, 0.65656566,
//    0.66666667, 0.67676768, 0.68686869, 0.6969697 , 0.70707071, 0.71717172,
//    0.72727273, 0.73737374, 0.74747475, 0.75757576, 0.76767677, 0.77777778,
//    0.78787879, 0.7979798 , 0.80808081, 0.81818182, 0.82828283, 0.83838384,
//    0.84848485, 0.85858586, 0.86868687, 0.87878788, 0.88888889, 0.8989899 ,
//    0.90909091, 0.91919192, 0.92929293, 0.93939394, 0.94949495, 0.95959596,
//    0.96969697, 0.97979798, 0.98989899, 1.        ,
//  };
//  double R[] = {
//    1.24928187, 1.25037791, 1.25164662, 1.25307592, 1.25465372, 1.25636793,
//    1.25820646, 1.26015722, 1.26220783, 1.26432432, 1.26643709, 1.2684732 ,
//    1.27035972, 1.27204779, 1.27356159, 1.27493914, 1.27621846, 1.27743755,
//    1.27861728, 1.27972182, 1.28070361, 1.28151512, 1.2821169 , 1.28252687,
//    1.28278764, 1.28294196, 1.28303254, 1.28308818, 1.28309865, 1.2830469 ,
//    1.28291589, 1.28268856, 1.28234789, 1.28187859, 1.2812746 , 1.28053274,
//    1.27964986, 1.2786228 , 1.27744867, 1.276135  , 1.27470289, 1.27317436,
//    1.27157144, 1.26991054, 1.26818273, 1.26637194, 1.26446206, 1.26243818,
//    1.26029643, 1.25803897, 1.25566806, 1.2531884 , 1.25063723, 1.24807477,
//    1.24556173, 1.2431587 , 1.24088867, 1.23868208, 1.23645539, 1.23412506,
//    1.23161335, 1.22890566, 1.22602622, 1.22299979, 1.21985393, 1.21665898,
//    1.21351899, 1.21053889, 1.20780586, 1.20524866, 1.20271309, 1.20004418,
//    1.19709218, 1.19385834, 1.1905151 , 1.18724412, 1.18422655, 1.18154187,
//    1.17904582, 1.1765643 , 1.17392753, 1.17107626, 1.16806833, 1.16496727,
//    1.16183612, 1.15872996, 1.15569703, 1.15278535, 1.15004111, 1.14745534,
//    1.14495522, 1.14246444, 1.13992525, 1.13736378, 1.13482978, 1.132373  ,
//    1.13004319, 1.12789008, 1.12596345, 1.12431302,
//  };
//  double Z[] = {
//   -1.08698226, -1.08196971, -1.07700781, -1.07209476, -1.06722871, -1.06240785,
//   -1.05763035, -1.05289438, -1.048198  , -1.0435299 , -1.03886333, -1.0341701 ,
//   -1.02942201, -1.02460028, -1.01971465, -1.01478027, -1.00981228, -1.00482584,
//   -0.99983182, -0.99482694, -0.98980503, -0.98475992, -0.97968661, -0.97458858,
//   -0.96947299, -0.96434696, -0.95921765, -0.95409065, -0.94896713, -0.9438475 ,
//   -0.93873217, -0.93362155, -0.92851607, -0.92341684, -0.91832872, -0.91325773,
//   -0.90820991, -0.90319127, -0.89820775, -0.89326079, -0.88834607, -0.88345886,
//   -0.87859442, -0.87374979, -0.86893001, -0.86414237, -0.85939418, -0.85469247,
//   -0.85004179, -0.84544538, -0.84090641, -0.83642654, -0.83198715, -0.82755527,
//   -0.82309763, -0.81858101, -0.81399271, -0.80937046, -0.80475963, -0.80020557,
//   -0.79575083, -0.79140744, -0.78716867, -0.78302749, -0.77897453, -0.77496393,
//   -0.77092106, -0.76677053, -0.76244948, -0.75800658, -0.75354893, -0.74918416,
//   -0.7450163 , -0.7410448 , -0.73715056, -0.73320813, -0.72909236, -0.72474766,
//   -0.72027146, -0.71578157, -0.71139291, -0.70714603, -0.70300264, -0.69892065,
//   -0.69485832, -0.69077974, -0.68665405, -0.68245053, -0.67813961, -0.67372661,
//   -0.66925716, -0.66477913, -0.66032943, -0.6558955 , -0.65145088, -0.6469691 ,
//   -0.64242369, -0.63778818, -0.63303609, -0.62814097,
//  };
  // Extended plates:
  double t[] = {
    0.        , 0.01010101, 0.02020202, 0.03030303, 0.04040404, 0.05050505,
    0.06060606, 0.07070707, 0.08080808, 0.09090909, 0.1010101 , 0.11111111,
    0.12121212, 0.13131313, 0.14141414, 0.15151515, 0.16161616, 0.17171717,
    0.18181818, 0.19191919, 0.2020202 , 0.21212121, 0.22222222, 0.23232323,
    0.24242424, 0.25252525, 0.26262626, 0.27272727, 0.28282828, 0.29292929,
    0.3030303 , 0.31313131, 0.32323232, 0.33333333, 0.34343434, 0.35353535,
    0.36363636, 0.37373737, 0.38383838, 0.39393939, 0.4040404 , 0.41414141,
    0.42424242, 0.43434343, 0.44444444, 0.45454545, 0.46464646, 0.47474747,
    0.48484848, 0.49494949, 0.50505051, 0.51515152, 0.52525253, 0.53535354,
    0.54545455, 0.55555556, 0.56565657, 0.57575758, 0.58585859, 0.5959596 ,
    0.60606061, 0.61616162, 0.62626263, 0.63636364, 0.64646465, 0.65656566,
    0.66666667, 0.67676768, 0.68686869, 0.6969697 , 0.70707071, 0.71717172,
    0.72727273, 0.73737374, 0.74747475, 0.75757576, 0.76767677, 0.77777778,
    0.78787879, 0.7979798 , 0.80808081, 0.81818182, 0.82828283, 0.83838384,
    0.84848485, 0.85858586, 0.86868687, 0.87878788, 0.88888889, 0.8989899 ,
    0.90909091, 0.91919192, 0.92929293, 0.93939394, 0.94949495, 0.95959596,
    0.96969697, 0.97979798, 0.98989899, 1.        ,
  };
  double R[] = {
    1.23400147, 1.23596982, 1.23692438, 1.23760647, 1.23873964, 1.24046625,
    1.24223112, 1.24357829, 1.24473493, 1.24609806, 1.24768076, 1.24927123,
    1.25072248, 1.25216181, 1.25378298, 1.25571142, 1.25800448, 1.26071586,
    1.26375963, 1.26674464, 1.26923947, 1.270992  , 1.27246039, 1.27427511,
    1.27662659, 1.27883727, 1.28012797, 1.27985427, 1.27885712, 1.2788943 ,
    1.28076395, 1.28304269, 1.28410951, 1.28383   , 1.28313288, 1.28288408,
    1.28294436, 1.28269267, 1.28188045, 1.28045658, 1.27875059, 1.27715353,
    1.27600109, 1.27507859, 1.27385662, 1.27199829, 1.26963499, 1.26696543,
    1.26414096, 1.26124713, 1.2583607 , 1.25547296, 1.25249254, 1.24935156,
    1.24611617, 1.24289388, 1.23972915, 1.23652666, 1.23317353, 1.22961774,
    1.22589572, 1.22205141, 1.21817485, 1.21443592, 1.21101198, 1.20792728,
    1.20483227, 1.20134014, 1.19741473, 1.19333087, 1.18934956, 1.18552628,
    1.18181437, 1.17817041, 1.17456403, 1.17096752, 1.16735816, 1.16372798,
    1.1600751 , 1.15645747, 1.15298346, 1.14971821, 1.14660428, 1.14356418,
    1.1405544 , 1.13757905, 1.1346467 , 1.13179194, 1.12907914, 1.12655402,
    1.12413273, 1.12168382, 1.11916745, 1.11664274, 1.11417392, 1.11183193,
    1.10969647, 1.10784783, 1.10636628, 1.10533213,
  };
  double Z[] = {
    -1.15621288, -1.15003049, -1.14365673, -1.13723457, -1.13090344, -1.12468749,
    -1.11847249, -1.11216563, -1.10581805, -1.09951633, -1.09326264, -1.08700972,
    -1.080725  , -1.07443785, -1.06819323, -1.06203269, -1.05599439, -1.05011605,
    -1.04437945, -1.03864409, -1.03275332, -1.02660506, -1.02031329, -1.01404485,
    -1.00788248, -1.00174305, -0.99552402, -0.98913856, -0.98267291, -0.97632012,
    -0.97016517, -0.96404269, -0.95776537, -0.95131811, -0.94480896, -0.93833879,
    -0.93189358, -0.92543579, -0.91904105, -0.91275869, -0.90652847, -0.90027248,
    -0.89392434, -0.88753274, -0.88121219, -0.87504236, -0.86901791, -0.86312083,
    -0.85731909, -0.85156118, -0.84579663, -0.84002768, -0.83430755, -0.8286761 ,
    -0.82309595, -0.81750594, -0.81188219, -0.80628355, -0.80077929, -0.7954017 ,
    -0.79012948, -0.78493671, -0.77976693, -0.77451069, -0.76905359, -0.76338046,
    -0.75771843, -0.75231866, -0.74720265, -0.74218826, -0.73710256, -0.73190848,
    -0.72663647, -0.72131539, -0.71596723, -0.71061258, -0.70526774, -0.69993643,
    -0.69461799, -0.68927613, -0.68384451, -0.67828401, -0.67263256, -0.66694068,
    -0.66123602, -0.65551422, -0.64976824, -0.6439824 , -0.63813111, -0.63219691,
    -0.6262177 , -0.62025164, -0.61431515, -0.60837962, -0.60241435, -0.59639522,
    -0.59030675, -0.584134  , -0.57786206, -0.57147602,
  };
  // Find indices in t that bound s.
  int idx_tlo, idx_tup;
  if (s < 1e-8) {
    idx_tlo = 0;
    idx_tup = 0;
  }
  else if (fabs(s-1.0) < 1e-8) {
    idx_tlo = npts-1;
    idx_tup = npts-1;
  }
  else {
    for (int i=0; i<npts-1; i++) {
      if (t[i] <= s && s < t[i+1]) {
        idx_tlo = i;
        idx_tup = i+1;
        break;
      }
    }
  }
  // Interpolate the value of R and Z.
  double Dt = t[idx_tup]-t[idx_tlo];
  if (idx_tlo == idx_tup) {
    RZ[0] = R[idx_tlo];
    RZ[1] = Z[idx_tlo];
  }
  else {
    RZ[0] = ((s-t[idx_tlo])/Dt)*R[idx_tup] + ((t[idx_tup]-s)/Dt)*R[idx_tlo];
    RZ[1] = ((s-t[idx_tlo])/Dt)*Z[idx_tup] + ((t[idx_tup]-s)/Dt)*Z[idx_tlo];
  }
}

double rho_psi(double psi, double psi_axis, double psi_sep)
{
  // Normalized radial coordinate.
  return sqrt((psi-psi_axis) / (psi_sep - psi_axis));
}

double psi_rho(double rho, double psi_axis, double psi_sep)
{
  // Poloidal flux given the normalized radial coordinate.
  return pow(rho,2) * (psi_sep - psi_axis) + psi_axis;
}

struct gkyl_gk_block_geom*
create_asdex_lsn_gk_block_geom(void *ctx)
{
  struct gk_asdex_ctx *params = ctx;

  struct gkyl_gk_block_geom *bgeom = gkyl_gk_block_geom_new(params->cdim, params->num_blocks);

  /* Block layout and coordinates.

    z  
    ^  
    |

    |                     +------------------+------------------+
    |                     |                  |                  | 
    |                     |b4                |b3                |
    |                     |inner PF          |lower inner sol   |
    |                     |                  |                  | 
    |                     |%%%%%%%%%%%%%%%%%%|                  |
    |  +------------------+------------------+------------------+
    |  |                 $|                  |$                 |
    |  |                 $|                  |$                 |
    |  | b5              $|                  |$ b2              |
    |  + core            $+                  +$ upper sol       |
    |  |                 $|                  |$                 |
    |  |                 $|                  |$                 |
    |  +------------------+------------------+------------------+
    |                     |%%%%%%%%%%%%%%%%%%|                  |
    |                     |                  |                  | 
    |                     | b0               |b1                |
    |                     | outer PF         |lower outer sol   |
    |                     |                  |                  |
    0                     +------------------+------------------+

       0 -------------------------------------------------------- -> x

    Edges that touch coincide are physically connected unless
    otherwise indicated by a special symbol. Edges with a special
    symbol such as o,x,%, or % are instead connected to the other
    edge with the same symbol. Edges that do not coincide with
    another edge are a physical boundary.
  */  

  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };
  // Copy eqdsk file into efit_inp.
  memcpy(efit_inp.filepath, params->geqdsk_file, sizeof(params->geqdsk_file));

  // Theta limits are actually set internally by the code.
  double theta_min = -1.0, theta_max = 1.0;

  double psi_sep  = params->psi_sep; // Psi at the separatrix.
  double psi_axis = params->psi_axis; // Psi at the magnetic axis.
  double psi_min_core = params->psi_min_core; // Minimum psi the core.
  double psi_max_sol  = params->psi_max_sol ; // Maximum psi the SOL.
  double psi_min_pf   = params->psi_min_pf  ; // Minimum psi the private flux.

  // Number of cells.
  int Npsi_sol        = params->Npsi_sol       ;
  int Npsi_pf         = params->Npsi_pf        ;
  int Npsi_core       = params->Npsi_core      ;
  int Ntheta_divertor = params->Ntheta_divertor;
  int Ntheta_sol      = params->Ntheta_sol     ;

  // Block 0: outer private flux (PF) region.
  gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_min_pf, theta_min },
      .upper = { psi_sep, theta_max },
      .cells = { Npsi_pf, Ntheta_divertor },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_PF_LO_R,
          .rleft = 1.1,
          .rright = 1.7,
          .rmin = 1.1,
          .rmax = 1.7,
          .zmin = -1.3,
          .zmax = -0.9,
          .zmin_left = -1.18,
          .zmax_right = -1.18,
          .plate_spec = true,
          .plate_func_lower = divertor_plate_func_out,
          .plate_func_upper = divertor_plate_func_in,
        }
      },

      .connections[0] = { // x-direction.
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL },  // Physical boundary.
        { .bid = 1, .dir = 0, .edge = GKYL_LOWER_POSITIVE },
      },
      .connections[1] = { // z-direction.
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // Physical boundary.
        { .bid = 4, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // Block 1: lower outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 1, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_sep, theta_min },
      .upper = { psi_max_sol, theta_max },
      .cells = { Npsi_sol, Ntheta_divertor },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO,
          .rclose = 2.5,
          .rleft = 0.8,
          .rright = 2.5,
          .rmin = 0.8,
          .rmax = 2.5,
          .zmin = -1.3,
          .zmax = 1.0,
          .plate_spec = true,
          .plate_func_lower = divertor_plate_func_out,
          .plate_func_upper = divertor_plate_func_in,
        }
      },
      
      .connections[0] = { // x-direction.
        { .bid = 0, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 1, .dir = 0, .edge = GKYL_PHYSICAL }, // Physical boundary.
      },
      .connections[1] = { // z-direction.
        { .bid = 1, .dir = 1, .edge = GKYL_PHYSICAL}, // Physical boundary.
        { .bid = 2, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // Block 2: mid SOL.
  gkyl_gk_block_geom_set_block(bgeom, 2, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_sep, theta_min },
      .upper = { psi_max_sol, theta_max },
      .cells = { Npsi_sol, Ntheta_sol },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID,
          .rclose = 2.5,
          .rleft = 0.8,
          .rright = 2.5,
          .rmin = 0.8,
          .rmax = 2.5,
          .zmin = -1.3,
          .zmax = 1.0,
          .plate_spec = true,
          .plate_func_lower = divertor_plate_func_out,
          .plate_func_upper = divertor_plate_func_in,
        }
      },
      
      .connections[0] = { // x-direction.
        { .bid = 5, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 2, .dir = 0, .edge = GKYL_PHYSICAL }, // Physical boundary.
      },
      .connections[1] = { // z-direction.
        { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 3, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // Block 3: lower inner SOL.
  gkyl_gk_block_geom_set_block(bgeom, 3, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_sep, theta_min },
      .upper = { psi_max_sol, theta_max },
      .cells = { Npsi_sol, Ntheta_divertor },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP,
          .rclose = 2.5,
          .rleft = 0.8,
          .rright = 2.5,
          .rmin = 0.8,
          .rmax = 2.5,
          .zmin = -1.3,
          .zmax = 1.0,
          .plate_spec = true,
          .plate_func_lower = divertor_plate_func_out,
          .plate_func_upper = divertor_plate_func_in,
        }
      },
      
      .connections[0] = { // x-direction.
        { .bid = 4, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 3, .dir = 0, .edge = GKYL_PHYSICAL }, // Physical boundary.
      },
      .connections[1] = { // z-direction.
        { .bid = 2, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 3, .dir = 1, .edge = GKYL_PHYSICAL}, // Physical boundary.
      }
    }
  );

  // Block 4: inner private flux (PF) region.
  gkyl_gk_block_geom_set_block(bgeom, 4, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_min_pf, theta_min },
      .upper = { psi_sep, theta_max },
      .cells = { Npsi_pf, Ntheta_divertor },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_PF_LO_L,
          .rleft = 1.1,
          .rright = 1.7,
          .rmin = 1.1,
          .rmax = 1.7,
          .zmin = -1.3,
          .zmax = -0.9,
          .zmin_left = -1.18,
          .zmax_right = -1.18,
          .plate_spec = true,
          .plate_func_lower = divertor_plate_func_out,
          .plate_func_upper = divertor_plate_func_in,
        }
      },

      .connections[0] = { // x-direction.
        { .bid = 4, .dir = 0, .edge = GKYL_PHYSICAL },  // Physical boundary.
        { .bid = 3, .dir = 0, .edge = GKYL_LOWER_POSITIVE },
      },
      .connections[1] = { // z-direction.
        { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE}, // Physical boundary.
        { .bid = 4, .dir = 1, .edge = GKYL_PHYSICAL},
      }
    }
  );

  // Block 5: core region.
  gkyl_gk_block_geom_set_block(bgeom, 5, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_min_core, theta_min },
      .upper = { psi_sep, theta_max },
      .cells = { Npsi_core, Ntheta_sol },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_GEOMETRY_TOKAMAK_CORE,
          .rclose = 2.0,
          .rleft = 0.8,
          .rright = 2.5,
          .rmin = 0.8,
          .rmax = 2.5,
          .zmin = -1.3,
          .zmax = 1.0,
        }
      },

      .connections[0] = { // x-direction.
        { .bid = 5, .dir = 0, .edge = GKYL_PHYSICAL },  // Physical boundary.
        { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE },
      },
      .connections[1] = { // z-direction.
        { .bid = 5, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 5, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  return bgeom;
}

// Velocity space mappings.
void mapc2p_vel_elc(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double vpar_max_elc = app->vpar_max_elc;
  double mu_max_elc = app->mu_max_elc;

  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_elc*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_elc*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_elc*2.0*pow(cvpar,2);

  // Quadratic map in mu.
  vp[1] = mu_max_elc*pow(cmu,2);
}

void mapc2p_vel_ion(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;

  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_ion*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_ion*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_ion*2.0*pow(cvpar,2);

  // Quadratic map in mu.
  vp[1] = mu_max_ion*pow(cmu,2);
}

double
init_profile(double psi, double f_min, double f_max, void *ctx)
{
  // Profile in D. Michels, et al. Phys. Plasmas 29, 032307 (2022) eqn 17:
  struct gk_asdex_ctx *params = ctx;
  double psi_min = params->psi_min_core;
  double psi_max = params->psi_max_sol;
  double psi_axis = params->psi_axis;
  double psi_sep = params->psi_sep;

  double rho_min = rho_psi(psi_min, psi_axis, psi_sep);
  double rho_max = rho_psi(psi_max, psi_axis, psi_sep);
  double rho = rho_psi(psi, psi_axis, psi_sep);

  double c1 = (f_max-f_min)/2.0;
  double c2 = M_PI/(rho_max-rho_min);
  double c3 = M_PI/2 - c2*rho_min;
  double c4 = (f_max+f_min)/2.0;

  double f = -1.0;
  if (rho <= rho_min)
    f = f_max;
  else if (rho_min < rho && rho <= rho_max)
    f = c1*sin(c2*rho + c3) + c4;
  else
    f = f_min;

  return f;
}

void
init_dens(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_asdex_ctx *params = ctx;
  double psi = xn[0], theta = xn[1];

  // Mimics the SOL profile in D. Michels, et al. Phys. Plasmas 29, 032307
  // (2022), figure 6 experimental.
  double den_min = 0.2e19;
  double den_max = 2.0e19;

  fout[0] = init_profile(psi, den_min, den_max, ctx);
}

void
init_temp_elc(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_asdex_ctx *params = ctx;
  double psi = xn[0], theta = xn[1];

  // Mimics the SOL profile in D. Michels, et al. Phys. Plasmas 29, 032307
  // (2022), figure 8 experimental.
  double eV = GKYL_ELEMENTARY_CHARGE; // Elementary charge.
  double T_min = 17.0*eV;
  double T_max = 300.0*eV;

  fout[0] = init_profile(psi, T_min, T_max, ctx);
}

void
init_temp_ion(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct gk_asdex_ctx *params = ctx;
  double psi = xn[0], theta = xn[1];

  // Mimics the SOL profile in D. Michels, et al. Phys. Plasmas 29, 032307
  // (2022), figure 7 GRILLIX w/ neutrals.
  double eV = GKYL_ELEMENTARY_CHARGE; // Elementary charge.
  double T_min = 17.0*eV;
  double T_max = 300.0*eV;

  fout[0] = init_profile(psi, T_min, T_max, ctx);
}

void
init_upar(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 0.0;
}

void
init_source_dens(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *params = ctx;
  double x = xn[0], z = xn[1];

  double lambda_src = params->lambda_src;
  double psi_src = params->psi_src;
  double ndot_src = params->ndot_src;

  double source_floor = 1e-10;
  if (x < psi_src + 3*lambda_src)
    source_floor = 1e-2;

  double src_prof = exp(-pow(x-psi_src,2)/(2*pow(lambda_src,2)));
  fout[0] = ndot_src * fmax(src_prof, source_floor);
}

void
init_source_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
init_source_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *params = ctx;
  double x = xn[0], z = xn[1];

  double lambda_src = params->lambda_src;
  double psi_src = params->psi_src;
  double Te_src = params->Te_src;
  double eV = GKYL_ELEMENTARY_CHARGE;

  if (x < psi_src + 3*lambda_src)
    fout[0] = Te_src;
  else
    fout[0] = 2.0*eV;
}

void
init_source_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *params = ctx;
  double x = xn[0], z = xn[1];

  double lambda_src = params->lambda_src;
  double psi_src = params->psi_src;
  double Ti_src = params->Ti_src;
  double eV = GKYL_ELEMENTARY_CHARGE;

  if (x < psi_src + 3*lambda_src)
    fout[0] = Ti_src;
  else
    fout[0] = 2.0*eV;
}

void
diffusion_D_func(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;

  fout[0] = 0.3; // Diffusivity [m^2/s].
}

struct gk_asdex_ctx
create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.

  double eps0 = GKYL_EPSILON0; // Permittivity of free space.
  double eV = GKYL_ELEMENTARY_CHARGE; // Elementary charge.
  double mi = 2.014*GKYL_PROTON_MASS; // Ion mass.
  double me = GKYL_ELECTRON_MASS; // Electron mass.
  double qi = eV; // Ion charge.
  double qe = -eV; // Electron charge.

  double Te = 150.0*eV; // Electron temperature.
  double Ti = 150.0*eV; // Ion temperature.
  double B0 = (1.937854e+00+3.937710e+00)/2.0; // B field amplitude.
  double n0 = 1.0e19; // Particle density.

  // Derived parameters.
  double vt_ion = sqrt(Ti/mi);
  double vt_elc = sqrt(Te/me);
  double c_s = sqrt(Te/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;

  // Collision parameters.
  double nu_frac = 1.0;  
  double logLambda_elc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nu_elc = nu_frac*logLambda_elc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*pow(eps0,2)*sqrt(me)*pow(Te,3.0/2.0));

  double logLambda_ion = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nu_ion = nu_frac*logLambda_ion*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*pow(eps0,2)*sqrt(mi)*pow(Ti,3.0/2.0));

  // Location of the numerical equilibrium.
  char geqdsk_file[128] = "gyrokinetic/data/eqdsk/asdex.geqdsk";

  // Position space parameters.
  double num_blocks = 6;
  double R_axis = (1.61640+1.70022)/2.0; // R of the magnetic axis.
  double Z_axis = (-0.0013+0.1001)/2.0; // Z of the magnetic axis.
  double R_sep_OZA = 2.1389435; // Separatrix major at outboard Z axis.
  double R_sep_omp = 2.1334876; // Separatrix major at the OMP.
  double psi_axis = -9.276977e-02; // Psi at the magnetic axis.
  // Get the separatrix psi.
  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };
  // Copy eqdsk file into efit_inp.
  memcpy(efit_inp.filepath, geqdsk_file, sizeof(geqdsk_file));
  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);
  double psi_sep = efit->psisep;
  double Rxpt = efit->Rxpt[0], Zxpt = efit->Zxpt[0];
  gkyl_efit_release(efit);
  // Here rho = sqrt((psi-psi_axis) / (psi_sep - psi_axis)).
  double rho_min_core = 0.9;
  double rho_max_sol = 1.04;
  double psi_min_core = psi_rho(rho_min_core, psi_axis, psi_sep);
  double psi_max_sol = psi_rho(rho_max_sol, psi_axis, psi_sep);
//  double psi_min_pf = 0.1446524024;
  double psi_min_pf = 0.14775;

  // Number of cells.
  int Npsi_sol = 8;
  int Npsi_pf = 2;
  int Npsi_core = 18;
  int Ntheta_divertor = 4;
  int Ntheta_sol = 8;
  int Nvpar = 16; // Number of cells in vpar.
  int Nmu = 8; // Number of cells in mu.

  // Adjust psi_min_core to ensure that dx_core = dx_sol.
  // we need ((psi_sep-shift_fac_core * psi_min_core)/Npsi_core) / ((psi_max_sol-psi_sep)/Npsi_sol) = 1
  double shift_fac_core = (-Npsi_core*psi_max_sol + Npsi_core*psi_sep + Npsi_sol*psi_sep)/(Npsi_sol*psi_min_core);
  psi_min_core *= shift_fac_core;
  double shift_fac_pf = (Npsi_pf*psi_min_core + Npsi_core*psi_sep - Npsi_pf*psi_sep)/(Npsi_core*psi_min_pf);
  psi_min_pf *= shift_fac_pf;
  printf("  shift_fac_core = %9e\n",shift_fac_core);
  printf("  shift_fac_pf = %9e\n",shift_fac_pf);

  double Lx_core = psi_sep - psi_min_core;
  // z location of the X-point on psi=psi_min.
  double z_xpt_psi_sep_lo = -2.8469;
  double z_xpt_psi_sep_up =  2.8486;

  // Source parameters.
  double psi_src = psi_min_core;
  double lambda_src = psi_rho(0.915, psi_axis, psi_sep) - psi_min_core;
  double Lc_src = 67.0; // Connection length in near SOL.
  double n_sep = 0.75e19;
  double Te_sep = 70.0*eV;
  double cs_sep = sqrt(Te_sep/mi);
  double ndot_src = 2.0*n_sep*cs_sep/Lc_src;
  double Te_src = 2.*Te;
  double Ti_src = 2.*Ti;

  // Physical velocity space limits
  double vpar_max_elc = 6.0*vt_elc;
  double mu_max_elc = me*pow(4.0*vt_elc,2)/(2.0*B0);

  double vpar_max_ion = 6.0*vt_ion;
  double mu_max_ion = mi*pow(4.0*vt_ion,2)/(2.0*B0);

  // Computational velocity space limits.
  double vpar_min_ion_c = -1.0/sqrt(2.0);
  double vpar_max_ion_c = 1.0/sqrt(2.0);
  double mu_min_ion_c = 0.;
  double mu_max_ion_c = 1.;
  // Computational velocity space limits.
  double vpar_min_elc_c = -1.0/sqrt(2.0);
  double vpar_max_elc_c = 1.0/sqrt(2.0);
  double mu_min_elc_c = 0.;
  double mu_max_elc_c = 1.;

  // Longest radial chord in the core in m. This is the one touching the x-point, calculated graphically.
  double Lx_core_m_max = sqrt(pow(1.5667-1.4433,2)+pow(-0.5603-(-0.9245),2));
  // Longest radial core in the sol in m. This is the one touching the x-point and the HFS wall, calculated graphically.
  double Lx_sol_m_max = sqrt(pow(1.1511-1.4433,2)+pow(-0.5513-(-0.9245),2));
  // rho_s with 120 eV is about 5.4e-4 m. 

  printf("  X-point @ (R,Z) = (%.9e,%9e)\n",Rxpt,Zxpt);
  printf("  psi_axis = %.13e\n",psi_axis);
  printf("  psi_sep = %.13e\n",psi_sep);
  printf("  psi_min_core = %.13e\n",psi_min_core);
  printf("  psi_max_sol = %.13e\n",psi_max_sol);
  printf("  psi_min_pf = %.13e\n",psi_min_pf);
  printf("  Npsi_sol        = %d\n",Npsi_sol       );
  printf("  Npsi_pf         = %d\n",Npsi_pf        );
  printf("  Npsi_core       = %d\n",Npsi_core      );
  printf("  Ntheta_divertor = %d\n",Ntheta_divertor);
  printf("  Ntheta_sol      = %d\n",Ntheta_sol     );

  double t_end = 1.0e-7;
  double num_frames = 1;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_asdex_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .charge_elc = qe, 
    .charge_ion = qi, 
    .mass_elc = me, 
    .mass_ion = mi,
    .Te = Te, 
    .Ti = Ti, 
    .c_s = c_s, 
    .nu_elc = nu_elc, 
    .nu_ion = nu_ion, 
    .B0 = B0, 
    .n0 = n0, 
    .num_blocks = num_blocks,
    .psi_axis = psi_axis,
    .psi_sep = psi_sep,
    .psi_min_core = psi_min_core,
    .psi_max_sol = psi_max_sol,
    .psi_min_pf = psi_min_pf,
    .Lx_core = Lx_core, 
    .z_xpt_psi_sep_lo = z_xpt_psi_sep_lo,
    .z_xpt_psi_sep_up = z_xpt_psi_sep_up,
    .lambda_src = lambda_src,
    .psi_src = psi_src,
    .ndot_src = ndot_src,
    .Te_src = Te_src,
    .Ti_src = Ti_src,
    // Physical velocity space limits
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
    // Computational velocity space limits
    .vpar_min_elc_c = vpar_min_elc_c,
    .vpar_max_elc_c = vpar_max_elc_c,
    .mu_min_elc_c = mu_min_elc_c,
    .mu_max_elc_c = mu_max_elc_c,
    .vpar_min_ion_c = vpar_min_ion_c,
    .vpar_max_ion_c = vpar_max_ion_c,
    .mu_min_ion_c = mu_min_ion_c,
    .mu_max_ion_c = mu_max_ion_c,
    .Npsi_sol        = Npsi_sol       ,
    .Npsi_pf         = Npsi_pf        ,
    .Npsi_core       = Npsi_core      ,
    .Ntheta_divertor = Ntheta_divertor,
    .Ntheta_sol      = Ntheta_sol     ,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells_v = {Nvpar, Nmu},
    .t_end = t_end, 
    .num_frames = num_frames, 
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  // Copy eqdsk file into ctx.
  memcpy(ctx.geqdsk_file, geqdsk_file, sizeof(geqdsk_file));
  return ctx;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  struct gk_asdex_ctx ctx = create_ctx(); // Context for init functions.
                    
  // Construct block geometry
  struct gkyl_gk_block_geom *bgeom = create_asdex_lsn_gk_block_geom(&ctx);

  int cells_v[ctx.vdim];
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells_v[d]);

  struct gkyl_gyrokinetic_projection elc_ic = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
    .ctx_density = &ctx,
    .density = init_dens,
    .ctx_upar = &ctx,
    .upar = init_upar,
    .ctx_temp = &ctx,
    .temp = init_temp_elc,
  };

  struct gkyl_gyrokinetic_multib_species_pb elc_sol = {
    .polarization_density = ctx.n0,
    .projection = elc_ic,
  };

  // Electrons.
  struct gkyl_gyrokinetic_multib_species_pb elc_blocks[ctx.num_blocks];

  elc_sol.block_id = 0;
  elc_blocks[0] = elc_sol;

  elc_sol.block_id = 1;
  elc_blocks[1] = elc_sol;

  elc_sol.block_id = 2;
  elc_blocks[2] = elc_sol;

  elc_sol.block_id = 3;
  elc_blocks[3] = elc_sol;

  elc_sol.block_id = 4;
  elc_blocks[4] = elc_sol;

  elc_blocks[5] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 5,
    .polarization_density = ctx.n0,

    .projection = elc_ic,

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = init_source_dens,
        .ctx_upar = &ctx,
        .upar = init_source_upar,
        .ctx_temp = &ctx,
        .temp = init_source_temp_elc,
      }, 
      .diagnostics = {
        .num_diag_moments = 1,
        .diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP, },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//        .time_integrated = true,
      },
    },

  };

  struct gkyl_gyrokinetic_bc elc_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
  };

  struct gkyl_gyrokinetic_multib_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .vdim = ctx.vdim,
    .lower = { ctx.vpar_min_elc_c, ctx.mu_min_elc_c},
    .upper = { ctx.vpar_max_elc_c, ctx.mu_max_elc_c},
    .cells = { cells_v[0], cells_v[1] },

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
      .no_by = true,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm.
      .temp_ref = ctx.Te, // Temperature used to calculate coulomb logarithm.
      .bmag_ref = ctx.B0,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
//      .write_diagnostics = true,
    },

    .num_physical_bcs = 10,
    .bcs = elc_phys_bcs,

    .blocks = elc_blocks,
    .duplicate_across_blocks = false,

    .num_diag_moments = 4,
    .diag_moments = { GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_BIMAXWELLIAN },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,
    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },
  };

  // Ions.
  struct gkyl_gyrokinetic_projection ion_ic = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
    .ctx_density = &ctx,
    .density = init_dens,
    .ctx_upar = &ctx,
    .upar = init_upar,
    .ctx_temp = &ctx,
    .temp = init_temp_ion,
  };

  struct gkyl_gyrokinetic_multib_species_pb ion_sol = {
    .polarization_density = ctx.n0,
    .projection = ion_ic,
  };

  struct gkyl_gyrokinetic_multib_species_pb ion_blocks[ctx.num_blocks];
  ion_sol.block_id = 0;
  ion_blocks[0] = ion_sol;

  ion_sol.block_id = 1;
  ion_blocks[1] = ion_sol;

  ion_sol.block_id = 2;
  ion_blocks[2] = ion_sol;

  ion_sol.block_id = 3;
  ion_blocks[3] = ion_sol;

  ion_sol.block_id = 4;
  ion_blocks[4] = ion_sol;

  ion_blocks[5] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 5,
    .polarization_density = ctx.n0,

    .projection = ion_ic,

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = init_source_dens,
        .ctx_upar = &ctx,
        .upar = init_source_upar,
        .ctx_temp = &ctx,
        .temp = init_source_temp_ion,
      }, 
      .diagnostics = {
        .num_diag_moments = 1,
        .diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP, },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//        .time_integrated = true,
      },
    },

  };

  struct gkyl_gyrokinetic_bc ion_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
  };

  struct gkyl_gyrokinetic_multib_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .vdim = ctx.vdim,
    .lower = { ctx.vpar_min_ion_c, ctx.mu_min_ion_c},
    .upper = { ctx.vpar_max_ion_c, ctx.mu_max_ion_c},
    .cells = { cells_v[0], cells_v[1] },

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
      .no_by = true,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm.
      .temp_ref = ctx.Ti, // Temperature used to calculate coulomb logarithm.
      .bmag_ref = ctx.B0,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
//      .write_diagnostics = true,
    },

    .num_physical_bcs = 10,
    .bcs = ion_phys_bcs,

    .blocks = ion_blocks,
    .duplicate_across_blocks = false,

    .num_diag_moments = 4,
    .diag_moments = { GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_BIMAXWELLIAN },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,
    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },
  };

  // Field object.
  struct gkyl_gyrokinetic_multib_field_pb field_blocks[1];
  field_blocks[0] = (struct gkyl_gyrokinetic_multib_field_pb) {
    // No block specific field info for this simulation
    .time_rate_diagnostics = true,
  };

  struct gkyl_gyrokinetic_bc field_phys_bcs[] = {
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 1, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 3, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 4, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 5, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_NEUMANN, .value = {0.0} },
  };

  struct gkyl_gyrokinetic_multib_field field = {
    .blocks = field_blocks, 
    .duplicate_across_blocks = true,

    .num_physical_bcs = ctx.num_blocks,
    .bcs = field_phys_bcs,
  };

  struct gkyl_gyrokinetic_multib *app_inp = gkyl_calloc(1,sizeof(struct gkyl_gyrokinetic_multib));
  app_inp->cdim = ctx.cdim;
  app_inp->poly_order = 1;
  app_inp->basis_type = app_args.basis_type;
  app_inp->use_gpu = app_args.use_gpu;
  app_inp->cfl_frac = 1.0;
  app_inp->gk_block_geom = bgeom;
  app_inp->num_species = 2;
  app_inp->species[0] = elc;
  app_inp->species[1] = ion;
  app_inp->field = field;
  app_inp->comm = comm;


  // Set app output name from the executable name (argv[0]).
  snprintf(app_inp->name, sizeof(app_inp->name), "%s", app_args.app_name);
  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_type = GKYL_GK_MULTIB,
    .multib_app_inp = *app_inp,
    .time_stepping = {
      .t_end = ctx.t_end,
      .num_frames = ctx.num_frames,
      .write_phase_freq = ctx.write_phase_freq,
      .int_diag_calc_num = ctx.int_diag_calc_num,
      .dt_failure_tol = ctx.dt_failure_tol,
      .num_failures_max = ctx.num_failures_max,
      .is_restart = app_args.is_restart,
      .restart_frame = app_args.restart_frame,
      .num_steps = app_args.num_steps,
    },
  };

  gkyl_gyrokinetic_run_simulation(&run_inp);

  gkyl_gk_block_geom_release(bgeom);
  gkyl_gyrokinetic_comms_release(comm);
  gkyl_free(app_inp);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif
  
  return 0;
}
