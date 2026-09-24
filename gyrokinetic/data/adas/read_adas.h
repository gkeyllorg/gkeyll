#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

typedef struct adas_field {
  FILE *logData;
  FILE *logT;
  FILE *logN;
  long NT;
  long NN; 
  int Zmax;
  struct gkyl_array fld;
  double Eiz[GKYL_MAX_CHARGE_STATE];
} adas_field;

// Functions to extract ADAS data and project onto DG data
static inline void
array_from_numpy(FILE *fp, long sz, int Zmax, int charge_state, struct gkyl_array *arr)
{
  int zi = charge_state;
  double array[Zmax][sz];
  long res_sz = fread(array, 1, sizeof(double[Zmax][sz]), fp);
		   
  for (int i=0; i<sz; ++i) {
    double *arr_d = (double*) gkyl_array_fetch(arr, i);
    arr_d[0] = array[zi][i]; 
  }
}

static inline void
minmax_from_numpy(FILE *fp, long sz, double minmax[2])
{
  double array[sz];
  long res_sz = fread(array, 1, sizeof(double[sz]), fp);
  double min = array[0];
  double max = array[sz-1];
  minmax[0] = min;
  minmax[1] = max;
}

void
read_adas_field_iz(enum gkyl_ion_type type_ion, struct adas_field *data);
 
void
read_adas_field_recomb(enum gkyl_ion_type type_ion, struct adas_field *data);
