#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_dg_iz.h>
#include <gkyl_dg_recomb.h>
#include <read_adas.h>

void
read_adas_field_iz(enum gkyl_ion_type type_ion, struct adas_field *data) {
  char fname[4000];
  const char *fmt = "%s/%s";
  if (type_ion == GKYL_ION_H) {
    data->NT = 29, data->NN = 24;
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/ioniz_h.npy");
    data->logData = fopen(fname,"rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logT_h.npy");
    data->logT = fopen(fname, "rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logN_h.npy");
    data->logN = fopen(fname, "rb");    
    data->Zmax = 1;
    data->Eiz[0] = 13.6; 
  }
  else if (type_ion == GKYL_ION_HE) {
    data->NT = 30;
    data->NN = 24;
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/ioniz_he.npy");
    data->logData = fopen(fname,"rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logT_he.npy");
    data->logT = fopen(fname, "rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logN_he.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 2;
    data->Eiz[0] = 24.6;
    data->Eiz[1] = 54.4;
  }
  else if (type_ion == GKYL_ION_LI) {
    data->NT = 25;
    data->NN = 16;
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/ioniz_li.npy");
    data->logData = fopen(fname,"rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logT_li.npy");
    data->logT = fopen(fname, "rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logN_li.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 3;
    data->Eiz[0] = 5.4;
    data->Eiz[1] = 75.6;
    data->Eiz[2] = 122.4;
  }
  else if (type_ion == GKYL_ION_BE) {
    data->NT = 25;
    data->NN = 24;
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/ioniz_be.npy");
    data->logData = fopen(fname,"rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logT_be.npy");
    data->logT = fopen(fname, "rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logN_be.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 4;
    data->Eiz[0] = 9.3;
    data->Eiz[1] = 18.2;
    data->Eiz[3] = 153.9;
    data->Eiz[4] = 217.7;
  }
  else if (type_ion == GKYL_ION_B) {
    data->NT = 48;
    data->NN = 26;
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/ioniz_b.npy");
    data->logData = fopen(fname,"rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logT_b.npy");
    data->logT = fopen(fname, "rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logN_b.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 5;
    data->Eiz[0] = 8.3;
    data->Eiz[1] = 25.2;
    data->Eiz[2] = 37.9;
    data->Eiz[3] = 259.4;
    data->Eiz[5] = 340.2;
  }
  else if (type_ion == GKYL_ION_C) {
    data->NT = 30;
    data->NN = 24;
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/ioniz_c.npy");
    data->logData = fopen(fname,"rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logT_c.npy");
    data->logT = fopen(fname, "rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logN_c.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 6;
    data->Eiz[0] = 11.3;
    data->Eiz[1] = 24.4;
    data->Eiz[2] = 47.9;
    data->Eiz[3] = 64.5;
    data->Eiz[4] = 392.1;
    data->Eiz[5] = 490.0;
  }
  else if (type_ion == GKYL_ION_N) {
    data->NT = 30;
    data->NN = 24;
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/ioniz_n.npy");
    data->logData = fopen(fname,"rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logT_n.npy");
    data->logT = fopen(fname, "rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logN_n.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 7;
    data->Eiz[0] = 14.5;
    data->Eiz[1] = 29.6;
    data->Eiz[2] = 47.5;
    data->Eiz[3] = 77.5;
    data->Eiz[4] = 97.9;
    data->Eiz[5] = 552.1;
    data->Eiz[6] = 667.0;
  }
  else if (type_ion == GKYL_ION_O) {
    data->NT = 30;
    data->NN = 24;
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/ioniz_o.npy");
    data->logData = fopen(fname,"rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logT_o.npy");
    data->logT = fopen(fname, "rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logN_o.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 8;
    data->Eiz[0] = 13.6;
    data->Eiz[1] = 35.1;
    data->Eiz[2] = 54.9;
    data->Eiz[3] = 77.4;
    data->Eiz[4] = 113.9;
    data->Eiz[5] = 138.1;
    data->Eiz[6] = 739.3;
    data->Eiz[7] = 871.4;
  }
  else if (type_ion == GKYL_ION_AR) {
    data->NT = 48;
    data->NN = 26;
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/ioniz_ar.npy");
    data->logData = fopen(fname,"rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logT_ar.npy");
    data->logT = fopen(fname, "rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logN_ar.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 18;
    data->Eiz[0] = 15.8;
    data->Eiz[1] = 27.6;
    data->Eiz[2] = 40.9;
    data->Eiz[3] = 52.3;
    data->Eiz[4] = 75.;
    data->Eiz[5] = 91.;
    data->Eiz[6] = 124.3;
    data->Eiz[7] = 143.5;
    data->Eiz[8] = 422.4;
    data->Eiz[9] = 478.7;
    data->Eiz[10] = 539.;
    data->Eiz[11] = 618.3;
    data->Eiz[12] = 686.1;
    data->Eiz[13] = 755.7;
    data->Eiz[14] = 854.8;
    data->Eiz[15] = 918.;
    data->Eiz[16] = 4120.7;
    data->Eiz[17] = 4426.2;
  }
  else
    fprintf(stderr, "Incorrect ion type for ionization.");
}
 
void
read_adas_field_recomb(enum gkyl_ion_type type_ion, struct adas_field *data) {
  char fname[4000];
  const char *fmt = "%s/%s";
  if (type_ion == GKYL_ION_H) {
    data->NT = 29, data->NN = 24;
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/recomb_h.npy");
    data->logData = fopen(fname,"rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logT_h.npy");
    data->logT = fopen(fname, "rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logN_h.npy");
    data->logN = fopen(fname, "rb");    
    data->Zmax = 1;
   }
  else if (type_ion == GKYL_ION_HE) {
    data->NT = 30;
    data->NN = 24;
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/recomb_he.npy");
    data->logData = fopen(fname,"rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logT_he.npy");
    data->logT = fopen(fname, "rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logN_he.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 2;
   }
  else if (type_ion == GKYL_ION_LI) {
    data->NT = 25;
    data->NN = 16;
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/recomb_li.npy");
    data->logData = fopen(fname,"rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logT_li.npy");
    data->logT = fopen(fname, "rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logN_li.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 3;
  }
  else if (type_ion == GKYL_ION_BE) {
    data->NT = 25;
    data->NN = 24;
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/recomb_be.npy");
    data->logData = fopen(fname,"rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logT_be.npy");
    data->logT = fopen(fname, "rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logN_be.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 4;
  }
  else if (type_ion == GKYL_ION_B) {
    data->NT = 48;
    data->NN = 26;
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/recomb_b.npy");
    data->logData = fopen(fname,"rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logT_b.npy");
    data->logT = fopen(fname, "rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logN_b.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 5;
  }
  else if (type_ion == GKYL_ION_C) {
    data->NT = 30;
    data->NN = 24;
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/recomb_c.npy");
    data->logData = fopen(fname,"rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logT_c.npy");
    data->logT = fopen(fname, "rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logN_c.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 6;
  }
  else if (type_ion == GKYL_ION_N) {
    data->NT = 30;
    data->NN = 24;
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/recomb_n.npy");
    data->logData = fopen(fname,"rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logT_n.npy");
    data->logT = fopen(fname, "rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logN_n.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 7;
  }
  else if (type_ion == GKYL_ION_O) {
    data->NT = 30;
    data->NN = 24;
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/recomb_o.npy");
    data->logData = fopen(fname,"rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logT_o.npy");
    data->logT = fopen(fname, "rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logN_o.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 8;
  }
  else if (type_ion == GKYL_ION_AR) {
    data->NT = 48;
    data->NN = 26;
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/recomb_ar.npy");
    data->logData = fopen(fname,"rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logT_ar.npy");
    data->logT = fopen(fname, "rb");
    snprintf(fname, sizeof fname, fmt, GKYL_SHARE_DIR, "adas/logN_ar.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 18;
  }
  else
    fprintf(stderr, "Incorrect ion type for recombination.");
}
