#include <gkyl_alloc.h>
#include <gkyl_dynvec.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_ref_count.h>
#include <gkyl_util.h>

#include <errno.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

/* Size by which vector grows each time it is reallocated */
static const size_t DYNVEC_ALLOC_SZ = 1024;

struct gkyl_dynvec_tag {
  enum gkyl_elem_type type; // type of data stored in vector
  size_t elemsz, ncomp; // size of elements, number of 'components'

  size_t cloc; // current location to insert data into
  size_t csize; // current number of elements allocated

  size_t esznc; // elemsz*ncomp
  void *data; // pointer to data
  double *tm_mesh; // time stamps
  
  struct gkyl_ref_count ref_count;  
};

static void
dynvec_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dynvec_tag *dv = container_of(ref, struct gkyl_dynvec_tag, ref_count);
  gkyl_free(dv->data);
  gkyl_free(dv->tm_mesh);
  gkyl_free(dv);
}

gkyl_dynvec
gkyl_dynvec_new(enum gkyl_elem_type type, size_t ncomp)
{
  struct gkyl_dynvec_tag *dv = gkyl_malloc(sizeof(struct gkyl_dynvec_tag));
  dv->type = type;
  dv->elemsz = gkyl_elem_type_size[type];
  dv->ncomp = ncomp;
  dv->esznc = dv->elemsz*dv->ncomp;
  dv->csize = DYNVEC_ALLOC_SZ;
  dv->cloc = 0;
  
  dv->data = gkyl_calloc(dv->csize, dv->esznc);
  dv->tm_mesh = gkyl_calloc(dv->csize, sizeof(double));
  
  dv->ref_count = gkyl_ref_count_init(dynvec_free);
  
  return dv;
}

int gkyl_dynvec_elem_type(gkyl_dynvec vec) { return vec->type; }
int gkyl_dynvec_ncomp(gkyl_dynvec vec) { return vec->ncomp; }

void
gkyl_dynvec_reserve_more(gkyl_dynvec dv, size_t rsize)
{
  int n = gkyl_int_div_up(rsize, DYNVEC_ALLOC_SZ);
  dv->csize = n*DYNVEC_ALLOC_SZ + dv->csize;
  dv->data = gkyl_realloc(dv->data, dv->csize*dv->esznc);
  dv->tm_mesh = gkyl_realloc(dv->tm_mesh, dv->csize*sizeof(double));  
}

void
gkyl_dynvec_append(gkyl_dynvec dv, double tm, const void *data)
{
  size_t loc = dv->cloc;
  if (loc >= dv->csize) {
    dv->csize += DYNVEC_ALLOC_SZ;
    dv->data = gkyl_realloc(dv->data, dv->csize*dv->esznc);
    dv->tm_mesh = gkyl_realloc(dv->tm_mesh, dv->csize*sizeof(double));
  }

  // set data
  dv->tm_mesh[loc] = tm;
  memcpy((char*)dv->data + dv->esznc*loc, data, dv->esznc);
  
  dv->cloc += 1;
}

bool
gkyl_dynvec_get(const gkyl_dynvec dv, size_t idx, void *data)
{
  if (idx >= dv->cloc) return false;
  memcpy(data, (char*)dv->data + dv->esznc*idx, dv->esznc);
  return true;  
}

double
gkyl_dynvec_get_tm(const gkyl_dynvec dv, size_t idx)
{
  if (idx >= dv->cloc) return 0.0;
  return dv->tm_mesh[idx];
}

bool
gkyl_dynvec_getlast(const gkyl_dynvec dv, void *data)
{
  size_t loc = dv->cloc;
  if (loc == 0) return false;
  return gkyl_dynvec_get(dv, loc-1, data);
}

double
gkyl_dynvec_getlast_tm(const gkyl_dynvec dv)
{
  size_t loc = dv->cloc;
  return loc == 0 ? 0.0 : dv->tm_mesh[loc-1];
}

size_t
gkyl_dynvec_size(const gkyl_dynvec vec)
{
  return vec->cloc;
}

size_t
gkyl_dynvec_capacity(const gkyl_dynvec vec)
{
  return vec->csize;
}

void
gkyl_dynvec_clear(gkyl_dynvec dv)
{
  dv->cloc = 0;
  dv->csize = DYNVEC_ALLOC_SZ;
  dv->data = gkyl_realloc(dv->data, dv->csize*dv->esznc);
  dv->tm_mesh = gkyl_realloc(dv->tm_mesh, dv->csize*sizeof(double));
}

void
gkyl_dynvec_clear_all_but(gkyl_dynvec dv, size_t num)
{
  if (num>dv->cloc) return;
  
  size_t cloc = dv->cloc;
  dv->csize = DYNVEC_ALLOC_SZ;

  void *data = gkyl_malloc(num*dv->esznc);
  double *tm_mesh = gkyl_malloc(sizeof(double[num]));

  size_t low = num>cloc ? 0 : cloc-num; // lower index to copy from
  size_t ncpy = num>cloc ? cloc : num; // number of elemetns to copy
  dv->cloc = ncpy;

  memcpy(tm_mesh, dv->tm_mesh+low, ncpy*sizeof(double));
  memcpy(data, (char*)dv->data+low*dv->esznc, ncpy*dv->esznc);
  
  dv->data = gkyl_realloc(dv->data, dv->csize*dv->esznc);
  dv->tm_mesh = gkyl_realloc(dv->tm_mesh, dv->csize*sizeof(double));

  memcpy(dv->data, data, ncpy*dv->esznc);
  memcpy(dv->tm_mesh, tm_mesh, ncpy*sizeof(double));

  gkyl_free(data);
  gkyl_free(tm_mesh);
}

gkyl_dynvec
gkyl_dynvec_acquire(const gkyl_dynvec vec)
{
  gkyl_ref_count_inc(&vec->ref_count);
  return (struct gkyl_dynvec_tag*) vec;
}

static int
gkyl_dynvec_write_mode(const gkyl_dynvec vec,
  const char *fname, const char *mode)
{
  const char g0[5] = "gkyl0";

  FILE *fp = 0;
  with_file (fp, fname, mode) {
    fseek(fp, 0, SEEK_END);
    
    // Version 1 header
    fwrite(g0, sizeof(char[5]), 1, fp);
    uint64_t version = 1;
    fwrite(&version, sizeof(uint64_t), 1, fp);
    fwrite(&gkyl_file_type_int[GKYL_DYNVEC_DATA_FILE], sizeof(uint64_t), 1, fp);
    uint64_t meta_size = 0; // THIS WILL CHANGE ONCE METADATA IS EMBEDDED
    fwrite(&meta_size, sizeof(uint64_t), 1, fp);
    
    uint64_t real_type = gkyl_array_data_type[vec->type];
    fwrite(&real_type, sizeof(uint64_t), 1, fp);

    uint64_t esznc = vec->esznc, size = gkyl_dynvec_size(vec);
    fwrite(&esznc, sizeof(uint64_t), 1, fp);
    fwrite(&size, sizeof(uint64_t), 1, fp); 

    fwrite(vec->tm_mesh, sizeof(double)*size, 1, fp);
    fwrite(vec->data, esznc*size, 1, fp);
  }

  return errno;
}

static int
gkyl_dynvec_write_mode_wmeta(const gkyl_dynvec vec, const char *fname,
  const struct gkyl_msgpack_data *meta, const char *mode)
{
  const char g0[5] = "gkyl0";

  FILE *fp = 0;
  with_file (fp, fname, mode) {
    fseek(fp, 0, SEEK_END);
    
    // Version 1 header
    fwrite(g0, sizeof(char[5]), 1, fp);
    uint64_t version = 1;
    fwrite(&version, sizeof(uint64_t), 1, fp);
    fwrite(&gkyl_file_type_int[GKYL_DYNVEC_DATA_FILE], sizeof(uint64_t), 1, fp);
    uint64_t meta_size = meta->meta_sz;
    fwrite(&meta_size, sizeof(uint64_t), 1, fp);
    if (meta_size > 0)
      fwrite(meta->meta, meta_size, 1, fp);
    
    uint64_t real_type = gkyl_array_data_type[vec->type];
    fwrite(&real_type, sizeof(uint64_t), 1, fp);

    uint64_t esznc = vec->esznc, size = gkyl_dynvec_size(vec);
    fwrite(&esznc, sizeof(uint64_t), 1, fp);
    fwrite(&size, sizeof(uint64_t), 1, fp); 

    fwrite(vec->tm_mesh, sizeof(double)*size, 1, fp);
    fwrite(vec->data, esznc*size, 1, fp);
  }

  return errno;
}

int
gkyl_dynvec_write(const gkyl_dynvec vec, const char *fname)
{
  return gkyl_dynvec_write_mode(vec, fname, "w");
}

int
gkyl_dynvec_awrite(const gkyl_dynvec vec, const char *fname)
{
  return gkyl_dynvec_write_mode(vec, fname, "a");
}

int
gkyl_dynvec_write_wmeta(const gkyl_dynvec vec, const char *fname, const struct gkyl_msgpack_data *meta)
{
  return gkyl_dynvec_write_mode_wmeta(vec, fname, meta, "w");
}

// ncomp returned in 'ncomp'
static bool
gkyl_dynvec_read_ncomp_1(FILE *fp, struct gkyl_dynvec_etype_ncomp *enc)
{
  size_t frr;
  // Version 1 header
  char g0[6];
  if (1 != fread(g0, sizeof(char[5]), 1, fp))
    return false;
  g0[5] = '\0'; // add the NULL
  if (strcmp(g0, "gkyl0") != 0)
    return false;

  uint64_t version;
  frr = fread(&version, sizeof(uint64_t), 1, fp);
  if (version != 1)
    return false;

  uint64_t file_type;
  frr = fread(&file_type, sizeof(uint64_t), 1, fp);
  if (file_type != gkyl_file_type_int[GKYL_DYNVEC_DATA_FILE])
    return false;

  uint64_t meta_size;
  frr = fread(&meta_size, sizeof(uint64_t), 1, fp);

  // read ahead by specified bytes: meta-data is not read in this
  // method
  fseek(fp, meta_size, SEEK_CUR);

  uint64_t real_code = 0;
  if (1 != fread(&real_code, sizeof(uint64_t), 1, fp))
    return false;
  enc->type = gkyl_array_code_to_data_type[real_code];

  uint64_t esznc;
  if (1 != fread(&esznc, sizeof(uint64_t), 1, fp))
    return false;

  int real_type = gkyl_array_code_to_data_type[real_code];
  enc->ncomp = esznc/gkyl_elem_type_size[real_type];
  
  return true;
}


struct gkyl_dynvec_etype_ncomp
gkyl_dynvec_read_ncomp(const char *fname)
{
  struct gkyl_dynvec_etype_ncomp enc = {
    .type = GKYL_DOUBLE,
    .ncomp = 0
  };
  FILE *fp = 0;
  with_file(fp, fname, "r")
    gkyl_dynvec_read_ncomp_1(fp, &enc);
  return enc;
}

static bool
gkyl_dynvec_read_1(gkyl_dynvec vec, FILE *fp) {
  size_t frr;
  // Version 1 header
  char g0[6];
  if (1 != fread(g0, sizeof(char[5]), 1, fp))
    return false;
  g0[5] = '\0'; // add the NULL
  if (strcmp(g0, "gkyl0") != 0)
    return false;

  uint64_t version;
  frr = fread(&version, sizeof(uint64_t), 1, fp);
  if (version != 1)
    return false;

  uint64_t file_type;
  frr = fread(&file_type, sizeof(uint64_t), 1, fp);
  if (file_type != gkyl_file_type_int[GKYL_DYNVEC_DATA_FILE])
    return false;

  uint64_t meta_size;
  frr = fread(&meta_size, sizeof(uint64_t), 1, fp);

  // read ahead by specified bytes: meta-data is not read in this
  // method
  fseek(fp, meta_size, SEEK_CUR);

  uint64_t real_type = 0;
  if (1 != fread(&real_type, sizeof(uint64_t), 1, fp))
    return false;
  if (real_type != gkyl_array_data_type[vec->type])
    return false;

  uint64_t esznc, size;
  if (1 != fread(&esznc, sizeof(uint64_t), 1, fp))
    return false;
  if (vec->esznc != esznc)
    return false;

  if (1 != fread(&size, sizeof(uint64_t), 1, fp))
    return false;

  // resize vector to allow storing new data
  gkyl_dynvec_reserve_more(vec, size);

  // read time-mesh data
  frr = fread(&vec->tm_mesh[vec->cloc], sizeof(double[size]), 1, fp);
  // read dynvec data
  frr = fread((char *)vec->data + vec->esznc * vec->cloc, size * esznc, 1, fp);

  // bump location so further inserts occurs after newly read data
  vec->cloc = vec->cloc + size;

  return true;
}

bool
gkyl_dynvec_read(gkyl_dynvec vec, const char *fname)
{
  bool status = false;
  FILE *fp = fopen(fname, "r");

  // keep reading till we have no more datasets
  while (1) {
    status = gkyl_dynvec_read_1(vec, fp);
    fpos_t curr_pos;
    fgetpos(fp, &curr_pos);
    
    char g0[6];
    if (1 != fread(g0, sizeof(char[5]), 1, fp))
      break;
    fsetpos(fp, &curr_pos);
  }
  fclose(fp);
  
  return status;
}

void
gkyl_dynvec_to_array(const gkyl_dynvec vec, struct gkyl_array *tm_mesh,
  struct gkyl_array *dyndata)
{
  int nv = gkyl_dynvec_size(vec);
  for (int i=0; i<nv; ++i) {
    double *tmm = gkyl_array_fetch(tm_mesh, i);
    tmm[0] = gkyl_dynvec_get_tm(vec, i);

    void *dd = gkyl_array_fetch(dyndata, i);
    gkyl_dynvec_get(vec, i, dd);
  }
}

void
gkyl_dynvec_release(gkyl_dynvec vec)
{
  if (vec)
    gkyl_ref_count_dec(&vec->ref_count);
}
