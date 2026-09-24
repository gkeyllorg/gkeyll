#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include <gkyl_util.h>
#include <gkyl_alloc.h>

#include <mpack.h>
#include <assert.h>

int
gkyl_search_str_int_pair_by_str(const struct gkyl_str_int_pair pairs[], const char *str, int def)
{
  for (int i=0; pairs[i].str != 0; ++i) {
    if (strcmp(pairs[i].str, str) == 0)
      return pairs[i].val;
  }
  return def;  
}

const char *
gkyl_search_str_int_pair_by_int(const struct gkyl_str_int_pair pairs[], int val, const char *def)
{
  for (int i=0; pairs[i].str != 0; ++i) {
    if (pairs[i].val == val)
      return pairs[i].str;
  }
  return def;  
}

int
gkyl_tm_trigger_check_and_bump(struct gkyl_tm_trigger *tmt, double tcurr)
{
  int status = 0;
  if (tcurr >= tmt->tcurr) {
    status = 1;
    tmt->curr += 1;
    tmt->tcurr += tmt->dt;
  }
  return status;
}

void
gkyl_exit(const char* msg)
{
  fprintf(stderr, "Error: %s\n", msg);
  exit(EXIT_FAILURE);
}

int
gkyl_compare_float(float a, float b, float eps)
{
  //if (isnanf(a) || isnanf(b)) return 0;
  
  float absa = fabs(a), absb = fabs(b), diff = fabs(a-b);

  if (a == b) return 1;
  if (a == 0 || b == 0 || (absa+absb < FLT_MIN)) return diff < eps;
  if (absa < eps) return diff < eps;
  if (absb < eps) return diff < eps;
  return diff/fminf(absa+absb, FLT_MAX) < eps;
}

int
gkyl_compare_double(double a, double b, double eps)
{
  if (isnan(a) || isnan(b)) return 0;
  
  double absa = fabs(a), absb = fabs(b), diff = fabs(a-b);
  if (a == b) return 1;
  if (a == 0 || b == 0 || (absa+absb < DBL_MIN)) return diff < eps;
  if (absa < eps) return diff < eps;
  if (absb < eps) return diff < eps;
  return diff/fmin(absa+absb, DBL_MAX) < eps;
}

struct timespec
gkyl_wall_clock(void)
{
  struct timespec tm = { 0 };
#ifdef GKYL_HAVE_CUDA
  cudaDeviceSynchronize();
#endif
  // we were using CLOCK_REALTIME here
  clock_gettime(CLOCK_MONOTONIC, &tm);
  return tm;
}

struct timespec
gkyl_time_diff(struct timespec start, struct timespec end)
{
  struct timespec tm;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    tm.tv_sec = end.tv_sec-start.tv_sec-1;
    tm.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  }
  else {
    tm.tv_sec = end.tv_sec-start.tv_sec;
    tm.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return tm;  
}

double
gkyl_time_diff_now_sec(struct timespec tm)
{
  return gkyl_time_sec(gkyl_time_diff(tm, gkyl_wall_clock()));
}
   
double
gkyl_time_sec(struct timespec tm)
{
  return tm.tv_sec + 1e-9*tm.tv_nsec;
}

double
gkyl_time_now(void)
{
  return gkyl_time_sec( gkyl_wall_clock() );
}

pcg32_random_t
gkyl_pcg32_init(bool nd_seed)
{
  pcg32_random_t rng;
  int rounds = 5;

  if (nd_seed)
    // seed with external entropy -- the time and some program addresses
    // (which will actually be somewhat random on most modern systems).
    pcg32_srandom_r(&rng, time(NULL) ^ (intptr_t)&printf, 
      (intptr_t)&rounds);
  else
    // seed with a fixed constant
    pcg32_srandom_r(&rng, 42u, 54u);

  return rng;
}

uint32_t
gkyl_pcg32_rand_uint32(pcg32_random_t* rng)
{
  return pcg32_random_r(rng);
}

double
gkyl_pcg32_rand_double(pcg32_random_t* rng)
{
  return ldexp(pcg32_random_r(rng), -32);
}

static void
pcg64_srandom_r(pcg64_random_t* rng, uint64_t seed1, uint64_t seed2,
  uint64_t seq1,  uint64_t seq2)
{
  uint64_t mask = ~0ull >> 1;
  // stream for each generators *must* be distinct
  if ((seq1 & mask) == (seq2 & mask)) 
    seq2 = ~seq2;
  pcg32_srandom_r(rng->gen,   seed1, seq1);
  pcg32_srandom_r(rng->gen+1, seed2, seq2);
}

static int _dummy_global = 0; // just to provide address for use in seed

pcg64_random_t
gkyl_pcg64_init(bool nd_seed)
{
  pcg64_random_t rng;
  int rounds = 5;

  if (nd_seed)
    pcg64_srandom_r(&rng,
      time(NULL) ^ (intptr_t)&printf, ~time(NULL) ^ (intptr_t)&pcg32_random_r,
      (intptr_t)&rounds, (intptr_t)&_dummy_global);
  else
    pcg64_srandom_r(&rng, 42u, 42u, 54u, 54u);

  return rng;
}

uint64_t
gkyl_pcg64_rand_uint64(pcg64_random_t* rng)
{
  return ((uint64_t)(pcg32_random_r(rng->gen)) << 32) | pcg32_random_r(rng->gen+1);
}

double
gkyl_pcg64_rand_double(pcg64_random_t* rng)
{
  return ldexp(gkyl_pcg64_rand_uint64(rng), -64);
}

bool
gkyl_check_file_exists(const char *fname)
{
  return access(fname, F_OK) == 0;
}

int64_t
gkyl_file_size(const char *fname)
{
  struct stat st;
  stat(fname, &st);
  return st.st_size;
}

char*
gkyl_load_file(const char *fname, int64_t *sz)
{
  int64_t msz = gkyl_file_size(fname);
  char *buff = gkyl_malloc(msz);
  FILE *fp = fopen(fname, "r");
  int n = fread(buff, msz, 1, fp);
  *sz = msz;
  fclose(fp);
  return buff;
}

bool
gkyl_msgpack_map_elem_has_key(int nvals, const struct gkyl_msgpack_map_elem *elist,
  const char *key)
{
  bool has_key = false;
  for (int i=0; i<nvals; ++i) {
    if (strcmp(key, elist[i].key) == 0) {
      has_key = true;
      break;
    }
  }
  return has_key;
}

struct gkyl_msgpack_map_elem *
gkyl_msgpack_map_elem_clone(int nvals, const struct gkyl_msgpack_map_elem *elist_in)
{
  struct gkyl_msgpack_map_elem *elist_out = gkyl_malloc(nvals*sizeof(struct gkyl_msgpack_map_elem));

  for (int i=0; i<nvals; ++i) {
    // Copy key name.
    size_t key_len = strlen(elist_in[i].key) + 1;
    elist_out[i].key = gkyl_malloc(key_len);
    strcpy(elist_out[i].key, elist_in[i].key);

    // Copy type.
    elist_out[i].elem_type = elist_in[i].elem_type;

    // Copy value.
    switch (elist_in[i].elem_type) {
      case GKYL_MP_BOOL:
        elist_out[i].bval = elist_in[i].bval;
        break;

      case GKYL_MP_UNSIGNED_INT:
        elist_out[i].uval = elist_in[i].uval;
        break;

      case GKYL_MP_INT:
        elist_out[i].ival = elist_in[i].ival;
        break;

      case GKYL_MP_FLOAT:
        elist_out[i].fval = elist_in[i].fval;
        break;

      case GKYL_MP_DOUBLE:
        elist_out[i].dval = elist_in[i].dval;
        break;

      case GKYL_MP_STRING:
        key_len = strlen(elist_in[i].cval) + 1;
        elist_out[i].cval = gkyl_malloc(key_len);
        strcpy(elist_out[i].cval, elist_in[i].cval);
        break;

      default:
        assert(false); // NYI.
        break;
    }
  }

  return elist_out;
}

struct gkyl_msgpack_map_elem *
gkyl_msgpack_map_elem_union(int numlist_union, int *nvals_union,
  const struct gkyl_msgpack_map_elem **elist_union, int *elist_out_len)
{
  int nvals_tot = 0; // Total number of elements.
  for (int j=0; j<numlist_union; ++j)
    nvals_tot += nvals_union[j];

  assert(nvals_tot > 0);

  elist_out_len[0] = nvals_tot;
  struct gkyl_msgpack_map_elem *elist_out = gkyl_malloc(elist_out_len[0]*sizeof(struct gkyl_msgpack_map_elem));

  int eidx = 0;
  for (int j=0; j<numlist_union; ++j) {

    int nvals_curr = nvals_union[j];
    const struct gkyl_msgpack_map_elem *elist_curr = elist_union[j];

    for (int i=0; i<nvals_curr; ++i) {
      // Copy key name.
      size_t slen = strlen(elist_curr[i].key) + 1;
      elist_out[eidx].key = gkyl_malloc(slen);
      strcpy(elist_out[eidx].key, elist_curr[i].key);

      // Copy type.
      elist_out[eidx].elem_type = elist_curr[i].elem_type;

      // Copy value.
      switch (elist_curr[i].elem_type) {
        case GKYL_MP_BOOL:
          elist_out[eidx].bval = elist_curr[i].bval;
          break;

        case GKYL_MP_UNSIGNED_INT:
          elist_out[eidx].uval = elist_curr[i].uval;
          break;

        case GKYL_MP_INT:
          elist_out[eidx].ival = elist_curr[i].ival;
          break;

        case GKYL_MP_FLOAT:
          elist_out[eidx].fval = elist_curr[i].fval;
          break;

        case GKYL_MP_DOUBLE:
          elist_out[eidx].dval = elist_curr[i].dval;
          break;

        case GKYL_MP_STRING:
          slen = strlen(elist_curr[i].cval) + 1;
          elist_out[eidx].cval = gkyl_malloc(slen);
          strcpy(elist_out[eidx].cval, elist_curr[i].cval);
          break;

	default:
	  assert(false); // NYI.
          break;
      }

      eidx++;
    }
  }

  return elist_out;
}

void
gkyl_msgpack_map_elem_set_double(int nvals, struct gkyl_msgpack_map_elem *elist,
  const char *key, double value)
{
  for (int i=0; i<nvals; ++i) {
    if (strcmp(key, elist[i].key) == 0) {
      assert(elist[i].elem_type == GKYL_MP_DOUBLE);
      elist[i].dval = value;
      break;
    }
  }
}

void
gkyl_msgpack_map_elem_set_uint(int nvals, struct gkyl_msgpack_map_elem *elist,
  const char *key, unsigned int value)
{
  for (int i=0; i<nvals; ++i) {
    if (strcmp(key, elist[i].key) == 0) {
      assert(elist[i].elem_type == GKYL_MP_UNSIGNED_INT);
      elist[i].uval = value;
      break;
    }
  }
}

double
gkyl_msgpack_map_elem_get_double(int nvals, struct gkyl_msgpack_map_elem *elist,
  const char *key)
{
  for (int i=0; i<nvals; ++i) {
    if (strcmp(key, elist[i].key) == 0) {
      assert(elist[i].elem_type == GKYL_MP_DOUBLE);
      return elist[i].dval;
    }
  }
  return 0;
}

unsigned int
gkyl_msgpack_map_elem_get_uint(int nvals, struct gkyl_msgpack_map_elem *elist,
  const char *key)
{
  for (int i=0; i<nvals; ++i) {
    if (strcmp(key, elist[i].key) == 0) {
      assert(elist[i].elem_type == GKYL_MP_UNSIGNED_INT);
      return elist[i].uval;
    }
  }
  return 0;
}

char *
gkyl_msgpack_map_elem_get_string(int nvals, struct gkyl_msgpack_map_elem *elist,
  const char *key)
{
  for (int i=0; i<nvals; ++i) {
    if (strcmp(key, elist[i].key) == 0) {
      assert(elist[i].elem_type == GKYL_MP_STRING);
      return elist[i].cval;
    }
  }
  return 0;
}

void
gkyl_msgpack_map_elem_release_string(int nvals, struct gkyl_msgpack_map_elem *elist,
  const char *key)
{
  for (int i=0; i<nvals; ++i) {
    if (strcmp(key, elist[i].key) == 0) {
      assert(elist[i].elem_type == GKYL_MP_STRING);
      MPACK_FREE(elist[i].cval);
      break;
    }
  }
}

void
gkyl_msgpack_map_elem_release(int nvals, struct gkyl_msgpack_map_elem *elist_in)
{
  for (int i=0; i<nvals; ++i) {
    gkyl_free(elist_in[i].key);

    if (elist_in[i].elem_type == GKYL_MP_STRING)
      gkyl_free(elist_in[i].cval);
  }

  gkyl_free(elist_in);
}

struct gkyl_msgpack_data *
gkyl_msgpack_create(int nvals, const struct gkyl_msgpack_map_elem *elist)
{
  struct gkyl_msgpack_data *mdata = gkyl_malloc(sizeof *mdata);
  mdata->meta_sz = 0;
  mdata->meta = 0;

  mpack_writer_t writer;
  mpack_writer_init_growable(&writer, &mdata->meta, &mdata->meta_sz);

  mpack_build_map(&writer);  
  for (int i=0; i<nvals; ++i) {
    mpack_write_cstr(&writer, elist[i].key);

    switch (elist[i].elem_type) {
      case GKYL_MP_BOOL:
        mpack_write_bool(&writer, elist[i].bval);
        break;

      case GKYL_MP_UNSIGNED_INT:
        mpack_write_u64(&writer, elist[i].uval);
        break;

      case GKYL_MP_INT:
        mpack_write_i64(&writer, elist[i].ival);
        break;

      case GKYL_MP_FLOAT:
        mpack_write_float(&writer, elist[i].fval);
        break;

      case GKYL_MP_DOUBLE:
        mpack_write_double(&writer, elist[i].dval);
        break;

      case GKYL_MP_STRING:
        mpack_write_cstr(&writer, elist[i].cval);
        break;

      default:
        assert(false); // NYI.
        break;
    }
  }

  mpack_complete_map(&writer);

  int status = mpack_writer_destroy(&writer);

  if (status != mpack_ok) {
    MPACK_FREE(mdata->meta); // we need to use free here as mpack does its own malloc
    gkyl_free(mdata);
    mdata = 0;
  }

  return mdata;
}

struct gkyl_msgpack_data *
gkyl_msgpack_create_union(int numlist_union, int *nvals_union, const struct gkyl_msgpack_map_elem **elist_union)
{
  struct gkyl_msgpack_data *mdata = gkyl_malloc(sizeof *mdata);
  mdata->meta_sz = 0;
  mdata->meta = 0;

  mpack_writer_t writer;
  mpack_writer_init_growable(&writer, &mdata->meta, &mdata->meta_sz);

  mpack_build_map(&writer);  

  for (int j=0; j<numlist_union; ++j) {

    int nvals = nvals_union[j];
    const struct gkyl_msgpack_map_elem *elist = elist_union[j];

    for (int i=0; i<nvals; ++i) {
      mpack_write_cstr(&writer, elist[i].key);
  
      switch (elist[i].elem_type) {
        case GKYL_MP_BOOL:
          mpack_write_bool(&writer, elist[i].bval);
          break;
  
        case GKYL_MP_UNSIGNED_INT:
          mpack_write_u64(&writer, elist[i].uval);
          break;
  
        case GKYL_MP_INT:
          mpack_write_i64(&writer, elist[i].ival);
          break;
  
        case GKYL_MP_FLOAT:
          mpack_write_float(&writer, elist[i].fval);
          break;
  
        case GKYL_MP_DOUBLE:
          mpack_write_double(&writer, elist[i].dval);
          break;
  
        case GKYL_MP_STRING:
          mpack_write_cstr(&writer, elist[i].cval);
          break;

        default:
          assert(false); // NYI.
          break;
      }
    }
  }

  mpack_complete_map(&writer);

  int status = mpack_writer_destroy(&writer);

  if (status != mpack_ok) {
    MPACK_FREE(mdata->meta); // we need to use free here as mpack does its own malloc
    gkyl_free(mdata);
    mdata = 0;
  }

  return mdata;
}

static void
msgpack_copy_value(mpack_reader_t* r, mpack_writer_t* w)
{
  // Copy the current entry at the top of the reader stack into the writer
  // stack. If it's a map (most cases in Gkeyll so far), loop over the
  // elements.
  mpack_tag_t tag = mpack_read_tag(r);

  switch (mpack_tag_type(&tag)) {
    case mpack_type_nil:
      mpack_write_nil(w);
      break;
    case mpack_type_bool:
      mpack_write_bool(w, mpack_tag_bool_value(&tag));
      break;
    case mpack_type_int:
      mpack_write_int(w, mpack_tag_int_value(&tag));
      break;
    case mpack_type_uint:
      mpack_write_uint(w, mpack_tag_uint_value(&tag));
      break;
    case mpack_type_float:
      mpack_write_float(w, mpack_tag_float_value(&tag));
      break;
    case mpack_type_double:
      mpack_write_double(w, mpack_tag_double_value(&tag));
      break;
    case mpack_type_str: {
      uint32_t len = mpack_tag_str_length(&tag);
      mpack_start_str(w, len);

      char buffer[64];
      uint32_t remaining = len;
      while (remaining > 0) {
        uint32_t chunk = remaining > sizeof(buffer)? sizeof(buffer) : remaining;
        mpack_read_bytes(r, buffer, chunk);
        mpack_write_bytes(w, buffer, chunk);
        remaining -= chunk;
      }

      mpack_done_str(r);
      mpack_finish_str(w);
      break;
    }
    case mpack_type_map: {
      uint32_t count = mpack_tag_map_count(&tag);
      mpack_start_map(w, count);
      for (uint32_t i=0; i<count; i++) {
        msgpack_copy_value(r, w); // Copy key.
        msgpack_copy_value(r, w); // Copy value.
      }
      mpack_done_map(r);
      mpack_finish_map(w);
      break;
    }
    case mpack_type_array: {
      uint32_t count = mpack_tag_array_count(&tag);
      mpack_start_array(w, count);
      for (uint32_t i=0; i<count; i++) {
        msgpack_copy_value(r, w);
      }
      mpack_done_array(r);
      mpack_finish_array(w);
      break;
    }
    default:
      assert(false); // NYI.
      break;
  }
}

struct gkyl_msgpack_data *
gkyl_msgpack_clone(struct gkyl_msgpack_data *mdata_in)
{
  struct gkyl_msgpack_data *mdata_out = gkyl_malloc(sizeof *mdata_out);
  mdata_out->meta_sz = 0;
  mdata_out->meta = 0;

  mpack_reader_t reader;
  mpack_reader_init_data(&reader, mdata_in->meta, mdata_in->meta_sz);
  mpack_writer_t writer;
  mpack_writer_init_growable(&writer, &mdata_out->meta, &mdata_out->meta_sz);

  msgpack_copy_value(&reader, &writer);

  // Check copy was successful.
  if (
      (mpack_reader_destroy(&reader) != mpack_ok || mpack_writer_destroy(&writer) != mpack_ok) ||
      (!(mdata_in->meta_sz == mdata_out->meta_sz && memcmp(mdata_in->meta, mdata_out->meta, mdata_in->meta_sz) == 0))
     ) {
    fprintf(stderr, "gkyl_msgpack_clone: error copying MessagePack.\n");
    MPACK_FREE(mdata_out->meta); // we need to use free here as mpack does its own malloc
    gkyl_free(mdata_out);
    mdata_out = 0;
  }

  return mdata_out;
}

void
gkyl_msgpack_to_map_elem_list(struct gkyl_msgpack_data* mpack_in, int nvals,
  struct gkyl_msgpack_map_elem *elist)
{
  mpack_tree_t tree;
  mpack_tree_init_data(&tree, mpack_in->meta, mpack_in->meta_sz);
  mpack_tree_parse(&tree);
  mpack_node_t root = mpack_tree_root(&tree);

  for (int i=0; i<nvals; ++i) {
    mpack_node_t node = mpack_node_map_cstr(root, elist[i].key);
  
    switch (elist[i].elem_type) {
      case GKYL_MP_BOOL:
        elist[i].bval = mpack_node_bool(node);
        break;
  
      case GKYL_MP_UNSIGNED_INT:
        elist[i].uval = mpack_node_uint(node);
        break;
  
      case GKYL_MP_INT:
        elist[i].ival = mpack_node_int(node);
        break;
  
      case GKYL_MP_FLOAT:
        elist[i].fval = mpack_node_float(node);
        break;
  
      case GKYL_MP_DOUBLE:
        elist[i].dval = mpack_node_double(node);
        break;
  
      case GKYL_MP_STRING:
        elist[i].cval = mpack_node_cstr_alloc(node, mpack_node_strlen(node)+1);
        break;

      default:
        assert(false); // NYI.
        break;
    }
  }

  mpack_tree_destroy(&tree);
}

void
gkyl_msgpack_data_release(struct gkyl_msgpack_data *mdata)
{
  if (!mdata) return;
  if (mdata->meta_sz > 0)
    MPACK_FREE(mdata->meta);
  gkyl_free(mdata);
}

