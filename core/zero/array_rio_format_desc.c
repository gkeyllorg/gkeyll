#include <gkyl_array_rio_format_desc.h>

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <gkyl_util.h>

size_t
gkyl_base_hdr_size(size_t meta_sz)
{
  size_t sz = 0;
  // magic string
  sz += 5; // "gkyl0"
  // version
  sz += sizeof(uint64_t);
  // file_type
  sz += sizeof(uint64_t);
  // metadata
  sz += sizeof(uint64_t) + meta_sz;

  return sz;
}

size_t
gkyl_file_type_1_partial_hrd_size(int ndim)
{
  size_t sz = 0;
  // real_type
  sz += sizeof(uint64_t);
  // ndim
  sz += sizeof(uint64_t);
  // cells
  sz += sizeof(uint64_t[ndim]);
  // lower, upper
  sz += sizeof(double[2*ndim]);
  return sz;
}

size_t
gkyl_file_type_1_hrd_size(int ndim)
{
  size_t sz = gkyl_file_type_1_partial_hrd_size(ndim);
  // esznc
  sz += sizeof(uint64_t);
  // size (total number of cells)
  sz += sizeof(uint64_t);
  return sz;
}

size_t
gkyl_file_type_2_hrd_size(void)
{
  size_t sz = 0;
  // real_type
  sz += sizeof(uint64_t);
  // esznc
  sz += sizeof(uint64_t);
  // size (total number of cells)
  sz += sizeof(uint64_t);
  return sz;
}

size_t
gkyl_file_type_3_hrd_size(int ndim)
{
  size_t sz = gkyl_file_type_1_hrd_size(ndim);
  sz += sizeof(uint64_t); // nrange
  return sz;
}

size_t
gkyl_file_type_3_range_hrd_size(int ndim)
{
  size_t sz = 0;
  // loidx and upidx
  sz += sizeof(uint64_t[2*ndim]);
  sz += sizeof(uint64_t);
  return sz;
}

int
gkyl_get_gkyl_file_type(const char *fname)
{
  int file_type = -1;
  FILE *fp = 0;

  with_file(fp, fname, "r") {
    size_t frr;
    char g0[6];
    frr = fread(g0, sizeof(char[5]), 1, fp); // no trailing '\0'
    g0[5] = '\0';                            // add the NULL
    if (strcmp(g0, "gkyl0") != 0) {
      file_type = -1;
      goto finish_with_file;
    }
  
    uint64_t version;
    frr = fread(&version, sizeof(uint64_t), 1, fp);
    if (version != 1) {
      file_type = -1;
      goto finish_with_file;
    }
    
    uint64_t file_type_u64;
    frr = fread(&file_type_u64, sizeof(uint64_t), 1, fp);
    if (1 != frr) {
      file_type = -1;
      goto finish_with_file;      
    }

    file_type = file_type_u64;
    
    finish_with_file:
    ;
  }
  return file_type;
}
