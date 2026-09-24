#include <assert.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_format_desc.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_elem_type_priv.h>

// Error message strings
static const char *array_rio_status_msg[] = {
  [GKYL_ARRAY_RIO_SUCCESS] = "Success",
  [GKYL_ARRAY_RIO_BAD_VERSION] = "Incorrect header version",
  [GKYL_ARRAY_RIO_FOPEN_FAILED] = "File open failed",
  [GKYL_ARRAY_RIO_FREAD_FAILED] = "Data read failed",
  [GKYL_ARRAY_RIO_DATA_MISMATCH] = "Data mismatch",
  [GKYL_ARRAY_RIO_META_FAILED] = "Metadata output failed"
};

const char*
gkyl_array_rio_status_msg(enum gkyl_array_rio_status status)
{
  return array_rio_status_msg[status];
}

static void
sub_array_write_priv(const struct gkyl_range *range,
  const struct gkyl_array *arr, FILE *fp)
{
#define _F(loc) gkyl_array_cfetch(arr, loc)

  // construct skip iterator to allow writing (potentially) in chunks
  // rather than element by element or requiring a copy of data
  struct gkyl_range_skip_iter skip;
  gkyl_range_skip_iter_init(&skip, range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &skip.range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&skip.range, iter.idx);
    fwrite(_F(start), arr->esznc*skip.delta, 1, fp);
  }
#undef _F
}

int
gkyl_header_meta_write_fp(const struct gkyl_array_header_info *hdr, FILE *fp)
{
  const char g0[5] = "gkyl0";

  // Version 1 header
  fwrite(g0, sizeof(char[5]), 1, fp);
  uint64_t version = 1;
  fwrite(&version, sizeof(uint64_t), 1, fp);
  fwrite(&hdr->file_type, sizeof(uint64_t), 1, fp);
  uint64_t meta_size = hdr->meta_size;
  fwrite(&meta_size, sizeof(uint64_t), 1, fp);
  if (meta_size > 0)
    fwrite(hdr->meta, meta_size, 1, fp);

  return GKYL_ARRAY_RIO_SUCCESS;
}

int
gkyl_grid_sub_array_header_write_fp(const struct gkyl_rect_grid *grid,
  const struct gkyl_array_header_info *hdr, FILE *fp)
{
  gkyl_header_meta_write_fp(hdr, fp);  
  
  // Version 0 format is used for rest of the header
  uint64_t real_type = gkyl_array_data_type[hdr->etype];
  fwrite(&real_type, sizeof(uint64_t), 1, fp);
  gkyl_rect_grid_write(grid, 0, fp);

  fwrite(&hdr->esznc, sizeof(uint64_t), 1, fp);
  fwrite(&hdr->tot_cells, sizeof(uint64_t), 1, fp);

  return GKYL_ARRAY_RIO_SUCCESS;
}

int
gkyl_header_meta_read_fp(struct gkyl_array_header_info *hdr, FILE *fp)
{
  size_t frr;
  hdr->meta_size = 0;

  char g0[6];
  frr = fread(g0, sizeof(char[5]), 1, fp); // no trailing '\0'
  g0[5] = '\0';                            // add the NULL
  if (strcmp(g0, "gkyl0") != 0)
    return GKYL_ARRAY_RIO_BAD_VERSION;
  
  uint64_t version;
  frr = fread(&version, sizeof(uint64_t), 1, fp);
  if (version != 1)
    return GKYL_ARRAY_RIO_BAD_VERSION;

  uint64_t file_type;
  frr = fread(&file_type, sizeof(uint64_t), 1, fp);

  uint64_t meta_size;
  frr = fread(&meta_size, sizeof(uint64_t), 1, fp);
  if (1 != frr)
    return GKYL_ARRAY_RIO_FREAD_FAILED;

  hdr->meta = 0;
  if (meta_size > 0) {
    hdr->meta = gkyl_malloc(meta_size);
    if (1 != fread(hdr->meta, meta_size, 1, fp)) {
      gkyl_free(hdr->meta);
      return GKYL_ARRAY_RIO_FREAD_FAILED;
    }
  }

  hdr->file_type = file_type;
  hdr->esznc = 0;
  hdr->tot_cells = 0;
  hdr->meta_size = meta_size;

  return GKYL_ARRAY_RIO_SUCCESS;  
}

static int
grid_sub_array_header_read_fp(struct gkyl_rect_grid *grid,
  struct gkyl_array_header_info *hdr, bool read_meta, FILE *fp)
{
  size_t frr;
  hdr->meta_size = 0;

  char g0[6];
  frr = fread(g0, sizeof(char[5]), 1, fp); // no trailing '\0'
  g0[5] = '\0';                            // add the NULL
  if (strcmp(g0, "gkyl0") != 0)
    return GKYL_ARRAY_RIO_BAD_VERSION;
  
  uint64_t version;
  frr = fread(&version, sizeof(uint64_t), 1, fp);
  if (version != 1)
    return GKYL_ARRAY_RIO_BAD_VERSION;

  uint64_t file_type;
  frr = fread(&file_type, sizeof(uint64_t), 1, fp);
  if (1 != frr)
    return GKYL_ARRAY_RIO_FREAD_FAILED;

  uint64_t meta_size;
  frr = fread(&meta_size, sizeof(uint64_t), 1, fp);
  if (1 != frr)
    return GKYL_ARRAY_RIO_FREAD_FAILED;

  if (meta_size > 0) {
    if (read_meta) {
      hdr->meta = gkyl_malloc(meta_size);
      if (1 != fread(hdr->meta, meta_size, 1, fp)) {
        gkyl_free(hdr->meta);
        return GKYL_ARRAY_RIO_FREAD_FAILED;
      }
    }
    else {
      fseek(fp, meta_size, SEEK_CUR);
    }
  }

  uint64_t real_type = 0;
  if (1 != fread(&real_type, sizeof(uint64_t), 1, fp))
    return GKYL_ARRAY_RIO_FREAD_FAILED;

  gkyl_rect_grid_read(grid, fp);

  uint64_t esznc = 0;
  if (1 != fread(&esznc, sizeof(uint64_t), 1, fp))
    return GKYL_ARRAY_RIO_FREAD_FAILED;;

  uint64_t tot_cells = 0;
  if (1 != fread(&tot_cells, sizeof(uint64_t), 1, fp))
    return GKYL_ARRAY_RIO_FREAD_FAILED;;

  uint64_t nrange = 1;
  if (file_type == gkyl_file_type_int[GKYL_MULTI_RANGE_DATA_FILE])
    if (1 != fread(&nrange, sizeof(uint64_t), 1, fp))
      return GKYL_ARRAY_RIO_FREAD_FAILED;

  hdr->file_type = file_type;
  hdr->etype = gkyl_array_code_to_data_type[real_type];
  hdr->esznc = esznc;
  hdr->tot_cells = tot_cells;
  hdr->meta_size = meta_size;
  hdr->nrange = nrange;

  return GKYL_ARRAY_RIO_SUCCESS;
}

int
gkyl_grid_sub_array_header_read_fp(struct gkyl_rect_grid *grid,
  struct gkyl_array_header_info *hdr, FILE *fp)
{
  return grid_sub_array_header_read_fp(grid, hdr, true, fp);
}

void
gkyl_grid_sub_array_header_release(struct gkyl_array_header_info *hdr)
{
  if (hdr->meta_size>0) {
    gkyl_free(hdr->meta);
    hdr->meta_size = 0;
  }
}

enum gkyl_array_rio_status
gkyl_grid_sub_array_header_read(struct gkyl_rect_grid *grid,
  struct gkyl_array_header_info *hdr, const char *fname)
{
  enum gkyl_array_rio_status status = GKYL_ARRAY_RIO_FOPEN_FAILED;
  FILE *fp = 0;
  with_file(fp, fname, "r") {
    status = gkyl_grid_sub_array_header_read_fp(grid, hdr, fp);
  }
  return status;
}

void
gkyl_array_header_info_release(struct gkyl_array_header_info *info)
{
  if (info->meta_size > 0)
    gkyl_free(info->meta);
}

enum gkyl_array_rio_status
gkyl_grid_sub_array_write(const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  const struct gkyl_msgpack_data *meta,
  const struct gkyl_array *arr, const char *fname)
{
  enum gkyl_array_rio_status status = GKYL_ARRAY_RIO_FOPEN_FAILED;
  FILE *fp = 0;
  int err;
  with_file (fp, fname, "w") {
    
    status = gkyl_grid_sub_array_header_write_fp(grid,
      &(struct gkyl_array_header_info) {
        .file_type = gkyl_file_type_int[GKYL_FIELD_DATA_FILE],
        .etype = arr->type,
        .esznc = arr->esznc,
        .tot_cells = range->volume,
        .meta_size = meta ? meta->meta_sz : 0,
        .meta = meta ? meta->meta : 0 
      },
      fp
    );

    if (status == GKYL_ARRAY_RIO_SUCCESS)
      sub_array_write_priv(range, arr, fp);
  }
  return status;
}

static enum gkyl_array_rio_status
grid_sub_array_read_ft_1(const struct gkyl_rect_grid *grid,
  struct gkyl_array_header_info *hdr, const struct gkyl_range *range,
  struct gkyl_array *arr, FILE *fp)
{
  size_t loc = gkyl_base_hdr_size(hdr->meta_size)
    + gkyl_file_type_1_hrd_size(grid->ndim);
  fseek(fp, loc, SEEK_SET);

  struct gkyl_range blk_rng;
  gkyl_range_init_from_shape1(&blk_rng, grid->ndim, grid->cells);

  struct gkyl_range inter;
  int not_empty = gkyl_range_intersect(&inter, &blk_rng, range);

  if (not_empty) {
    uint64_t sz = hdr->tot_cells;
    gkyl_mem_buff buff = gkyl_mem_buff_new(sz*hdr->esznc);

    if (1 != fread(gkyl_mem_buff_data(buff), sz*hdr->esznc, 1, fp)) {
      gkyl_mem_buff_release(buff);
      return GKYL_ARRAY_RIO_FREAD_FAILED;
    }
    
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &inter);
    while (gkyl_range_iter_next(&iter)) {
      
      char *out = gkyl_array_fetch(arr, gkyl_range_idx(range, iter.idx));
      const char *inp = gkyl_mem_buff_data(buff) + hdr->esznc*gkyl_range_idx(&blk_rng, iter.idx);
      memcpy(out, inp, hdr->esznc);
    }    

    gkyl_mem_buff_release(buff);
  }
    
  return GKYL_ARRAY_RIO_SUCCESS;
}

static enum gkyl_array_rio_status
grid_sub_array_read_ft_3(const struct gkyl_rect_grid *grid,
  struct gkyl_array_header_info *hdr, const struct gkyl_range *range,
  struct gkyl_array *arr, FILE *fp)
{
  size_t rng_sz = gkyl_file_type_3_range_hrd_size(grid->ndim);
  size_t loc = gkyl_base_hdr_size(hdr->meta_size) + gkyl_file_type_3_hrd_size(grid->ndim);

  gkyl_mem_buff buff = gkyl_mem_buff_new(10); // will be reallocated

  for (int r=0; r<hdr->nrange; ++r) {
    uint64_t sz, loidx[GKYL_MAX_DIM], upidx[GKYL_MAX_DIM];
    fseek(fp, loc, SEEK_SET);

    // read lower, upper indices and number of elements stored
    if (1 != fread(loidx, sizeof(uint64_t[grid->ndim]), 1, fp))
      return GKYL_ARRAY_RIO_FREAD_FAILED;
    if (1 != fread(upidx, sizeof(uint64_t[grid->ndim]), 1, fp))
      return GKYL_ARRAY_RIO_FREAD_FAILED;
    if (1 != fread(&sz, sizeof(uint64_t), 1, fp))
      return GKYL_ARRAY_RIO_FREAD_FAILED;

    // construct range of indices corresponding to data in block
    int loidx_i[GKYL_MAX_DIM]= { 0 } , upidx_i[GKYL_MAX_DIM] = { 0 };
    for (int d=0; d<grid->ndim; ++d) {
      loidx_i[d] = loidx[d];
      upidx_i[d] = upidx[d];
    }
    struct gkyl_range blk_rng; // block range
    gkyl_range_init(&blk_rng, grid->ndim, loidx_i, upidx_i);

    struct gkyl_range inter; // intersection
    int not_empty = gkyl_range_intersect(&inter, &blk_rng, range);
    if (not_empty) {

      buff = gkyl_mem_buff_resize(buff, sz*hdr->esznc);
      if (1 != fread(gkyl_mem_buff_data(buff), sz*hdr->esznc, 1, fp)) {
        gkyl_mem_buff_release(buff);
        return GKYL_ARRAY_RIO_FREAD_FAILED;
      }

      struct gkyl_range_iter iter;
      gkyl_range_iter_init(&iter, &inter);
      while (gkyl_range_iter_next(&iter)) {
        char *out = gkyl_array_fetch(arr, gkyl_range_idx(range, iter.idx));
        const char *inp = gkyl_mem_buff_data(buff) + hdr->esznc*gkyl_range_idx(&blk_rng, iter.idx);
        memcpy(out, inp, hdr->esznc);
      }
    }

    loc += rng_sz + sz*hdr->esznc;
  }

  gkyl_mem_buff_release(buff);
  
  return GKYL_ARRAY_RIO_SUCCESS;
}

enum gkyl_array_rio_status
gkyl_grid_sub_array_read(struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  struct gkyl_array *arr, const char *fname)
{
  enum gkyl_array_rio_status status = GKYL_ARRAY_RIO_FOPEN_FAILED;
  struct gkyl_array_header_info hdr;
  FILE *fp = 0;
  with_file (fp, fname, "r") {
    grid_sub_array_header_read_fp(grid, &hdr, false, fp);
    
    if (hdr.file_type == 1)
      status = grid_sub_array_read_ft_1(grid, &hdr, range, arr, fp);
    if (hdr.file_type == 3)
      status = grid_sub_array_read_ft_3(grid, &hdr, range, arr, fp);
  }
  return status;
}

struct gkyl_array*
gkyl_grid_array_new_from_file(struct gkyl_rect_grid *grid, const char* fname)
{
  struct gkyl_array *arr = 0;
  struct gkyl_array_header_info hdr;

  enum gkyl_array_rio_status status = GKYL_ARRAY_RIO_FREAD_FAILED;
  FILE *fp = 0;
  with_file (fp, fname, "r") {
    status = grid_sub_array_header_read_fp(grid, &hdr, false, fp);
  }

  if (status != GKYL_ARRAY_RIO_SUCCESS)
    return 0;

  size_t nc = hdr.esznc/gkyl_elem_type_size[hdr.etype];
  arr = gkyl_array_new(hdr.etype, nc, hdr.tot_cells);
  struct gkyl_range range;
  gkyl_range_init_from_shape1(&range, grid->ndim, grid->cells);

  status = gkyl_grid_sub_array_read(grid, &range, arr, fname);

  if (status != GKYL_ARRAY_RIO_SUCCESS) {
    gkyl_array_release(arr);
    arr = 0;
  }
    
  return arr;
}
