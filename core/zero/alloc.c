#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_util.h>

// Output command for use in debugging memory usage
#define GKYL_MEMMSG(fmt, ...) do {              \
      if (gkyl_mem_debug)                       \
        fprintf(stderr, fmt, __VA_ARGS__);      \
  } while (0);

// Output command for use in debugging CUDA memory usage
#define GKYL_CU_MEMMSG(fmt, ...) do {           \
      if (gkyl_cu_dev_mem_debug)                \
        fprintf(stderr, fmt, __VA_ARGS__);      \
    } while (0);

// by default, do not print memory allocation traces
static bool gkyl_mem_debug = false;
static bool gkyl_cu_dev_mem_debug = false;

void gkyl_mem_debug_set(bool flag)
{
  gkyl_mem_debug = flag;
}

void gkyl_cu_dev_mem_debug_set(bool flag)
{
  gkyl_cu_dev_mem_debug = flag;
}

// Compute first 'align' boundary after 'num'
#define align_up(num, align)                    \
    (((num) + ((align) - 1)) & ~((align) - 1))

static const size_t PTR_OFFSET_SZ = sizeof(uint16_t);

void*
gkyl_malloc_(const char *file, int line, const char *func, size_t size)
{
  void *mem = malloc(size);
  GKYL_MEMMSG("%p [%zu] 0.malloc: %s %s:%d\n", mem, size, file, func, line);  
  if (0 == mem) gkyl_exit("malloc failed!");
  return mem;
}

void*
gkyl_calloc_(const char *file, int line, const char *func, size_t num, size_t size)
{
  void *mem = calloc(num, size);
  GKYL_MEMMSG("%p [%zu] 0.calloc: %s %s:%d\n", mem, size, file, func, line);
  if (0 == mem) gkyl_exit("calloc failed!");
  return mem;
}

void*
gkyl_realloc_(const char *file, int line, const char *func, void *ptr, size_t new_size)
{
  void *mem = realloc(ptr, new_size);
  GKYL_MEMMSG("%p [%zu] 0.realloc: %s %s:%d\n", mem, new_size, file, func, line);  
  if (0 == mem) gkyl_exit("realloc failed!");
  return mem;
}

void
gkyl_free_(const char *file, int line, const char *func, void *ptr)
{
  GKYL_MEMMSG("%p 1.free: %s %s:%d\n", ptr, file, func, line);
  free(ptr);
}

void*
gkyl_aligned_alloc_(const char *file, int line, const char *func,
  size_t align, size_t size)
{
  void *ptr = 0;
  assert((align & (align-1)) == 0); // power of 2?

  if (align && size) {
    uint32_t hdr_size = PTR_OFFSET_SZ + (align-1);
    void *p = gkyl_calloc(size+hdr_size, 1);
    if (p) {
      ptr = (void *) align_up(((uintptr_t)p + PTR_OFFSET_SZ), align);
      *((uint16_t *)ptr - 1) = (uint16_t) ((uintptr_t) ptr - (uintptr_t) p);
    }
  }

  GKYL_MEMMSG("%p 0.aligned_alloc: %s %s:%d\n", ptr, file, func, line);
  return ptr;
}

void*
gkyl_aligned_realloc_(const char *file, int line, const char *func,
  void *ptr, size_t align, size_t old_sz, size_t new_sz)
{
  void *nptr = gkyl_aligned_alloc(align, new_sz);
  if (0 == nptr) {
    gkyl_exit("aligned_realloc failed!");
    return 0;
  }
  memcpy(nptr, ptr, old_sz < new_sz ? old_sz : new_sz);
  gkyl_aligned_free(ptr);

  GKYL_MEMMSG("%p 0.aligned_realloc: %s %s:%d\n", ptr, file, func, line);
  return nptr;
}

void
gkyl_aligned_free_(const char *file, int line, const char *func,
  void* ptr)
{
  assert(ptr);
  GKYL_MEMMSG("%p 1.aligned_free: %s %s:%d\n", ptr, file, func, line);
    
  uint16_t offset = *((uint16_t *)ptr - 1);
  gkyl_free((uint8_t *)ptr - offset);
}

struct gkyl_mem_buff_tag {
  bool on_gpu; // is this on GPU?
  size_t count; // size of memory in bytes
  char *data; // Allocated memory
};

gkyl_mem_buff
gkyl_mem_buff_new(size_t count)
{
  struct gkyl_mem_buff_tag *mem = gkyl_malloc(sizeof(*mem));
  mem->on_gpu = false;
  mem->count = count;
  mem->data = gkyl_malloc(count);
  return mem;
}

gkyl_mem_buff
gkyl_mem_buff_cu_new(size_t count)
{
  struct gkyl_mem_buff_tag *mem = gkyl_malloc(sizeof(*mem));
  mem->on_gpu = true;
  mem->count = count;
  mem->data = gkyl_cu_malloc(count);
  return mem;
}

gkyl_mem_buff
gkyl_mem_buff_resize(gkyl_mem_buff mem, size_t count)
{
  if (count > mem->count) {
    if (mem->on_gpu) {
      char *data_new = gkyl_cu_malloc(count);
      gkyl_cu_memcpy(data_new, mem->data, mem->count, GKYL_CU_MEMCPY_D2D);
      gkyl_cu_free(mem->data);
      mem->data = data_new;
    }
    else {
      mem->data = gkyl_realloc(mem->data, count);
    }
    mem->count = count;
  }
  return mem;
}

size_t
gkyl_mem_buff_size(gkyl_mem_buff mem)
{
  return mem->count;
}

char*
gkyl_mem_buff_data(gkyl_mem_buff mem)
{
  return mem->data;
}

void
gkyl_mem_buff_release(gkyl_mem_buff mem)
{
  if (mem->on_gpu)
    gkyl_cu_free(mem->data);
  else
    gkyl_free(mem->data);
  
  gkyl_free(mem);
}

// CUDA specific code

#ifdef GKYL_HAVE_CUDA

#include <cuda_runtime.h>

void*
gkyl_cu_malloc_(const char *file, int line, const char *func, size_t size)
{
  void *ptr;
  cudaError_t err = cudaMalloc(&ptr, size);
  if (err != cudaSuccess)
    gkyl_exit("cudaMalloc failed!");
  
  GKYL_CU_MEMMSG("%p 0.cudaMalloc: %s %s:%d\n", ptr, file, func, line);

  return ptr;
}

void*
gkyl_cu_malloc_host_(const char *file, int line, const char *func, size_t size)
{
  // Allocate pinned host memory.
  void *ptr;
  cudaError_t err = cudaMallocHost(&ptr, size);
  if (err != cudaSuccess)
    gkyl_exit("cudaMallocHost failed!");

  GKYL_CU_MEMMSG("%p 0.cudaMallocHost: %s %s:%d\n", ptr, file, func, line);

  return ptr;
}

void
gkyl_cu_free_(const char *file, int line, const char *func, void *ptr)
{
  GKYL_CU_MEMMSG("%p 1.cudaFree: %s %s:%d\n", ptr, file, func, line);
  cudaFree(ptr);
}

void
gkyl_cu_free_host_(const char *file, int line, const char *func, void *ptr)
{
  GKYL_CU_MEMMSG("%p 1.cudaFreeHost: %s %s:%d\n", ptr, file, func, line);
  cudaFreeHost(ptr);
}

void
gkyl_cu_memcpy(void *dst, const void *src, size_t count, enum gkyl_cu_memcpy_kind kind)
{
  cudaError_t err = cudaMemcpy(dst, src, count, kind);
  if (err != cudaSuccess) {
    char str[1024];
    sprintf(str, "\nCUDA error: %s\n", cudaGetErrorString(err));
    gkyl_exit(str);
  }
}

void
gkyl_cu_memcpy_async(void *dst, const void *src, size_t count, enum gkyl_cu_memcpy_kind kind, cudaStream_t stream)
{
  cudaError_t err = cudaMemcpyAsync(dst, src, count, kind, stream);
  if (err != cudaSuccess) {
    char str[1024];
    sprintf(str, "\nCUDA error: %s\n", cudaGetErrorString(err));
    gkyl_exit(str);
  }
}

void
gkyl_cu_memset(void *data, int val, size_t count)
{
  cudaError_t err = cudaMemset(data, val, count);
  if (err != cudaSuccess)
    gkyl_exit("gkyl_cu_memset failed!");
}

#else

// These non-CUDA functions will simply abort. When not using CUDA
// none of these methods should be called at all.

void*
gkyl_cu_malloc_(const char *file, int line, const char *func, size_t size)
{
  assert(false);
  return 0;
}

void*
gkyl_cu_malloc_host_(const char *file, int line, const char *func, size_t size)
{
  assert(false);
  return 0;
}

void
gkyl_cu_free_(const char *file, int line, const char *func, void *ptr)
{
  assert(false);
}

void
gkyl_cu_free_host_(const char *file, int line, const char *func, void *ptr)
{
  assert(false);
}

void
gkyl_cu_memcpy(void *dst, const void *src, size_t count, enum gkyl_cu_memcpy_kind kind)
{
  assert(false);  
}

void
gkyl_cu_memcpy_async(void *dst, const void *src, size_t count, enum gkyl_cu_memcpy_kind kind, int stream)
{
  assert(false);  
}

void
gkyl_cu_memset(void *data, int val, size_t count)
{
  assert(false);  
}

#endif // CUDA specific code
