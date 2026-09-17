#pragma once

#include <gkyl_util.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

/**
 * Set the global flag to turn on memory allocation/deallocation
 * tracing.
 *
 * @param flag Flag to set
 */
void gkyl_mem_debug_set(bool flag);

/**
 * Set the global flag to turn on cuda memory allocation/deallocation
 * tracing.
 *
 * @param flag Flag to set
 */
void gkyl_cu_dev_mem_debug_set(bool flag);

// The following allocators are implemented as macros to allow
// debugging memory leaks. In general, it is best to use
// valgrind. However, valgrind can be very slow and also it gets
// confused with CUDA symbols in the executable (even for host
// code). Further, it seems that cuda-memcheck is not really tracking
// memory leaks.

#define gkyl_malloc(size)                                       \
    gkyl_malloc_(__FILE__, __LINE__, __FUNCTION__, size)
#define gkyl_calloc(num, size)                                  \
    gkyl_calloc_(__FILE__, __LINE__, __FUNCTION__, num, size)
#define gkyl_realloc(ptr, new_size)                                     \
    gkyl_realloc_(__FILE__, __LINE__, __FUNCTION__, ptr, new_size)
#define gkyl_free(ptr) \
    gkyl_free_(__FILE__, __LINE__, __FUNCTION__, ptr)

// Allocate memory that is aligned to given byte boundary. You must
// only use gkyl_aligned_realloc() and gkyl_aligned_free() methods to
// reallocate or free memory returned by this methods.

#define gkyl_aligned_alloc(align, size)                                 \
    gkyl_aligned_alloc_(__FILE__, __LINE__, __FUNCTION__, align, size)
#define gkyl_aligned_realloc(ptr, align, old_sz, new_sz)                       \
    gkyl_aligned_realloc_(__FILE__, __LINE__, __FUNCTION__, ptr, align, old_sz, new_sz)
#define gkyl_aligned_free(ptr) \
    gkyl_aligned_free_(__FILE__, __LINE__, __FUNCTION__, ptr)

// The following allocators have the same calling/return behavior as
// standard C allocators. However, an error is signaled if allocation
// fails.

void* gkyl_malloc_(const char *file, int line, const char *func, size_t size);
void* gkyl_calloc_(const char *file, int line, const char *func, size_t num, size_t size);
void* gkyl_realloc_(const char *file, int line, const char *func, void *ptr, size_t new_size);
void gkyl_free_(const char *file, int line, const char *func, void *ptr);

/**
 * Allocate memory that is aligned to given byte boundary. You must
 * only use gkyl_aligned_realloc() and gkyl_aligned_free() methods to
 * reallocate or free memory returned by this method.
 *
 * @param align Alignment boundary. Must be power of 2.
 * @param size Number of bytes to allocate.
 */
void* gkyl_aligned_alloc_(const char *file, int line, const char *func, size_t align, size_t size);

/**
 * Reallocate memory that is aligned to given byte boundary. The
 * pointer passed to this must be allocated by
 * gkyl_aligned_alloc(). Returned pointer must be freed using
 * gkyl_aligned_free() method.
 *
 * @param ptr Pointer to reallocate
 * @param align Alignment boundary. Must be power of 2.
 * @param old_sz Old size of memory.
 * @param new_sz New size of memory.
 */
void* gkyl_aligned_realloc_(const char *file, int line, const char *func,
  void *ptr, size_t align, size_t old_sz, size_t new_sz);

/**
 * Free memory allocated by gkyl_aligned_alloc().
 *
 * @param ptr Memory to free.
 */
void gkyl_aligned_free_(const char *file, int line, const char *func, void *ptr);

// Represents a sized chunk of memory
typedef struct gkyl_mem_buff_tag* gkyl_mem_buff;

/** Allocate new memory buffer with count bytes */
gkyl_mem_buff gkyl_mem_buff_new(size_t count);

/** Allocate new memory buffer on GPU with count bytes */
gkyl_mem_buff gkyl_mem_buff_cu_new(size_t count);

/** Resize the mem_buff */
gkyl_mem_buff gkyl_mem_buff_resize(gkyl_mem_buff mem, size_t count);

/** Get size of mem_buff */
size_t gkyl_mem_buff_size(gkyl_mem_buff mem);

/** Get pointer to data buffer */
char* gkyl_mem_buff_data(gkyl_mem_buff mem);

/** Free buffer */
void gkyl_mem_buff_release(gkyl_mem_buff mem);

// CUDA specific code (NV: Nvidia)

#define gkyl_cu_malloc(size)                                    \
    gkyl_cu_malloc_(__FILE__, __LINE__, __FUNCTION__, size)
#define gkyl_cu_free(ptr)                                       \
    gkyl_cu_free_(__FILE__, __LINE__, __FUNCTION__, ptr)

#define gkyl_cu_malloc_host(size)                               \
    gkyl_cu_malloc_host_(__FILE__, __LINE__, __FUNCTION__, size)
#define gkyl_cu_free_host(ptr)                                  \
    gkyl_cu_free_host_(__FILE__, __LINE__, __FUNCTION__, ptr)

/** Allocate memory on NV-GPU */
void* gkyl_cu_malloc_(const char *file, int line, const char *func, size_t size);

/** Free memory on device */
void gkyl_cu_free_(const char *file, int line, const char *func, void *ptr);

/** Allocate pinned host memory on NV-GPU */
void* gkyl_cu_malloc_host_(const char *file, int line, const char *func, size_t size);

/** Free pinned host memory on device */
void gkyl_cu_free_host_(const char *file, int line, const char *func, void *ptr);

/** Copy data between host/device */
void gkyl_cu_memcpy(void *dst, const void *src, size_t count, enum gkyl_cu_memcpy_kind kind);

/** Copy data between host/device */
#ifdef GKYL_HAVE_CUDA
void gkyl_cu_memcpy_async(void *dst, const void *src, size_t count, enum gkyl_cu_memcpy_kind kind, cudaStream_t stream);
#else
void gkyl_cu_memcpy_async(void *dst, const void *src, size_t count, enum gkyl_cu_memcpy_kind kind, int stream);
#endif

/** Set memory on device */
void gkyl_cu_memset(void *data, int val, size_t count);
