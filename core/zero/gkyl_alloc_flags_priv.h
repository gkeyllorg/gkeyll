// private header to provide some commonly used flags in allocations
#pragma once

#include <stdint.h>

// flags and corresponding bit-masks
enum gkyl_alloc_flags { GKYL_IS_CU_ALLOC, GKYL_IS_ALLOC_ALIGNED, GKYL_IS_ALLOC_EXTERN };
static const uint32_t gkyl_alloc_flags_masks[] =
{ 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };

// NV-GPU flags
#define GKYL_SET_CU_ALLOC(flags) ((flags) |= gkyl_alloc_flags_masks[GKYL_IS_CU_ALLOC])
#define GKYL_CLEAR_CU_ALLOC(flags) ((flags) &= ~gkyl_alloc_flags_masks[GKYL_IS_CU_ALLOC])
#define GKYL_IS_CU_ALLOC(flags) (((flags) & gkyl_alloc_flags_masks[GKYL_IS_CU_ALLOC]) != 0)

// Alignment flags
#define GKYL_SET_ALLOC_ALIGNED(flags) ((flags) |= gkyl_alloc_flags_masks[GKYL_IS_ALLOC_ALIGNED])
#define GKYL_CLEAR_ALLOC_ALIGNED(flags) ((flags) &= ~gkyl_alloc_flags_masks[GKYL_IS_ALLOC_ALIGNED])
#define GKYL_IS_ALLOC_ALIGNED(flags) (((flags) & gkyl_alloc_flags_masks[GKYL_IS_ALLOC_ALIGNED]) != 0)

// Flag to indicate if there was no actual allocation but only external allocated memory was used
#define GKYL_SET_ALLOC_EXTERN(flags) ((flags) |= gkyl_alloc_flags_masks[GKYL_IS_ALLOC_EXTERN])
#define GKYL_CLEAR_ALLOC_EXTERN(flags) ((flags) &= ~gkyl_alloc_flags_masks[GKYL_IS_ALLOC_EXTERN])
#define GKYL_IS_ALLOC_EXTERN(flags) (((flags) & gkyl_alloc_flags_masks[GKYL_IS_ALLOC_EXTERN]) != 0)
