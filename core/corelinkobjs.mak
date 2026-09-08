# -*- makefile-gmake -*-

# Core include paths
CORE_INCS = -I../core/$(KERNELS_DIR)/array_average -I../core/$(KERNELS_DIR)/array_integrate -I../core/$(KERNELS_DIR)/basis -I../core/$(KERNELS_DIR)/bin_op -I../core/$(KERNELS_DIR)/dg_interpolate -I../core/$(KERNELS_DIR)/skin_surf_from_ghost -I../core/apps -I../core/zero -I../core/minus -I../core/minus/STC/include -I../zero

# Core link objects

CORE_LINK_OBJS = ../$(BUILD_DIR)/core/minus/*.o ../$(BUILD_DIR)/core/minus/lapack_lite/*.o ../$(BUILD_DIR)/core/$(KERNELS_DIR)/*/*.o ../$(BUILD_DIR)/core/zero/*.o ../$(BUILD_DIR)/core/apps/*.o
# ... with GPUs
CORE_LINK_CU_OBJS = ../$(BUILD_DIR)/core/unit/*.cu.o
