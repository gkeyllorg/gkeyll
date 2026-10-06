# -*- makefile-gmake -*-

# Gyrokinetic include objects
GYROKINETIC_INCS = -I../gyrokinetic/$(KERNELS_DIR)/ambi_bolt_potential -I../gyrokinetic/$(KERNELS_DIR)/bgk_gyrokinetic -I../gyrokinetic/$(KERNELS_DIR)/deflate_geo -I../gyrokinetic/$(KERNELS_DIR)/deflate_surf -I../gyrokinetic/$(KERNELS_DIR)/derived_geo -I../gyrokinetic/$(KERNELS_DIR)/dg_diffusion_gyrokinetic -I../gyrokinetic/$(KERNELS_DIR)/fem_parproj -I../gyrokinetic/$(KERNELS_DIR)/fem_poisson_perp -I../gyrokinetic/$(KERNELS_DIR)/gyrokinetic -I../gyrokinetic/$(KERNELS_DIR)/gyrokinetic_pol_density -I../gyrokinetic/$(KERNELS_DIR)/inflate_surf -I../gyrokinetic/$(KERNELS_DIR)/lbo_gyrokinetic -I../gyrokinetic/$(KERNELS_DIR)/neutral -I../gyrokinetic/$(KERNELS_DIR)/positivity_shift_gyrokinetic -I../gyrokinetic/$(KERNELS_DIR)/rad -I../gyrokinetic/$(KERNELS_DIR)/translate_dim -I../gyrokinetic/$(KERNELS_DIR)/twistshift_dg -I../gyrokinetic/apps -I../gyrokinetic/zero

# Gyrokinetic link objects

GYROKINETIC_LINK_OBJS = ../$(BUILD_DIR)/gyrokinetic/zero/*.o ../$(BUILD_DIR)/gyrokinetic/apps/*.o ../$(BUILD_DIR)/gyrokinetic/data/adas/*.o ../$(BUILD_DIR)/gyrokinetic/$(KERNELS_DIR)/*/*.o
# ... with GPUs
GYROKINETIC_LINK_CU_OBJS = 
