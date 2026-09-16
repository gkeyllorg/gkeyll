# -*- makefile-gmake -*-

# Vlasov include paths
VLASOV_INCS = -I../vlasov/$(KERNELS_DIR)/advection -I../vlasov/$(KERNELS_DIR)/bgk_vlasov -I../vlasov/$(KERNELS_DIR)/canonical_fluid -I../vlasov/$(KERNELS_DIR)/dg_diffusion_fluid -I../vlasov/$(KERNELS_DIR)/dg_gr_maxwell -I../vlasov/$(KERNELS_DIR)/dg_diffusion_gen -I../vlasov/$(KERNELS_DIR)/dg_diffusion_vlasov -I../vlasov/$(KERNELS_DIR)/euler -I../vlasov/$(KERNELS_DIR)/fpo -I../vlasov/$(KERNELS_DIR)/lbo_vlasov -I../vlasov/$(KERNELS_DIR)/maxwell -I../vlasov/$(KERNELS_DIR)/vlasov -I../vlasov/$(KERNELS_DIR)/mom_vlasov -I../vlasov/$(KERNELS_DIR)/flux_vlasov -I../vlasov/$(KERNELS_DIR)/positivity_shift_vlasov -I../vlasov/apps -I../vlasov/zero

# Vlasov link objects

VLASOV_LINK_OBJS = ../$(BUILD_DIR)/vlasov/$(KERNELS_DIR)/*/*.o ../$(BUILD_DIR)/vlasov/zero/*.o ../$(BUILD_DIR)/vlasov/apps/*.o
# ... with GPUs
VLASOV_LINK_CU_OBJS = 
