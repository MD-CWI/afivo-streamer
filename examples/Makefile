A5_DIR  := ..
INCDIRS	:= $(A5_DIR)/src
LIBDIRS := $(A5_DIR)/src $(A5_DIR)/external_libraries/silo/lib
LIBS	:= afivo silo
OUTDIR	:= output

include $(A5_DIR)/makerules.make

# All program sources (the .f90 files except modules)
PROGS	:= random_refinement_2d random_refinement_3d poisson_basic_2d poisson_basic_3d poisson_cyl poisson_cyl_dielectric poisson_benchmark_2d poisson_benchmark_3d drift_diffusion_2d drift_diffusion_3d implicit_diffusion_2d

.PHONY: all clean

all:	$(PROGS)

clean:
	$(RM) $(PROGS) *.o *.mod

$(OUTDIR):
		mkdir -p $@

%_2d.f90: %_Xd.f90
	@echo "*** $< has changed, regenerating $@ ***"
	@$(PPROCESS_2D)

%_3d.f90: %_Xd.f90
	@echo "*** $< has changed, regenerating $@ ***"
	@$(PPROCESS_3D)

# Dependency information
$(PROGS): $(A5_DIR)/src/libafivo.a | $(OUTDIR)

# Poisson programs and objects depend on m_gaussians
PROGS_POIS	:= $(filter poisson_%, $(PROGS))
OBJS_POIS	:= $(PROGS_POIS:%=%.o)
$(PROGS_POIS): 	m_gaussians.o
$(OBJS_POIS):	m_gaussians.o
