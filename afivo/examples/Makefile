A5_DIR  := ..
INCDIRS	:= $(A5_DIR)/src
LIBDIRS := $(A5_DIR)/src $(A5_DIR)/external_libraries/silo/lib
LIBS	:= afivo silo
OUTDIR	:= output

include $(A5_DIR)/makerules.make

# All program sources (the .f90 files except modules)
PROGS := random_refinement_2d random_refinement_3d poisson_basic_2d		\
poisson_basic_3d poisson_cyl poisson_cyl_dielectric poisson_benchmark_2d	\
poisson_benchmark_3d advection_2d advection_3d implicit_diffusion_2d		\
implicit_diffusion_3d simple_streamer_2d computational_domain_2d		\
computational_domain_3d boundary_conditions_2d boundary_conditions_3d		\
poisson_neumann_2d poisson_neumann_3d

.PHONY: all clean

# Keep 2d and 3d source files
.PRECIOUS: %_2d.f90 %_3d.f90

all:	$(PROGS)

clean:
	$(RM) $(PROGS) *.o *.mod

$(OUTDIR):
		mkdir -p $@

%_2d.f90: Xd_files/%_Xd.f90
	@echo "*** $< has changed, regenerating $@ ***"
	@$(PPROCESS_2D)

%_3d.f90: Xd_files/%_Xd.f90
	@echo "*** $< has changed, regenerating $@ ***"
	@$(PPROCESS_3D)

# Dependency information
$(PROGS): $(A5_DIR)/src/libafivo.a | $(OUTDIR)

# Poisson programs and objects depend on m_gaussians
PROGS_POIS	:= $(filter poisson_%, $(PROGS))
OBJS_POIS	:= $(PROGS_POIS:%=%.o)
$(PROGS_POIS): 	m_gaussians.o
$(OBJS_POIS):	m_gaussians.o
