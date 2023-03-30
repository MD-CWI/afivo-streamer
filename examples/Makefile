AF_DIR  := ..
LIBDIRS := $(AF_DIR)/external_libraries/silo/lib $(AF_DIR)/external_libraries/hypre/lib
LIBS	:= afivo silo HYPRE
OUTDIR	:= output

include $(AF_DIR)/src/makerules.make

.PHONY: all clean run run_2d run_3d

PROGS_XD := random_refinement poisson_basic poisson_benchmark advection		\
 computational_domain boundary_conditions poisson_neumann particles_to_grid	\
 ghostcell_benchmark particles_gravity poisson_coarse_solver implicit_diffusion	\
 poisson_dielectric poisson_helmholtz reaction_diffusion advection_v2		\
 electrode_example poisson_lsf_test stencil_test poisson_stencil		\
 dielectric_surface time_integration helmholtz_variable_stencil			\
 check_ghostcells check_prolongation

PROGS_1D := $(PROGS_XD:%=%_1d)

PROGS_2D := $(PROGS_XD:%=%_2d) poisson_cyl poisson_cyl_dielectric		\
simple_streamer poisson_cyl_analytic poisson_helmholtz_cyl solid_body_rotation	\
amr_solid_body_rotation euler_gas_dynamics dielectric_surface			\
electrode_dielectric electrode_example two_electrodes_2d incompressible_flow_2d	\
compressible_flow_wall_2d

PROGS_3D := $(PROGS_XD:%=%_3d) poisson_div_cleaning two_electrodes_3d

PROGS := $(PROGS_1D) $(PROGS_2D) $(PROGS_3D)

all:	$(PROGS)

clean:
	$(RM) $(PROGS) *.o *.mod

run:	$(PROGS)
	for prog in $(PROGS); do\
		echo "Running $$prog";\
		./$$prog || { echo "$$prog failed" ; exit 1; };\
	done

run_2d:
	for prog in $(PROGS_2D); do\
		echo "Running $$prog";\
		./$$prog || { echo "$$prog failed" ; exit 1; };\
	done

run_3d:
	for prog in $(PROGS_3D); do\
		echo "Running $$prog";\
		./$$prog || { echo "$$prog failed" ; exit 1; };\
	done

$(OUTDIR):
	mkdir -p $@

# Set flags
$(PROGS_1D): $(AF_DIR)/lib_1d/libafivo.a | $(OUTDIR)
$(PROGS_1D): INCDIRS+=$(AF_DIR)/lib_1d
$(PROGS_1D): LIBDIRS+=$(AF_DIR)/lib_1d
$(PROGS_1D): FFLAGS += -DNDIM=1

$(PROGS_2D): $(AF_DIR)/lib_2d/libafivo.a | $(OUTDIR)
$(PROGS_2D): INCDIRS+=$(AF_DIR)/lib_2d
$(PROGS_2D): LIBDIRS+=$(AF_DIR)/lib_2d
$(PROGS_2D): FFLAGS += -DNDIM=2

$(PROGS_3D): $(AF_DIR)/lib_3d/libafivo.a | $(OUTDIR)
$(PROGS_3D): INCDIRS+=$(AF_DIR)/lib_3d
$(PROGS_3D): LIBDIRS+=$(AF_DIR)/lib_3d
$(PROGS_3D): FFLAGS += -DNDIM=3

# How to create executables
%_1d: %.f90
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -I,$(INCDIRS)) \
	$(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
%_2d: %.f90
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -I,$(INCDIRS)) \
	$(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
%_3d: %.f90
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -I,$(INCDIRS)) \
	$(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

# Some Poisson programs and objects depend on m_gaussians and m_config
PROGS_POIS	:= $(filter poisson_%, $(PROGS))
$(PROGS_POIS): 	m_gaussians.o m_config.o
