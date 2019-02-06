AF_DIR  := ..
LIBDIRS := $(AF_DIR)/external_libraries/silo/lib $(AF_DIR)/external_libraries/hypre/lib
LIBS	:= afivo silo HYPRE
OUTDIR	:= output

include $(AF_DIR)/src/makerules.make

.PHONY: all clean run

PROGS_XD := random_refinement poisson_basic poisson_benchmark advection		\
implicit_diffusion computational_domain boundary_conditions poisson_neumann	\
particles_to_grid poisson_dielectric ghostcell_benchmark poisson_helmholtz	\
particles_gravity poisson_coarse_solver

PROGS_2D := $(PROGS_XD:%=%_2d) poisson_cyl poisson_cyl_dielectric	\
simple_streamer poisson_cyl_analytic poisson_helmholtz_cyl

PROGS_3D := $(PROGS_XD:%=%_3d)

PROGS := $(PROGS_2D) $(PROGS_3D)

all:	$(PROGS)

clean:
	$(RM) $(PROGS) *.o *.mod

run:
	for prog in $(PROGS); do\
		./$$prog || { echo "$$prog failed" ; exit 1; };\
	done


$(OUTDIR):
	mkdir -p $@

# Set flags for 2d and 3d version of programs
$(PROGS_2D): $(AF_DIR)/lib_2d/libafivo.a | $(OUTDIR)
$(PROGS_2D): INCDIRS+=$(AF_DIR)/lib_2d
$(PROGS_2D): LIBDIRS+=$(AF_DIR)/lib_2d
$(PROGS_2D): FFLAGS += -DNDIM=2

$(PROGS_3D): $(AF_DIR)/lib_3d/libafivo.a | $(OUTDIR)
$(PROGS_3D): INCDIRS+=$(AF_DIR)/lib_3d
$(PROGS_3D): LIBDIRS+=$(AF_DIR)/lib_3d
$(PROGS_3D): FFLAGS += -DNDIM=3

# How to create executables
%_2d: %.f90
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -I,$(INCDIRS)) \
	$(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
%_3d: %.f90
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -I,$(INCDIRS)) \
	$(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

# Poisson programs and objects depend on m_gaussians
PROGS_POIS	:= $(filter poisson_%, $(PROGS))
OBJS_POIS	:= $(PROGS_POIS:%=%.o)
$(PROGS_POIS): 	m_gaussians.o
$(OBJS_POIS):	m_gaussians.o
