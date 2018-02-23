SRC_ND := m_afivo_types.f90 m_morton.f90 m_vtk.f90 m_write_silo.f90	\
	m_flux_schemes.f90
SRC_XD := m_aX_types.f90 m_aX_core.f90 m_aX_output.f90 m_aX_ghostcell.f90	\
	m_aX_restrict.f90 m_aX_prolong.f90 m_aX_utils.f90 m_aX_multigrid.f90	\
	m_aX_interp.f90 m_aX_particles.f90 m_aX_all.f90

SRC_2D	:= $(SRC_XD:m_aX_%.f90=m_a2_%.f90)
SRC_3D	:= $(SRC_XD:m_aX_%.f90=m_a3_%.f90)
SRCS	:= $(SRC_ND) $(SRC_2D) $(SRC_3D)

OBJS_ND := $(SRC_ND:%.f90=%.o)
OBJS_2D := $(SRC_2D:%.f90=%.o)
OBJS_3D := $(SRC_3D:%.f90=%.o)
OBJS	:= $(OBJS_ND) $(OBJS_2D) $(OBJS_3D)

INCDIRS	:=
LIBDIRS :=
LIBS	:=

include ../makerules.make

.PHONY: srcs clean allclean

libafivo.a: $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

srcs:	$(SRCS)

clean:
	$(RM) libafivo.a *.o *.mod

allclean: clean
	$(RM) $(SRC_2D) $(SRC_3D)

m_a2_%.f90: m_aX_%.f90
	@echo "*** $< has changed, regenerating $@ ***"
	@$(PPROCESS_2D)

m_a3_%.f90: m_aX_%.f90
	@echo "*** $< has changed, regenerating $@ ***"
	@$(PPROCESS_3D)

# Requires Silo for compilation
m_write_silo.o:		INCDIRS += ../external_libraries/silo/include
m_write_silo.o:		FFLAGS += -Wno-implicit-interface

# Dependency information (generated with list_dependencies.sh)
m_a2_all.o: m_a2_core.mod
m_a2_all.o: m_a2_ghostcell.mod
m_a2_all.o: m_a2_interp.mod
m_a2_all.o: m_a2_multigrid.mod
m_a2_all.o: m_a2_output.mod
m_a2_all.o: m_a2_particles.mod
m_a2_all.o: m_a2_prolong.mod
m_a2_all.o: m_a2_restrict.mod
m_a2_all.o: m_a2_types.mod
m_a2_all.o: m_a2_utils.mod
m_a2_core.o: m_a2_ghostcell.mod
m_a2_core.o: m_a2_prolong.mod
m_a2_core.o: m_a2_restrict.mod
m_a2_core.o: m_a2_types.mod
m_a2_core.o: m_a2_utils.mod
m_a2_core.o: m_morton.mod
m_a2_ghostcell.o: m_a2_prolong.mod
m_a2_ghostcell.o: m_a2_types.mod
m_a2_interp.o: m_a2_types.mod
m_a2_interp.o: m_a2_utils.mod
m_a2_multigrid.o: m_a2_ghostcell.mod
m_a2_multigrid.o: m_a2_prolong.mod
m_a2_multigrid.o: m_a2_restrict.mod
m_a2_multigrid.o: m_a2_types.mod
m_a2_multigrid.o: m_a2_utils.mod
m_a2_output.o: m_a2_core.mod
m_a2_output.o: m_a2_interp.mod
m_a2_output.o: m_a2_types.mod
m_a2_output.o: m_vtk.mod
m_a2_output.o: m_write_silo.mod
m_a2_particles.o: m_a2_ghostcell.mod
m_a2_particles.o: m_a2_restrict.mod
m_a2_particles.o: m_a2_types.mod
m_a2_particles.o: m_a2_utils.mod
m_a2_prolong.o: m_a2_types.mod
m_a2_restrict.o: m_a2_types.mod
m_a2_types.o: m_afivo_types.mod
m_a2_utils.o: m_a2_types.mod
m_a3_all.o: m_a3_core.mod
m_a3_all.o: m_a3_ghostcell.mod
m_a3_all.o: m_a3_interp.mod
m_a3_all.o: m_a3_multigrid.mod
m_a3_all.o: m_a3_output.mod
m_a3_all.o: m_a3_particles.mod
m_a3_all.o: m_a3_prolong.mod
m_a3_all.o: m_a3_restrict.mod
m_a3_all.o: m_a3_types.mod
m_a3_all.o: m_a3_utils.mod
m_a3_core.o: m_a3_ghostcell.mod
m_a3_core.o: m_a3_prolong.mod
m_a3_core.o: m_a3_restrict.mod
m_a3_core.o: m_a3_types.mod
m_a3_core.o: m_a3_utils.mod
m_a3_core.o: m_morton.mod
m_a3_ghostcell.o: m_a3_prolong.mod
m_a3_ghostcell.o: m_a3_types.mod
m_a3_interp.o: m_a3_types.mod
m_a3_interp.o: m_a3_utils.mod
m_a3_multigrid.o: m_a3_ghostcell.mod
m_a3_multigrid.o: m_a3_prolong.mod
m_a3_multigrid.o: m_a3_restrict.mod
m_a3_multigrid.o: m_a3_types.mod
m_a3_multigrid.o: m_a3_utils.mod
m_a3_output.o: m_a3_core.mod
m_a3_output.o: m_a3_interp.mod
m_a3_output.o: m_a3_types.mod
m_a3_output.o: m_vtk.mod
m_a3_output.o: m_write_silo.mod
m_a3_particles.o: m_a3_ghostcell.mod
m_a3_particles.o: m_a3_restrict.mod
m_a3_particles.o: m_a3_types.mod
m_a3_particles.o: m_a3_utils.mod
m_a3_prolong.o: m_a3_types.mod
m_a3_restrict.o: m_a3_types.mod
m_a3_types.o: m_afivo_types.mod
m_a3_utils.o: m_a3_types.mod
