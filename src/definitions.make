SRC := m_morton.f90 m_vtk.f90 m_write_silo.f90 m_af_flux_schemes.f90		\
	m_af_types.f90 m_af_core.f90 m_af_output.f90 m_af_ghostcell.f90		\
	m_af_restrict.f90 m_af_prolong.f90 m_af_utils.f90 m_af_multigrid.f90	\
	m_mg_types.f90 m_af_interp.f90 m_af_particles.f90 m_af_all.f90	\
	m_coarse_solver.f90

OBJS := $(SRC:%.f90=%.o)

# Requires Silo for compilation
m_write_silo.o:		INCDIRS += ../external_libraries/silo/include
m_write_silo.o:		FFLAGS += -Wno-implicit-interface
m_coarse_solver.o:	FFLAGS += -Wno-implicit-interface

# Dependency information (generated with list_dependencies.sh)
m_af_all.o: m_af_core.mod
m_af_all.o: m_af_flux_schemes.mod
m_af_all.o: m_af_ghostcell.mod
m_af_all.o: m_af_interp.mod
m_af_all.o: m_af_multigrid.mod
m_af_all.o: m_af_output.mod
m_af_all.o: m_af_particles.mod
m_af_all.o: m_af_prolong.mod
m_af_all.o: m_af_restrict.mod
m_af_all.o: m_af_types.mod
m_af_all.o: m_af_utils.mod
m_af_core.o: m_af_ghostcell.mod
m_af_core.o: m_af_prolong.mod
m_af_core.o: m_af_restrict.mod
m_af_core.o: m_af_types.mod
m_af_core.o: m_af_utils.mod
m_af_core.o: m_morton.mod
m_af_ghostcell.o: m_af_prolong.mod
m_af_ghostcell.o: m_af_types.mod
m_af_interp.o: m_af_types.mod
m_af_interp.o: m_af_utils.mod
m_af_multigrid.o: m_af_ghostcell.mod
m_af_multigrid.o: m_af_prolong.mod
m_af_multigrid.o: m_af_restrict.mod
m_af_multigrid.o: m_af_types.mod
m_af_multigrid.o: m_af_utils.mod
m_af_multigrid.o: m_coarse_solver.mod
m_af_multigrid.o: m_mg_types.mod
m_af_output.o: m_af_core.mod
m_af_output.o: m_af_interp.mod
m_af_output.o: m_af_types.mod
m_af_output.o: m_vtk.mod
m_af_output.o: m_write_silo.mod
m_af_particles.o: m_af_ghostcell.mod
m_af_particles.o: m_af_restrict.mod
m_af_particles.o: m_af_types.mod
m_af_particles.o: m_af_utils.mod
m_af_prolong.o: m_af_types.mod
m_af_restrict.o: m_af_types.mod
m_af_utils.o: m_af_types.mod
m_coarse_solver.o: m_af_types.mod
m_coarse_solver.o: m_mg_types.mod
m_mg_types.o: m_af_types.mod
