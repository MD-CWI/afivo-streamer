FC 	:= gfortran
FFLAGS	:= -Wall -O2 -fcheck=all -g -pg -ffpe-trap=invalid,zero,overflow\
	-pedantic -finit-real=nan -std=f2008 -fopenmp -cpp
# FFLAGS	:= -Ofast -fopenmp -cpp
OBJS	:= m_morton.f90 m_afivo_constants.f90 m_afivo_2d.f90 m_vtk.f90 m_mg2d.f90 m_gmres.f90

INCDIRS	:=
LIBDIRS :=
LIBS	:=

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all: test_morton test_base test_unstr_vtk test_drift_diff test_mg \
	test_gmres test_mg_elec test_omp test_base_3d

clean:
	$(RM) *.o *.mod test_m_morton m_afivo_2d.f90 m_afivo_3d.f90

m_afivo_2d.f90: m_afivo_Xd.f90
	sed 's|$$D|2|g' < $< | cpp -traditional-cpp -P | cat -s > $@

m_afivo_3d.f90: m_afivo_Xd.f90
	sed 's|$$D|3|g' < $< | cpp -traditional-cpp -P | cat -s > $@


# Dependency information
m_afivo_2d.f90:		m_afivo_Xd.f90
m_afivo_3d.f90:		m_afivo_Xd.f90
m_afivo_2d.o:		m_vtk.o m_morton.o m_afivo_constants.o
m_afivo_3d.o:		m_vtk.o m_morton.o m_afivo_constants.o
m_mg2d.o:		m_afivo_2d.o
test_morton:		m_morton.o
test_base:		m_afivo_2d.o m_vtk.o m_morton.o
test_base_3d:		m_afivo_3d.o m_vtk.o m_morton.o
test_drift_diff:	m_afivo_2d.o m_vtk.o m_morton.o
test_mg:		m_afivo_2d.o m_vtk.o m_morton.o m_mg2d.o m_gmres.o
test_mg_elec:		m_afivo_2d.o m_vtk.o m_morton.o m_mg2d.o m_gmres.o
test_gmres:		m_gmres.o
test_unstr_vtk:		m_vtk.o
test_omp:		test_mod_omp.o
