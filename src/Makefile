FC 	:= gfortran
FFLAGS	:= -Wall -O2 -fcheck=all -g -pg -ffpe-trap=invalid,zero,overflow\
	-pedantic -std=f2008
# FFLAGS	:= -Ofast
OBJS	:= m_morton.f90 m_afivo.f90 m_vtk.f90 m_mg2d.f90 m_gmres.f90

INCDIRS	:=
LIBDIRS :=
LIBS	:=

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all: test_morton test_base test_unstr_vtk test_drift_diff test_mg test_gmres

clean:
	$(RM) *.o *.mod test_m_morton

# Dependency information
m_afivo.o:		m_vtk.o m_morton.o
m_mg2d.o:		m_afivo.o
test_morton:		m_morton.o
test_base:		m_afivo.o m_vtk.o m_morton.o
test_drift_diff:	m_afivo.o m_vtk.o m_morton.o
test_mg:		m_afivo.o m_vtk.o m_morton.o m_mg2d.o m_gmres.o
test_gmres:		m_gmres.o
test_unstr_vtk:		m_vtk.o
