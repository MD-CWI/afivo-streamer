FC 	:= mpif90
FFLAGS	:= -Wall -O2 -fcheck=all -g -ffpe-trap=invalid,zero,overflow
OBJS	:= m_morton.f90 m_afivo.f90 m_vtk.f90

INCDIRS	:=
LIBDIRS :=
LIBS	:=

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all: test_morton test_base test_unstr_vtk test_diffusion

clean:
	$(RM) *.o *.mod test_m_morton

# Dependency information
m_afivo.o:		m_vtk.o m_morton.o
test_morton:		m_morton.o
test_base:		m_afivo.o m_vtk.o m_morton.o
test_diffusion:		m_afivo.o m_vtk.o m_morton.o
test_unstr_vtk:		m_vtk.o
