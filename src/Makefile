FC 	:= mpif90
FFLAGS	:= -Wall -O2 -fcheck=all -g -ffpe-trap=invalid,zero,overflow
OBJS	:= m_morton.f90 m_write_silo.f90 m_afivo.f90 m_vtk.f90

INCDIRS	:= ../silo/include $(HOME)/opt/src/Lib_VTK_IO/mod
LIBDIRS := ../silo/lib $(HOME)/opt/src/Lib_VTK_IO/static
LIBS	:= silo vtkio

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all: test_morton test_silo test_base test_unstr_vtk

clean:
	$(RM) *.o *.mod test_m_morton

# Dependency information
test_morton:		m_morton.o
test_silo:		m_write_silo.o
test_base:		m_afivo.o m_write_silo.o
test_unstr_vtk:		m_vtk.o
