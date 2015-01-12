FC 	:= gfortran
FFLAGS	:= -Wall -O2 -std=f2008 -fopenmp -cpp

ifeq ($(DEBUG), 1)
	FFLAGS += -fcheck=all -g -pg -ffpe-trap=invalid,zero,overflow \
	-pedantic -finit-real=nan
endif

SRCS	:= m_morton.f90 m_afivo_constants.f90 m_afivo_2d.f90 \
	m_afivo_3d.f90 m_vtk.f90 m_mg_2d.f90 m_mg_3d.f90
OBJS	:= $(SRCS:%.f90=%.o)

INCDIRS	:=
LIBDIRS :=
LIBS	:=

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

.PHONY: srcs clean

libafivo.a: $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

srcs:	$(SRCS)

clean:
	$(RM) libafivo.a *.o *.mod m_afivo_2d.f90 m_afivo_3d.f90 \
	m_mg_2d.f90 m_mg_3d.f90

m_afivo_2d.f90: m_afivo_Xd.f90
	sed 's|$$D|2|g' < $< | cpp -traditional-cpp -P | cat -s > $@

m_afivo_3d.f90: m_afivo_Xd.f90
	sed 's|$$D|3|g' < $< | cpp -traditional-cpp -P | cat -s > $@

m_mg_2d.f90: m_mg_Xd.f90
	sed 's|$$D|2|g' < $< | cpp -traditional-cpp -P | cat -s > $@

m_mg_3d.f90: m_mg_Xd.f90
	sed 's|$$D|3|g' < $< | cpp -traditional-cpp -P | cat -s > $@

# Dependency information
m_afivo_2d.o:		m_vtk.o m_morton.o m_afivo_constants.o
m_afivo_3d.o:		m_vtk.o m_morton.o m_afivo_constants.o
m_mg_2d.o:		m_afivo_2d.o
m_mg_3d.o:		m_afivo_3d.o
