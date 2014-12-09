FC 	:= gfortran
FFLAGS	:= -Wall -O2 -std=f2008 -fopenmp -cpp

ifeq ($(DEBUG), 1)
	FFLAGS += -fcheck=all -g -pg -ffpe-trap=invalid,zero,overflow \
	-pedantic -finit-real=nan
endif

INCDIRS	:= ../src
LIBDIRS := ../src
LIBS	:= afivo

TESTS	:= test_morton test_base test_unstr_vtk test_drift_diff test_mg \
	test_base_3d test_mg2 test_streamer_2d

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all:	$(TESTS)

clean:
	$(RM) $(TESTS) *.o *.mod

# Dependency information
