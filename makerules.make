# Disable built in rules
.SUFFIXES:

COMPILER = gfortran

ifeq ($(COMPILER), gfortran)
	FC 	:= gfortran
	FFLAGS	:= -O2 -std=f2008 -fopenmp -Wall -Wextra -Wimplicit-interface -Wshadow
	ifeq ($(DEBUG), 1)
		FFLAGS += -O0 -fcheck=all -g -pg -ffpe-trap=invalid,zero,overflow \
		-pedantic -finit-real=snan
	endif
	ifeq ($(PROF), 1)
		FFLAGS += -pg
	endif

else ifeq ($(COMPILER), ifort)
	FC 	:= ifort
	FFLAGS	:= -warn all -O2 -stand f08 -openmp -assume realloc-lhs
	ifeq ($(DEBUG), 1)
		FFLAGS += -check all,noarg_temp_created -g -p -debug all
	endif
endif

# How to get .o object files from .f90 source files
%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get executables from .o object files
%: %.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

# How to preprocess dimension-independent files to 2D
PPROCESS_2D = sed 's|$$D|2|g' < $< | cpp -traditional-cpp -P | cat -s > $@

# How to preprocess dimension-independent files to 3D
PPROCESS_3D = sed 's|$$D|3|g' < $< | cpp -traditional-cpp -P | cat -s > $@
