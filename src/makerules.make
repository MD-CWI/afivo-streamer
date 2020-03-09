# Disable built in rules
.SUFFIXES:

COMPILER = gfortran

ifeq ($(COMPILER), gfortran)
	FC := gfortran
	FFLAGS := -O2 -g -std=f2008 -fopenmp -Wall -Wextra -Wimplicit-interface	\
	-Wshadow -Wno-unused-dummy-argument -cpp
	ifeq ($(DEBUG), 1)
		FFLAGS += -O0 -fcheck=all -pg				\
		-ffpe-trap=invalid,zero,overflow -pedantic -finit-real=snan
	endif
else ifeq ($(COMPILER), ifort)
	FC 	:= ifort
	FFLAGS := -warn all -O2 -g -stand f08 -openmp -assume realloc-lhs -fpp
	ifeq ($(DEBUG), 1)
		FFLAGS += -check all,noarg_temp_created -p -debug all
	endif
endif

ifeq ($(PROF), gprof)
  FFLAGS += -pg
else ifeq ($(PROF), gperftools)
  LIBS += profiler
endif

# How to get .o object files from .f90 source files
%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get .mod files from .f90 source files (remake only if they have been
# removed, otherwise assume they are up to date)
%.mod: %.f90 %.o
	@test -f $@ || $(FC) -c -o $(@:.mod=.o) $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get executables from .o object files
%: %.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
