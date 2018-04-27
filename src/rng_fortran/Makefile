FC := gfortran
FFLAGS := -O3 -flto -g -std=f2008 -Wall -Wextra -fopenmp
OBJS := m_random.o
TESTS := example performance_test parallel_test

.PHONY:	all clean

all: 	$(TESTS)

clean:
	$(RM) $(TESTS) $(OBJS) $(OBJS:.o=.mod)

# Dependency information
$(TESTS): $(OBJS)

# How to get .o object files from .f90 source files
%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS)

# How to get executables from .o object files
%: %.o
	$(FC) -o $@ $^ $(FFLAGS)
