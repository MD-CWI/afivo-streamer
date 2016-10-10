FC := gfortran
FFLAGS := -O2 -g -std=f2008 -Wall -Wextra
OBJS := m_find_index.o m_lookup_table.o
TESTS := test_find_index_simple test_find_index_performance	\
test_lookup_table_performance

.PHONY:	all clean

all: 	$(TESTS)

clean:
	$(RM) $(TESTS) $(OBJS) $(OBJS:.o=.mod)

# Dependency information
$(TESTS): m_lookup_table.o m_find_index.o
m_lookup_table.o: m_find_index.o

# How to get .o object files from .f90 source files
%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS)

# How to get executables from .o object files
%: %.o
	$(FC) -o $@ $^ $(FFLAGS)
