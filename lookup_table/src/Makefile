FC 	:= gfortran
FFLAGS	:= -O2 -g -std=f2008 -Wall -Wextra

OBJS	:= m_find_index.o m_lookup_table.o 

.PHONY:	all test clean

all: 	test_lookup_table

clean:
	$(RM) test_lookup_table m_lookup_table.o m_lookup_table.mod \
		m_find_index.o m_find_index.mod

# Dependency information
test_lookup_table:	$(OBJS)
m_lookup_table.o:	m_find_index.o


# How to get .o object files from .f90 source files
%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS)

# How to get executables from .o object files
%: %.o
	$(FC) -o $@ $^ $(FFLAGS)
