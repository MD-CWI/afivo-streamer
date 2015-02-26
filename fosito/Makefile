FC 	:= gfortran
FFLAGS	:= -Wall -O2
OBJS	:= m_units_constants.o m_config.o m_lookup_table.o\
	m_random.o m_mrgrnk.o m_linked_list.o m_find_index.o

ifeq ($(DEBUG), 1)
	FFLAGS += -fcheck=array-temps,bounds,do,mem,pointer\
	-g -ffpe-trap=invalid,zero,overflow
endif

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS)

.PHONY: all test clean

all: 	libfosito.a

libfosito.a: $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

test: 	libfosito.a
	$(MAKE) -C test

clean:
	$(RM) -f *.o *.mod
	$(MAKE) -C test clean

# Dependencies
m_lookup_table.o:	m_find_index.o