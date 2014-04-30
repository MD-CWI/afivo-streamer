FC 	:= gfortran
FFLAGS	:= -Wall -fcheck=all -ffpe-trap=invalid,zero,overflow -g -O3
OBJS	:= m_units_constants.o m_config.o m_lookup_table.o\
	m_random.o m_mrgrnk.o
TESTS	:= test_m_config test_m_lookup_table test_m_random

LIBS	:= fosito

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS)

%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -l,$(LIBS))

.PHONY: all test clean

all: 	libfosito.a

libfosito.a: $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

$(TESTS): libfosito.a

test: 	$(TESTS)
	$(foreach test, $(TESTS), ./$(test);)

clean:
	$(RM) -f *.o *.mod $(TESTS)
