FC 	:= gfortran
FFLAGS	:= -Wall -fcheck=all -ffpe-trap=invalid,zero,overflow -g -O0
TESTS	:= test_m_config test_m_lookup_table test_m_random test_m_linked_list\
	test_m_find_index

# Need these libraries for the tests
BASEDIR	:= ../../
LIBS	:= fosito
INCDIRS := $(addprefix $(BASEDIR), $(LIBS))
LIBDIRS := $(INCDIRS)

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all test clean

all: 	$(TESTS)

test:	$(TESTS)
	$(foreach test, $(TESTS), ./$(test);)

clean:
	$(RM) *.o *.mod $(TESTS)
