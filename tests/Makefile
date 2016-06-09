# Need these libraries for the tests
INCDIRS	:= ../afivo/src ../src
LIBDIRS := ../afivo/src ../afivo/silo/lib ../src .
LIBS	:= afivo silo streamer

include ../afivo/makerules.make

TESTS	:= test_m_config test_m_lookup_table test_m_random test_m_linked_list\
	test_m_find_index test_m_photons test_photoionization

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

# Dependency information
test_m_config.o:	../src/libstreamer.a
$(TESTS):		../src/libstreamer.a
