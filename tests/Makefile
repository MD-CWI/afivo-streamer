# Need these libraries for the tests
INCDIRS	:= ../afivo/src ../src
LIBDIRS := ../afivo/src ../afivo/external_libraries/silo/lib ../src .
LIBS	:= afivo silo streamer

include ../afivo/makerules.make

TESTS	:= test_photoionization

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
