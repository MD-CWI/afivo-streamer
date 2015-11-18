INCDIRS	:= ../src
LIBDIRS := ../src ../silo/lib
LIBS	:= afivo silo

include ../makerules.make

STEPS	:= step_1 step_2 step_3 step_4 implicit

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all:	$(STEPS)

clean:
	$(RM) $(STEPS) *.o *.mod

# Dependency information
$(STEPS): 		../src/libafivo.a
