INCDIRS	:= ../src
LIBDIRS := ../src ../silo/lib
LIBS	:= afivo silo

include ../makerules.make

STEPS	:= step_1 step_2 step_3 step_4

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

step_1.f90: all_steps.f90
	sed 's|\$$STEP|1|g' < $< | cpp -traditional-cpp -P | cat -s > $@
step_2.f90: all_steps.f90
	sed 's|\$$STEP|2|g' < $< | cpp -traditional-cpp -P | cat -s > $@
step_3.f90: all_steps.f90
	sed 's|\$$STEP|3|g' < $< | cpp -traditional-cpp -P | cat -s > $@
step_4.f90: all_steps.f90
	sed 's|\$$STEP|4|g' < $< | cpp -traditional-cpp -P | cat -s > $@
