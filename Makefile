INCDIRS	:= ../src
LIBDIRS := ../src ../silo/lib
LIBS	:= afivo silo

include ../makerules.make

PROGS	:= streamer_2d

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all:	$(PROGS)

clean:
	$(RM) $(PROGS) *.o *.mod

# Dependency information
$(PROGS): 		../src/libafivo.a
