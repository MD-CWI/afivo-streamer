INCDIRS	:= ../afivo/src ../fosito
LIBDIRS := ../afivo/src ../afivo/silo/lib ../fosito
LIBS	:= afivo silo fosito

include ../afivo/makerules.make

PROGS	:= ../streamer_2d test_m_photons

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

../%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all:	$(PROGS)

clean:
	$(RM) $(PROGS) *.o *.mod

# Dependency information
$(PROGS): ../afivo/src/libafivo.a ../fosito/libfosito.a
test_m_photons:	m_photons.o
