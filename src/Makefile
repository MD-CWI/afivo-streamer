INCDIRS	:= ../afivo/src ../fosito
LIBDIRS := ../afivo/src ../afivo/silo/lib ../fosito
LIBS	:= afivo silo fosito

include ../afivo/makerules.make

PROGS	:= ../streamer_2d ../streamer_3d ../streamer_cyl\
	../guiding_3d test_m_photons test_photoionization

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
$(PROGS): 			../afivo/src/libafivo.a ../fosito/libfosito.a
test_m_photons:		m_photons.o
test_photoionization:	m_photons.o
../streamer_2d: 	m_transport_data.o m_photons.o m_geom.o m_streamer.o
../streamer_3d: 	m_transport_data.o m_photons.o m_geom.o m_streamer.o
../guiding_3d: 		m_transport_data.o m_photons.o m_geom.o m_streamer.o
../streamer_cyl: 	m_transport_data.o m_photons.o m_geom.o m_streamer.o
m_streamer.o:		m_transport_data.o
