INCDIRS	:= ../afivo/src
LIBDIRS := ../afivo/src ../afivo/silo/lib .
LIBS	:= afivo silo streamer

include ../afivo/makerules.make

OBJS	:= m_units_constants.o m_config.o m_lookup_table.o\
	m_random.o m_merge_sort_ranking.o m_linked_list.o m_find_index.o \
	m_photons.o m_streamer.o m_geometry.o m_transport_data.o
PROGS	:= ../streamer_2d ../streamer_3d ../streamer_cyl

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

../%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all: 	$(PROGS) 

libstreamer.a: $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

clean:
	$(RM) $(PROGS) *.o *.mod libstreamer.a

# Dependency information
$(PROGS): 		../afivo/src/libafivo.a libstreamer.a
m_lookup_table.o:	m_find_index.o
m_streamer.o:		m_transport_data.o
