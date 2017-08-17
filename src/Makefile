INCDIRS	:= ../afivo/src
LIBDIRS := ../afivo/src ../afivo/external_libraries/silo/lib .
LIBS	:= afivo silo streamer
VPATH	:= config_fortran:lookup_table_fortran:rng_fortran

include ../afivo/makerules.make

OBJS := m_units_constants.o m_config.o m_lookup_table.o m_random.o		\
	m_find_index.o m_photoi_mc.o m_streamer.o m_geometry.o m_transport_data.o	\
	m_field_2d.o m_field_3d.o m_init_cond_2d.o m_init_cond_3d.o		\
	m_photoi_helmh_2d.o m_photoi_helmh_3d.o m_photoi_2d.o m_photoi_3d.o

PROGS	:= ../streamer_2d ../streamer_3d

.PHONY: all clean

all: 	$(PROGS)

libstreamer.a: $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

clean:
	$(RM) $(PROGS) *.o *.mod libstreamer.a

%_2d.f90: %_Xd.f90
	@echo "*** $< has changed, regenerating $@ ***"
	@$(PPROCESS_2D)

%_3d.f90: %_Xd.f90
	@echo "*** $< has changed, regenerating $@ ***"
	@$(PPROCESS_3D)

# Keep intermediate files for inspection
.PRECIOUS: %_2d.f90 %_3d.f90

# How to get executables from .o object files
../%: %.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

# Dependency information
$(PROGS): 		../afivo/src/libafivo.a libstreamer.a
streamer_2d.o streamer_3d.o: m_streamer.o m_geometry.o
streamer_2d.o: m_field_2d.o m_init_cond_2d.o m_photoi_2d.o
streamer_3d.o: m_field_3d.o m_init_cond_3d.o m_photoi_3d.o
m_field_2d.o m_init_cond_2d.o: m_streamer.o
m_field_3d.o m_init_cond_3d.o: m_streamer.o
m_photoi_helmh_2d.o m_photoi_helmh_3d.o: m_streamer.o
m_lookup_table.o:	m_find_index.o
m_streamer.o:		m_transport_data.o m_config.o m_random.o
m_photoi_mc.o:		m_lookup_table.o m_units_constants.o m_random.o m_streamer.o
m_photoi_2d.o: m_photoi_mc.o m_photoi_helmh_2d.o
m_photoi_3d.o: m_photoi_mc.o m_photoi_helmh_3d.o

# Hide some incorrect warnings
m_photoi_helmh_2d.o m_photoi_helmh_3d.o: FFLAGS += -Wno-unused-function
m_photoi_2d.o m_photoi_3d.o: FFLAGS += -Wno-unused-function

