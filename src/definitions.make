OBJS := m_units_constants.o m_config.o m_lookup_table.o m_random.o		\
	m_find_index.o m_photoi_mc.o m_streamer.o m_geometry.o			\
	m_transport_data.o m_field.o m_init_cond.o m_photoi_helmh.o m_photoi.o	\
	m_chemistry.o

# Dependency information
streamer.o: m_streamer.o m_geometry.o
streamer.o: m_field.o m_init_cond.o m_photoi.o
m_field.o m_init_cond.o: m_streamer.o m_units_constants.o
m_photoi_helmh.o: m_streamer.o m_units_constants.o
m_lookup_table.o:	m_find_index.o
m_streamer.o:		m_transport_data.o m_lookup_table.o m_config.o m_random.o
m_photoi_mc.o: m_lookup_table.o m_units_constants.o m_random.o m_streamer.o
m_photoi.o: m_photoi_mc.o m_photoi_helmh.o
m_chemistry.o: m_streamer.o

# Hide some incorrect warnings
m_photoi_helmh.o: FFLAGS += -Wno-unused-function
m_photoi.o: FFLAGS += -Wno-unused-function

