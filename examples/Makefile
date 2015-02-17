INCDIRS	:= ../src
LIBDIRS := ../src
LIBS	:= afivo

include ../makerules.make

TESTS	:= test_morton test_base_2d test_unstr_vtk test_drift_diff_2d \
	test_mg_2d test_base_3d test_mg2_2d test_mg2_3d test_streamer_2d \
	test_drift_diff_3d test_mgb_2d test_mg3_2d test_streamer_3d

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all:	$(TESTS)

clean:
	$(RM) $(TESTS) *.o *.mod

# Dependency information
$(TESTS): 		../src/libafivo.a
test_mg3_2d:		m_mg_diel.o
test_streamer_2d.o: 	m_mg_diel.o
test_streamer_2d: 	m_mg_diel.o
test_mg3_2d.o:		m_mg_diel.o
test_streamer_3d.o: 	m_mg_diel.o
test_streamer_3d: 	m_mg_diel.o
