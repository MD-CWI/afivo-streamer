INCDIRS	:= ../src
LIBDIRS := ../src ../silo/lib
LIBS	:= afivo silo

include ../makerules.make

TESTS	:= test_morton test_base_2d test_drift_diff_2d \
	test_mg_2d test_base_3d test_mg2_2d test_mg2_3d \
	test_drift_diff_3d test_mgb_2d test_mg3_2d \
	test_mg_cyl test_mg_cyl_diel test_reduction_2d \
	test_fractal test_mg_benchmark_2d

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
