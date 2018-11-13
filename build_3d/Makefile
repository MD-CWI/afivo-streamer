NDIM := 3
AF_DIR := ../afivo
INCDIRS := $(AF_DIR)/lib_$(NDIM)d
LIBDIRS := $(AF_DIR)/lib_$(NDIM)d $(AF_DIR)/external_libraries/silo/lib .
LIBS := afivo silo streamer
PROG := ../streamer_$(NDIM)d

.PHONY: all clean
all: $(PROG)

vpath %.f90 ../src ../src/config_fortran ../src/lookup_table_fortran ../src/rng_fortran
include ../src/definitions.make
include $(AF_DIR)/src/makerules.make

$(OBJS): FFLAGS += -DNDIM=$(NDIM)
$(OBJS): $(AF_DIR)/lib_$(NDIM)d/libafivo.a
$(PROG): FFLAGS += -DNDIM=$(NDIM)
$(PROG): libstreamer.a

$(PROG): streamer.f90
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -I,$(INCDIRS)) \
	$(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

libstreamer.a: $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

$(LIB): $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

clean:
	$(RM) $(PROG) *.o *.mod libstreamer.a

