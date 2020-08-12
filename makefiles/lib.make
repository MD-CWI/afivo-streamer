# This is a template to build afivo-streamer libraries in 1D/2D/3D

ifndef NDIM
$(error NDIM is not set)
endif

ifndef MAIN_DIR
$(error MAIN_DIR is not set)
endif

AFIVO_DIR := $(MAIN_DIR)/afivo
LIB := libstreamer.a

 # Used in the compilation rules
INCDIRS := $(AFIVO_DIR)/lib_$(NDIM)d

.PHONY: all clean allclean always_recompile
all: $(LIB)

# Where to find the source files
vpath %.f90 $(MAIN_DIR)/src $(MAIN_DIR)/src/config_fortran \
	$(MAIN_DIR)/src/lookup_table_fortran $(MAIN_DIR)/src/rng_fortran

# Dependencies
include $(MAIN_DIR)/src/definitions.make

# Compilation rules
include  $(MAIN_DIR)/makefiles/makerules.make

FFLAGS += -DNDIM=$(NDIM)

$(OBJS): $(AFIVO_DIR)/lib_$(NDIM)d/libafivo.a

$(LIB): $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

clean:
	$(RM) *.o *.mod $(LIB)

allclean: clean
	$(MAKE) -C $(AFIVO_DIR) clean

$(LIB): $(AFIVO_DIR)/lib_$(NDIM)d/libafivo.a

$(AFIVO_DIR)/lib_$(NDIM)d/libafivo.a: always_recompile
	$(MAKE) -C $(AFIVO_DIR) lib_$(NDIM)d
