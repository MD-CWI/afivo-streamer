A5_DIR  := ..
INCDIRS	:= $(A5_DIR)/src
LIBDIRS := $(A5_DIR)/src $(A5_DIR)/silo/lib
LIBS	:= afivo silo

include $(A5_DIR)/makerules.make

# All program sources (the .f90 files except modules)
SRCS	:= $(filter-out $(wildcard m_*.f90), $(wildcard *.f90))
PROGS	:= $(SRCS:%.f90=%)

.PHONY: all clean

all:	$(PROGS)

clean:
	$(RM) $(PROGS) *.o *.mod

$(OUTDIR):
		mkdir -p $@

# Dependency information
$(PROGS): $(A5_DIR)/src/libafivo.a
