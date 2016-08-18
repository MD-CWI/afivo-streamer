A5_DIR  := ..
INCDIRS	:= $(A5_DIR)/src
LIBDIRS := $(A5_DIR)/src $(A5_DIR)/external_libraries/silo/lib
LIBS	:= afivo silo
OUTDIR	:= output

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
$(PROGS): $(A5_DIR)/src/libafivo.a | $(OUTDIR)

# Poisson programs and objects depend on m_gaussians
PROGS_POIS	:= $(filter poisson_%, $(PROGS))
OBJS_POIS	:= $(PROGS_POIS:%=%.o)
$(PROGS_POIS): 	m_gaussians.o
$(OBJS_POIS):	m_gaussians.o
