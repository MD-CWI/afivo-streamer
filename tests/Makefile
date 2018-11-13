AF_DIR  := ..
LIBDIRS := $(AF_DIR)/external_libraries/silo/lib
LIBS	:= afivo silo

include $(AF_DIR)/src/makerules.make

# All program sources (the .f90 files except modules)
SRCS := $(filter-out $(wildcard m_*.f90), $(wildcard *.f90))
TESTS_2D := $(SRCS:%.f90=%_2d)
TESTS_3D := $(SRCS:%.f90=%_3d)
TESTS := $(TESTS_2D) $(TESTS_3D)

ANSDIR	:= answers
RESDIR	:= results
RESULTS	:= $(TESTS:%=$(RESDIR)/%)
ANSWERS	:= $(TESTS:%=$(ANSDIR)/%)

# phonytest helps to always regenerate test results
.PHONY: all test score clean phonytest

all:	$(TESTS)

test:	$(RESULTS)

score: 	$(TESTS)
	@$(MAKE) $(RESULTS) | cut -d\  -f 1 | egrep '(PASSED|FAILED)' | sort | uniq -c

# Rule for generating (correct) answers to tests
$(ANSWERS):	$(TESTS) | $(ANSDIR)
	@echo Storing answer $@
	@./$(@:$(ANSDIR)/%=%) > $@

# Rule to generate results and compare to answer
$(RESULTS): $(TESTS) phonytest | $(RESDIR)
	@./$(@:$(RESDIR)/%=%) > $@
	@if cmp -s $@ $(@:$(RESDIR)/%=$(ANSDIR)/%) ; \
	then echo PASSED $(@:$(RESDIR)/%=%) ; \
	else echo FAILED $(@:$(RESDIR)/%=%) ; \
	fi

clean:
	$(RM) $(TESTS) $(RESULTS) *.o *.mod

$(RESDIR) $(ANSDIR):
	@mkdir -p $@

# How to create executables
%_2d: %.f90
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -I,$(INCDIRS)) \
	$(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
%_3d: %.f90
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -I,$(INCDIRS)) \
	$(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

$(TESTS_2D): $(AF_DIR)/lib_2d/libafivo.a | $(OUTDIR)
$(TESTS_2D): INCDIRS+=$(AF_DIR)/lib_2d
$(TESTS_2D): LIBDIRS+=$(AF_DIR)/lib_2d
$(TESTS_2D): FFLAGS += -DNDIM=2

$(TESTS_3D): $(AF_DIR)/lib_3d/libafivo.a | $(OUTDIR)
$(TESTS_3D): INCDIRS+=$(AF_DIR)/lib_3d
$(TESTS_3D): LIBDIRS+=$(AF_DIR)/lib_3d
$(TESTS_3D): FFLAGS += -DNDIM=3
