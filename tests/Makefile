A5_DIR  := ..
INCDIRS	:= $(A5_DIR)/src
LIBDIRS := $(A5_DIR)/src $(A5_DIR)/silo/lib
LIBS	:= afivo silo

include $(A5_DIR)/makerules.make

# All program sources (the .f90 files except modules)
SRCS	:= $(filter-out $(wildcard m_*.f90), $(wildcard *.f90))
TESTS	:= $(SRCS:%.f90=%)

ANSDIR	:= answers
RESDIR	:= results
RESULTS	:= $(TESTS:%=$(RESDIR)/%)
ANSWERS	:= $(TESTS:%=$(ANSDIR)/%)

# phonytest helps to always regenerate test results
.PHONY: all score clean phonytest

all:	$(RESULTS)

score: 	$(TESTS) | $(ANSWERS)
		@$(MAKE) $(RESULTS) | cut -d\  -f 1 | egrep '(PASSED|FAILED)' | sort | uniq -c

# Rule for generating (correct) answers to tests
$(ANSWERS):	$(TESTS) | $(ANSDIR)
		@echo Storing answer $@
		@$(PWD)/$(@:$(ANSDIR)/%=%) > $@

# Rule to generate results and compare to answer
$(RESULTS): $(TESTS) phonytest | $(RESDIR) $(ANSWERS)
		@$(PWD)/$(@:$(RESDIR)/%=%) > $@
		@if cmp -s $@ $(@:$(RESDIR)/%=$(ANSDIR)/%) ; \
        then echo PASSED $(@:$(RESDIR)/%=%) ; \
        else echo FAILED $(@:$(RESDIR)/%=%) ; \
        fi

clean:
	$(RM) $(TESTS) $(RESULTS) *.o *.mod

$(RESDIR) $(ANSDIR):
		@mkdir -p $@

# Dependency information
$(TESTS): $(A5_DIR)/src/libafivo.a
