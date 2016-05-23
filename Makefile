SRC_DIRS	:= src examples
CREATE_DIRS	:= silo

# Directories with altered names (useful for cleaning)
CLEANSRC	:= $(SRC_DIRS:%=clean-%)

# Hide entering/leaving directory messages
MAKEFLAGS	:= --no-print-directory

.PHONY:	all test doc clean $(SRC_DIRS) $(EXT_LIBS) $(CLEANSRC)

all:	$(SRC_DIRS)

tests:	$(SRC_DIRS)
		@$(MAKE) -C $@

doc:
		$(MAKE) srcs -C src
		doxygen

clean:	$(CLEANSRC)

$(SRC_DIRS): | $(CREATE_DIRS)
		@$(MAKE) -C $@

$(CLEANSRC):
		$(MAKE) -C $(@:clean-%=%) clean

silo:
		./build_silo.sh

# Dependency information
$(SRC_DIRS):
