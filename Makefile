SRC_DIRS	:= src examples
CREATE_DIRS	:= silo

# Directories with altered names (useful for cleaning)
CLEANSRC	:= $(SRC_DIRS:%=clean-%)

# phonytest ensures that tests are always performed
.PHONY:	all test doc clean phonytest $(SRC_DIRS) $(EXT_LIBS) $(CLEANSRC)

all:	$(SRC_DIRS)

tests:	$(SRC_DIRS) phonytest
		@$(MAKE) -C $@

doc:
		$(MAKE) srcs -C src
		doxygen

clean:	$(CLEANSRC)

$(SRC_DIRS): | $(CREATE_DIRS)
		@echo "  *********** Build information ***********"
		@echo "  Debug is set to [$(DEBUG)], set it to 1 to enable a debug build."
		@echo "  For example: make clean; make DEBUG=1"
		@echo "  *****************************************"
		@$(MAKE) -C $@

$(CLEANSRC):
		$(MAKE) -C $(@:clean-%=%) clean

silo:
		./build_silo.sh

# Dependency information
$(SRC_DIRS):
