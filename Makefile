SRC_DIRS	:= src examples tests

# Directories with altered names (useful for cleaning)
CLEANSRC	:= $(SRC_DIRS:%=clean-%)

.PHONY:	all test doc clean $(SRC_DIRS) $(EXT_LIBS) $(CLEANSRC)

all:	$(SRC_DIRS)

doc:    $(SRC_DIRS)
	@doxygen

clean:	$(CLEANSRC)

$(SRC_DIRS):
	@$(MAKE) -C $@

$(CLEANSRC):
	$(MAKE) -C $(@:clean-%=%) clean

external_libraries/silo:
	@echo "Locally installing the Silo library"
	@cd external_libraries; ./build_silo.sh

# Dependency information
$(SRC_DIRS): external_libraries/silo
examples tests: src
