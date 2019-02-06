SRC_DIRS	:= lib_2d lib_3d examples tests

# Directories with altered names (useful for cleaning)
CLEANSRC	:= $(SRC_DIRS:%=clean-%)

.PHONY:	all test doc clean $(SRC_DIRS) $(EXT_LIBS) $(CLEANSRC)

all:	$(SRC_DIRS)

doc:    $(SRC_DIRS)
	@doxygen
	@echo "Done generating documentation, open it with e.g.:"
	@echo "firefox documentation/html/index.html"

clean:	$(CLEANSRC)

$(SRC_DIRS):
	@$(MAKE) -C $@

$(CLEANSRC):
	$(MAKE) -C $(@:clean-%=%) clean

external_libraries/silo:
	@echo "Locally installing the Silo library"
	@cd external_libraries; ./build_silo.sh

external_libraries/hypre:
	@echo "Locally installing the Hypre library"
	@cd external_libraries; ./build_hypre.sh

# Dependency information
$(SRC_DIRS): external_libraries/silo external_libraries/hypre
examples tests: lib_2d lib_3d
