AFIVO_DIR := afivo
2D_PROGS := bin_2d
3D_PROGS := bin_3d
SRC_DIRS := lib_2d lib_3d $(AFIVO_DIR)/lib_2d $(AFIVO_DIR)/lib_3d $(2D_PROGS)	\
$(3D_PROGS)

# Directories with altered names (useful for cleaning)
CLEANSRC := $(SRC_DIRS:%=clean-%)

.PHONY:	all 2d 3d doc clean $(SRC_DIRS) $(CLEANSRC)

all: 		2d 3d
2d:		$(2D_PROGS)
3d:		$(3D_PROGS)

doc:    	$(SRC_DIRS)
		@doxygen

clean: 		$(CLEANSRC)

$(SRC_DIRS):	| output
		$(MAKE) -C $@

$(CLEANSRC):
		$(MAKE) -C $(@:clean-%=%) clean

# Ensure the output folder exists
output:
		mkdir -p output

# Dependecy information
lib_2d: $(AFIVO_DIR)/lib_2d
lib_3d: $(AFIVO_DIR)/lib_3d
$(2D_PROGS): lib_2d
$(3D_PROGS): lib_3d


