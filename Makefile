2D_PROGS := bin_2d
3D_PROGS := bin_3d
SRC_DIRS := $(2D_PROGS) $(3D_PROGS)

# Directories with altered names (useful for cleaning)
CLEANSRC := $(SRC_DIRS:%=clean-%) clean-lib_2d clean-lib_3d clean-afivo/lib_2d	\
clean-afivo/lib_3d

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


