# List with all the program directories
PROGRAM_DIRS := programs/standard_2d programs/standard_3d		\
programs/dielectric_2d programs/stability_cyl programs/2d_sprite	\
programs/3d_sprite

.PHONY: all $(PROGRAM_DIRS)

# We should not build multiple programs in parallel
.NOTPARALLEL: $(PROGRAM_DIRS)

all: $(PROGRAM_DIRS)

$(PROGRAM_DIRS):
	@$(MAKE) -C $@
