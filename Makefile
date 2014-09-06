SRC_DIRS	:= src
EXT_DEPS	:= ext_libs/silo

# Directories with altered names (useful for cleaning)
CLEANSRC	:= $(SRC_DIRS:%=clean-%)

.PHONY:	all clean $(SRC_DIRS) $(EXT_LIBS) $(CLEANSRC)

all: 		$(SRC_DIRS) | $(CREATE_DIRS)

clean: 		$(CLEANSRC)

$(SRC_DIRS):
		$(MAKE) -C $@
$(CLEANSRC):
		$(MAKE) -C $(@:clean-%=%) clean

ext_libs/silo:
	cd ext_libs && bash build_silo.sh

# Dependecy information
$(SRC_DIRS):	| EXT_DEPS
