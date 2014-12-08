SRC_DIRS	:= src tests

# Directories with altered names (useful for cleaning)
CLEANSRC	:= $(SRC_DIRS:%=clean-%)

.PHONY:	all clean $(SRC_DIRS) $(EXT_LIBS) $(CLEANSRC)

all: 		$(SRC_DIRS) | $(CREATE_DIRS)

clean: 		$(CLEANSRC)

$(SRC_DIRS):
		@echo "\n*********** Build information ***********"
		@echo "  Debug is set to: [$(DEBUG)],"
		@echo "  Set it to 1 to enable a debug build."
		@echo "  For example: make clean; make DEBUG=1"
		@echo "*****************************************\n"
		$(MAKE) -C $@
$(CLEANSRC):
		$(MAKE) -C $(@:clean-%=%) clean

# Dependecy information
$(SRC_DIRS):
