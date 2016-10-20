SRC_DIRS	:= src afivo

# Directories with altered names (useful for cleaning)
CLEANSRC	:= $(SRC_DIRS:%=clean-%)

.PHONY:	all doc clean $(SRC_DIRS) $(CLEANSRC)

all: 		$(SRC_DIRS)

doc:    	$(SRC_DIRS)
		@doxygen

clean: 		$(CLEANSRC)

$(SRC_DIRS):	| output
		@echo "\n*********** Build information ***********"
		@echo "  Debug is set to: [$(DEBUG)],"
		@echo "  Set it to 1 to enable a debug build."
		@echo "  For example: make clean; make DEBUG=1"
		@echo "*****************************************\n"
		$(MAKE) -C $@

$(CLEANSRC):
		$(MAKE) -C $(@:clean-%=%) clean

# Ensure the output folder exists
output:
		mkdir -p output

# Dependecy information
src:		afivo
