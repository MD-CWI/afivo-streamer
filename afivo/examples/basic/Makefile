A5_DIR  := ../..
INCDIRS	:= $(A5_DIR)/src
LIBDIRS := $(A5_DIR)/src $(A5_DIR)/silo/lib
LIBS	:= afivo silo

include $(A5_DIR)/makerules.make

SRCS	:= $(wildcard *.f90)
TESTS	:= $(patsubst %.f90, %, $(SRCS))

.PHONY: all clean

all:	$(TESTS)

clean:
	$(RM) $(TESTS) *.o *.mod

# Dependency information
$(TESTS): $(A5_DIR)/src/libafivo.a
