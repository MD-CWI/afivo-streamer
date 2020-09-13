NDIM := 2
LIB := libafivo.a

.PHONY: all clean
all: $(LIB)

vpath %.f90 ../src
include ../src/definitions.make
include ../src/makerules.make

$(OBJS): FFLAGS += -DNDIM=$(NDIM)

$(LIB): $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

clean:
	$(RM) $(LIB) *.o *.mod
