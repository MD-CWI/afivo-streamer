SRC_ND	:= m_af_t.f90 m_morton.f90 m_vtk.f90 m_write_silo.f90
SRC_XD	:= m_aX_t.f90 m_aX_core.f90 m_aX_io.f90 m_aX_gc.f90 \
		m_aX_restrict.f90 m_aX_prolong.f90 m_aX_utils.f90 \
		m_aX_mg.f90 m_aX_interp.f90

SRC_2D	:= $(SRC_XD:m_aX_%.f90=m_a2_%.f90)
SRC_3D	:= $(SRC_XD:m_aX_%.f90=m_a3_%.f90)
SRCS	:= $(SRC_ND) $(SRC_2D) $(SRC_3D)
OBJS	:= $(SRCS:%.f90=%.o)

INCDIRS	:=
LIBDIRS :=
LIBS	:=

include ../makerules.make

.PHONY: srcs clean allclean

libafivo.a: $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

srcs:	$(SRCS)

clean:
	$(RM) libafivo.a *.o *.mod

allclean: clean
	$(RM) $(SRC_2D) $(SRC_3D)

m_a2_%.f90: m_aX_%.f90
	@echo "*** $< has changed, regenerating $@ ***"
	@sed 's|$$D|2|g' < $< | cpp -traditional-cpp -P | cat -s > $@

m_a3_%.f90: m_aX_%.f90
	@echo "*** $< has changed, regenerating $@ ***"
	@sed 's|$$D|3|g' < $< | cpp -traditional-cpp -P | cat -s > $@

# Dependency information
m_a2_core.o:		m_morton.o m_a2_t.o
m_a2_io.o:			m_vtk.o m_write_silo.o m_a2_t.o m_a2_interp.o
m_a2_gc.o:			m_a2_t.o m_a2_prolong.o
m_a2_restrict.o:	m_a2_t.o
m_a2_prolong.o:		m_a2_t.o
m_a2_interp.o:      m_a2_utils.o
m_a2_utils.o:		m_a2_t.o
m_a2_mg.o:			m_a2_t.o m_a2_utils.o
m_a2_t.o:			m_af_t.o

m_a3_core.o:		m_morton.o m_a3_t.o
m_a3_io.o:			m_vtk.o m_write_silo.o m_a3_t.o m_a3_interp.o
m_a3_gc.o:			m_a3_t.o m_a3_prolong.o
m_a3_restrict.o:	m_a3_t.o
m_a3_prolong.o:		m_a3_t.o
m_a3_interp.o:      m_a3_utils.o
m_a3_utils.o:		m_a3_t.o
m_a3_mg:			m_a3_t.o m_a3_utils.o
m_a3_t.o:			m_af_t.o

m_write_silo.o:		INCDIRS += ../silo/include
m_write_silo.o:		FFLAGS += -Wno-implicit-interface
