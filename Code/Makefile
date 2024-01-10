FILES=Makefile IPFEMPB.c IPFEMPB.h coefficient.c coefficient.h functions_aly.c functions_aly.h functions_general.c functions_general.h \
functions_LPBE.c functions_LPBE.h functions_NPBE.c functions_NPBE.h options.c options.h unit.h
OBJ=IPFEMPB.o coefficient.o functions_aly.o functions_general.o functions_LPBE.o functions_NPBE.o options.o

default: IPFEMPB

all: IPFEMPB

IPFEMPB:$(OBJ)
	$(LINKER) ${LDFLAGS} -o $@ ${OBJ} ${LIBS}

clean:~
	-/bin/rm -f *.o IPFEMPB

distclean: clean
	-/bin/rm -f $(OBJ) IPFEMPB

#include /***/phg-0.9.7-20220927/share/phg/Makefile.inc

.PHONY: default all clean distclean run git
