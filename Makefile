# put common definitions in here

TARGET1DMI      = stat-mi-1d
TARGET2DMI      = stat-mi-2d
TARGET3DMI      = stat-mi-3d
TARGET1DCMI     = stat-cmi-1d
TARGET2DCMI     = stat-cmi-2d
TARGET3DCMI      = stat-cmi-3d

CONFIGFILE      = 

BASEDIR         = $(HOME)/opt
CPP         	= g++
#CPPFLAGS        = -g3
CPPFLAGS        = -O2
FORTRAN         = gfortran
CC              = gcc
LD              = g++
CTT             = ${BASEDIR}/bin/ctt.py
LDFLAGS         = -L${BASEDIR}/lib -lm -ldranxor
RM              = rm
SHELL	        = /bin/bash
INCLUDES        = -I${BASEDIR}/include
EXEEXT          = .p4
EXEPREFIX       = ctx-
OBJS1DCMI       = digamma.o base1d.o main-cmi-1d.o cmutual.o utils.o
OBJS2DCMI       = digamma.o base2d.o main-cmi-2d.o cmutual.o utils.o
OBJS3DCMI       = digamma.o base3d.o main-cmi-3d.o cmutual.o utils.o
OBJS1DMI        = digamma.o base1d.o main-mi-1d.o mutual.o utils.o
OBJS2DMI        = digamma.o base2d.o main-mi-2d.o mutual.o utils.o
OBJS3DMI        = digamma.o base3d.o main-mi-3d.o mutual.o utils.o



.PHONY: clean all install


# all:  ${EXEPREFIX}${TARGET1DMI}${EXEEXT} ${EXEPREFIX}${TARGET2DMI}${EXEEXT} ${EXEPREFIX}${TARGET3DMI}${EXEEXT} ${EXEPREFIX}${TARGET1DCMI}${EXEEXT} ${EXEPREFIX}${TARGET2DCMI}${EXEEXT} ${EXEPREFIX}${TARGET3DCMI}${EXEEXT}
all:  ${EXEPREFIX}${TARGET1DMI}${EXEEXT} ${EXEPREFIX}${TARGET1DCMI}${EXEEXT} 

install:
# 	install ${EXEPREFIX}${TARGET2DMI}${EXEEXT} ${BASEDIR}/bin
	install ${EXEPREFIX}${TARGET1DMI}${EXEEXT} ${BASEDIR}/bin
# 	install ${EXEPREFIX}${TARGET2DCMI}${EXEEXT} ${BASEDIR}/bin
	install ${EXEPREFIX}${TARGET1DCMI}${EXEEXT} ${BASEDIR}/bin


${EXEPREFIX}${TARGET1DMI}${EXEEXT}: ${OBJS1DMI}
	${LD} ${OBJS1DMI} -o $@ ${CPPFLAGS} ${LDFLAGS}

${EXEPREFIX}${TARGET2DMI}${EXEEXT}: ${OBJS2DMI}
	${LD} ${OBJS2DMI} -o $@ ${CPPFLAGS} ${LDFLAGS}

${EXEPREFIX}${TARGET3DMI}${EXEEXT}: ${OBJS3DMI}
	${LD} ${OBJS3DMI} -o $@ ${CPPFLAGS} ${LDFLAGS}

${EXEPREFIX}${TARGET1DCMI}${EXEEXT}: ${OBJS1DCMI}
	${LD} ${OBJS1DCMI} -o $@ ${CPPFLAGS} ${LDFLAGS}

${EXEPREFIX}${TARGET2DCMI}${EXEEXT}: ${OBJS2DCMI}
	${LD} ${OBJS2DCMI} -o $@ ${CPPFLAGS} ${LDFLAGS}

${EXEPREFIX}${TARGET3DCMI}${EXEEXT}: ${OBJS3DCMI}
	${LD} ${OBJS3DCMI} -o $@ ${CPPFLAGS} ${LDFLAGS}


%.o : %.cxx
	${CPP} ${INCLUDES} ${CPPFLAGS} -c $<

%.cxx : %.ct
	${CTT} -n -i -u $<

clean:
	${RM} -f *.o ${EXEPREFIX}*${EXEEXT}

