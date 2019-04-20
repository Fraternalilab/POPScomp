#===============================================================================
# Makefile.am : for automake
# Copyright (C) 2009 Jens Kleinjung and Franca Fraternali
#===============================================================================

bin_PROGRAMS = pops

AM_CPPFLAGS = $(INTI_CFLAGS) $(XML_CFLAGS) $(ZLIB_CFLAGS)
AM_CFLAGS = -Wall

pops_SOURCES = \
arg.c arg.h argpdb.h atomgroup.c atomgroup.h atomgroup_data.h \
cJSON.c cJSON.h \
error.c error.h fs.c fs.h \
getpdb.c getpdb.h getpdbml.c getpdbml.h gettraj.c gettraj.h \
json.c json.h \
matrix.c matrix.h modstring.c modstring.h \
pattern.c pattern.h pdb_structure.h pops.c pops.h \
putsasa.c putsasa.h putbsasa.c putbsasa.h putsigma.c putsigma.h \
safe.c safe.h sasa.c sasa.h sasa_const.h sasa_data.h seq.h \
sigma.c sigma.h sigma_const.h sigma_data.h \
topol.c topol.h type.c type.h vector.c vector.h

pops_LDADD = $(INTI_LIBS) $(XML_LIBS) $(ZLIB_LIBS)

EXTRA_DIST = doxygen.cfg

CLEANFILES = $(TARGET) *.o 
DISTCLEANFILES = libtool config.cache config.log

