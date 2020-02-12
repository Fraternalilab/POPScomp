#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_INIT(pops, 3.1.1, jens@jkleinj.eu)

AC_PREREQ([2.63])
AC_CONFIG_SRCDIR([src/arg.c])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

# Checks for programs.
AC_PROG_CC_C_O
AC_PROG_CXX
AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB([m], [main], [],
				AC_MSG_ERROR([Could not find required C math library.]))

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([pow sqrt])

# Checks for OpenMP
AC_OPENMP

# enable OpenMP
sinclude(ax_openmp.m4)
AC_ARG_ENABLE(openmp, [  --enable-openmp         Turn on OpenMP],
                [ac_use_openmp="yes"], [ac_use_openmp="no"])

if test "$ac_use_openmp" = "yes"; then
   AX_OPENMP([], AC_MSG_ERROR([Could not find required OpenMP library.]))
   CFLAGS="$CFLAGS $OPENMP_CFLAGS"
   LIBS="$LIBS $OPENMP_LIBS"
   AC_DEFINE(OPENMP,1,[Compiling with OpenMP settings.])
fi

# enable debug
AC_ARG_ENABLE(debug, [  --enable-debug  Turn on debugging], [ac_use_debug="yes"], [ac_use_debug="no"])
if test "$ac_use_debug" = "yes"; then
    CFLAGS="$CFLAGS -ggdb"
    AC_DEFINE(DEBUG,1,[Define to compile with DEBUG settings.])
else
    AC_DEFINE(DEBUG,0,[Define to compile without DEBUG settings.])
fi

# enable profiling 
AC_ARG_ENABLE(profiling, [  --enable-profiling  Turn on profiling], [ac_use_profiling="yes"], [ac_use_profiling="no"])
if test "$ac_use_profiling" = "yes"; then
    CFLAGS="$CFLAGS -pg -O2"
    AC_DEFINE(PROFILING,1,[Define to compile with PROFILING settings.])
fi

# check for libxml2
PKG_CHECK_MODULES(XML, libxml-2.0 >= 2.4)

# check for zlib 
PKG_CHECK_MODULES([ZLIB], zlib >= 1.0)

# output files
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_CONFIG_HEADERS([src/config.h])
AC_OUTPUT

