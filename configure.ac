#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([PATMOL],[2.8],[patjnichols@gmail.com])
AC_CONFIG_SRCDIR([basis.c])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE([1.9 foreign])
# Checks for programs.

AC_PROG_CXX([mpic++ mpicxx g++])
AC_PROG_CC([mpicc gcc])
AX_PROG_CC_MPI
AC_LANG([C])
# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdint.h stdlib.h string.h sys/time.h])
AC_CHECK_HEADER([mpi.h],[main],[],AC_MSG_ERROR([could find mpi header]))
AC_CHECK_LIB([m],[main],[],AC_MSG_ERROR([could not find math library]))
AC_CHECK_LIB([mpi],[main],[],AC_MSG_ERROR([could not find mpi library]))

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_HEADER_TIME

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([clock_gettime gettimeofday memset pow rint sqrt strtol strtoul])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
