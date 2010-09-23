#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([warg], [0.0.1], [x.didelot@warwick.ac.uk])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR([src/movescaletree.h])
AC_CONFIG_HEADERS([src/config.h])
AM_INIT_AUTOMAKE

# Check whether the user gave CXXFLAGS
if test x"${CXXFLAGS+set}" = xset; then
  # the user set CXXFLAGS; don't override it.
  flags_were_set=true
else
  flags_were_set=false
fi 

# Check whether the user wants to use a particular MPI installation
AC_MSG_CHECKING(for MPI)
AC_ARG_WITH(mpi,
  [  --with-mpi=dir       give the directory prefix for the MPI implementation.]
  [                        []],
  MPIdir="$withval/bin", MPIdir="")
AC_MSG_RESULT($MPIdir)

MPICC="mpicc mpcc"
MPICXX="mpicxx mpic++ mpiCC mpCC"
if ( test ! "x$MPIdir" = "x") then
        AC_PATH_PROGS(CC, $MPICC, mpicc, $MPIdir)
        AC_PATH_PROGS(CXX, $MPICXX, mpicxx, $MPIdir)
fi


# Checks for programs.
# the following two lines look for MPI compilers
# AC_PROG_CXX($MPICXX)
# AC_PROG_CC($MPICC)
# the next two lines look for standard compilers
AC_PROG_CXX
AC_PROG_CC

# warn the user if MPI wasn't found
if( test "x$MPIdir" = "x" && test ! "x$CXX" = "xmpicxx" && test ! "x$CXX" = "xmpic++" && test ! "x$CXX" = "xmpiCC" && test ! "x$CXX" = "xmpCC" ) then
  AC_WARN( "Please ensure that CXX is set to an MPI c++ compiler.  If the MPI compiler is not in your PATH, its path may be given by the --with-mpi= configuration option.  Alternatively, ensure that all necessary MPI compiler flags have been added to CXXFLAGS" );
fi


# Set CXXFLAGS empty unless user specified otherwise
if $flags_were_set; then
  ABC=
else
  CXXFLAGS=
fi 

# Checks for libraries.
# find gsl, wherever it's hiding on the system
#AC_CHECK_LIB([m],[cos])
#AC_CHECK_LIB([gslcblas],[cblas_dgemm])
#AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
AX_PATH_GSL

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([floor pow])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
