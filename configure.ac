dnl -$- Autoconf -$-
dnl Process this file with autoconf to produce a configure script.

AC_INIT([Matrix Product Toolkit], [HEAD], [ianmcc@physics.uq.edu.au], 
        [mptoolkit], [http://physics.uq.edu.au/people/ianmcc/mptoolkit])

AC_CONFIG_SRCDIR([common/trace.h])
AC_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([m4])
m4_include([m4/macros.m4])
AC_PREFIX_DEFAULT("$HOME")

dnl debug flags -- needs to be called before any AC_PROG_* macros

AX_CHECK_ENABLE_DEBUG(no, [CONFIG_DEBUG], [CONFIG_NODEBUG], no)

dnl put DNDEBUG on the command line, rather than in config.h
if test "x$enable_debug" != "xyes" ; then
   CFLAGS="${CFLAGS} -DNDEBUG"
   CXXFLAGS="${CXXFLAGS} -DNDEBUG"
fi

dnl Checks for compilers
AC_PROG_CXX
dnl AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS

dnl backup cxx flags, if nothing has been specified

if test "x$CXXFLAGS" = "x" ; then
   AC_MSG_NOTICE([Setting default CXXFLAGS])
   CXXFLAGS="-O2"
   AC_LANG_PUSH(C++)
   AX_CHECK_COMPILE_FLAG("-march=native", [CXXFLAGS="${CXXFLAGS} -march=native"])
   AC_LANG_POP(C++)
fi

dnl C++11 support
AX_CXX_COMPILE_STDCXX_11(noext,optional)

dnl large file support
AC_SYS_LARGEFILE

dnl openmp support, defines variable OPENMP_CXXFLAGS and sets HAVE_OPENMP
AC_ARG_WITH([openmp], [AC_HELP_STRING([--with-openmp], [enable OpenMP parallelization @<:@default=check@:>@])],
	[],
	[with_open=check])
if test "x$with_openmp" != xno; then
   AC_LANG_PUSH(C++)
   AX_OPENMP([AC_SUBST([OPENMP_CXXFLAGS]) AC_SUBST([HAVE_OPENMP], [1])], [AC_SUBST([HAVE_OPENMP], []) 
                                        if test "x$with_openmp" = yes; then
                     AC_MSG_FAILURE([--with-openmp was given, but the test for OpenMP failed!])
                  fi
                 ])
   AC_LANG_POP(C++)
fi

dnl MPI support, if found sets HAVE_MPI and defines MPICXX and MPILIBS
AC_ARG_WITH([mpi], [AC_HELP_STRING([--with-mpi], [enable MPI @<:@default=check@:>@])],
	[],
	[with_mpi=check])
if test "x$with_mpi" != xno; then
   AC_LANG_PUSH(C++)
   AX_MPI([AC_SUBST([HAVE_MPI], [1])], [AC_SUBST([HAVE_MPI], []) if test "x$with_mpi" = yes; then
                  AC_MSG_FAILURE([--with-mpi was given, but the test for MPI failed!])
               fi
              ])
   AC_LANG_POP(C++)
fi

AC_ARG_WITH([pstream], AC_HELP_STRING([--without-pstream],[disables persistent streams support]))
test x"$with_pstream" != xno && AC_DEFINE([USE_PSTREAM], [], [Define if you want persistent streams])

dnl AC_CHECK_LIB(pthread, pthread_create)

dnl Checks for header files.
AC_CHECK_HEADERS(fcntl.h unistd.h)

dnl Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_OFF_T
AC_TYPE_SIZE_T

AC_CHECK_TYPE([long long], [])

   AC_DEFINE([HAVE_LONGLONG], [], [Do we have long long?])
   AC_ARG_WITH([longlong], AC_HELP_STRING([--without-longlong],
               [disables 'long long' support in the streaming lib]))
   test x"$with_longlong" != xno && AC_DEFINE([USE_LONGLONG], [], [If we need to use long long])
   AC_CHECK_SIZEOF(long long)


AC_DEFINE([CURRENTFORMAT_AUTOCONF], [], 
   [Define if we can determine the pstream binary format from autoconf])

AC_DEFINE([CONFIG_EXPOKIT], [],
   [Define if we want to include functions from EXPOKIT])

AC_DEFINE([DEFAULT_PAGE_SIZE], [65536],
	[Default page size for file I/O])

AC_DEFINE([DEFAULT_PAGE_CACHE_SIZE], [67108864],
	[Default page cache size])

AC_CHECK_SIZEOF(short)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(float)
AC_CHECK_SIZEOF(double)

dnl Checks for GnuMP
ACX_GMP([], [AC_MSG_ERROR(["GnuMP library not found; see http://gmplib.org, or set --with-gmp=DIR"])])
AC_CHECK_LIB(gmp, __gmpz_init, [], [AC_MSG_ERROR([Problem with GnuMP install])])

dnl BLAS/LAPACK
AX_BLAS([], [AC_MSG_ERROR([BLAS library not found, aborting.])])
AX_LAPACK([], [AC_MSG_ERROR([LAPACK library not found, aborting.])])

AX_BOOST_BASE([1.35], [], [AC_MSG_ERROR([boost libraries not found, aborting.])])
AX_BOOST_PROGRAM_OPTIONS

ACX_ARPACK

if test "x$acx_arpack_ok" != "xyes"; then
  AC_MSG_ERROR([ARPACK is required, but not configured or not found, aborting.])
fi

if test x"$BOOST_PROGRAM_OPTIONS_LIB" == x; then
  AC_MSG_ERROR([boost program options library not found, aborting.])
fi

dnl C++ features
ACX_FORTRAN_COMPLEX_RETURN
ACX_CXX_RESTRICT
ACX_CXX_NORETURN

dnl Checks for library functions.
AC_FUNC_MEMCMP
AC_TYPE_SIGNAL
AC_CHECK_FUNCS(gettimeofday)
AC_HEADER_TIOCGWINSZ
AC_CHECK_HEADERS(sys/ioctl.h)
AC_FUNC_STRERROR_R
AC_C_BIGENDIAN

dnl Checks for programs

AC_PROG_INSTALL

dnl shared library configuration
# Check for enable shared library build option
private_enable_shared=yes
AC_ARG_ENABLE(shared, AC_HELP_STRING([--disable-shared],[build static libraries (default is shared libraries)]),
 [
  if test "x$enableval" = "xno" 
  then
   private_enable_shared=no
  fi
 ]
)

AC_OUTPUT(Makefile, echo timestamp > stamp-h)
