dnl -$- Autoconf -$-
dnl Process this file with autoconf to produce a configure script.

AC_INIT([Matrix Product Toolkit], [(git)], [ianmcc@physics.uq.edu.au], 
        [mptoolkit], [https://people.smp.uq.edu.au/IanMcCulloch/mptoolkit/])

dnl use C++ mode by default
AC_LANG(C++)

AC_CONFIG_SRCDIR([common/trace.h])
AC_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([m4])
m4_include([m4/macros.m4])
AC_PREFIX_DEFAULT("$HOME")

CXXFLAGS_MPTK_SAVE="$CXXFLAGS"

dnl debug flags -- needs to be called before any AC_PROG_* macros
AX_CHECK_ENABLE_DEBUG(no, [CONFIG_DEBUG], [CONFIG_NODEBUG], no)

dnl Checks for compilers
AC_PROG_CXX
dnl AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS

dnl put DNDEBUG on the command line, rather than in config.h
if test "x$enable_debug" != "xyes" ; then
   CFLAGS="${CFLAGS} -DNDEBUG"
   CXXFLAGS="${CXXFLAGS} -DNDEBUG"

   dnl backup cxx flags, if nothing has been specified.
   if test "x$CXXFLAGS_MPTK_SAVE" = "x" ; then
      AC_MSG_NOTICE([Setting default CXXFLAGS])
      CXXFLAGS="${CXXFLAGS} -O2"
      AC_LANG_PUSH(C++)
      AX_CHECK_COMPILE_FLAG("-march=native", [CXXFLAGS="${CXXFLAGS} -march=native"])
      AX_CHECK_COMPILE_FLAG("-flto", [CXXFLAGS="${CXXFLAGS} -flto"])
      AX_CHECK_COMPILE_FLAG("-Wno-lto-type-mismatch", [CXXFLAGS="${CXXFLAGS} -Wno-lto-type-mismatch"])
      AC_LANG_POP(C++)
      AC_LANG_PUSH([Fortran 77])
      AX_CHECK_COMPILE_FLAG("-flto", [FFLAGS="${FFLAGS} -flto"])
      AC_LANG_POP([Fortran 77])
      AC_MSG_NOTICE([CXXFLAGS are $CXXFLAGS])
   fi
fi

dnl C++14 support
AX_CXX_COMPILE_STDCXX_14(noext,mandatory)

dnl we have C++11 as mandatory, but set it in config.h too
AC_DEFINE(HAVE_CXX11,1, [define if the compiler supports basic C++11 syntax])
AC_SUBST(HAVE_CXX11)

dnl we have C++14 as mandatory, but set it in config.h too
AC_DEFINE(HAVE_CXX14,1, [define if the compiler supports basic C++14 syntax])
AC_SUBST(HAVE_CXX14)

AX_CXX_VAR_PRETTYFUNC

dnl C++ compiler version and vendor information
AC_LANG_PUSH(C++)
AX_COMPILER_VENDOR
AC_DEFINE_UNQUOTED([CONFIG_COMPILER_VENDOR], ["$ax_cv_cxx_compiler_vendor"], [C++ compiler vendor])
AX_COMPILER_VERSION
AC_DEFINE_UNQUOTED([CONFIG_COMPILER_VERSION], ["$ax_cv_cxx_compiler_version"], [C++ compiler version])
AC_LANG_POP(C++)

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

dnl SLURM integration
AC_ARG_WITH([slurm], AC_HELP_STRING([--with-slurm@<:@=ARG@:>@],[SLURM queue integration for checkpointing, libslurm from standard location (ARG=yes), from the specified location (ARG=<path>) or disable it (ARG=no) @<:@ARG=yes@:>@]),
           [], [with_slurm=check])

if test "x$with_slurm" != xno; then
   AC_MSG_CHECKING([for SLURM library])
   slurm_lib="$with_slurm"
   if test "x$slurm_lib" == xcheck; then
      slurm_lib="-lslurm";
   fi
   save_LIBS="$LIBS"; LIBS="$slurm_lib $LIBS"
   AC_MSG_CHECKING([for slurm_get_end_time in "$slurm_lib"])
   AC_TRY_LINK_FUNC([slurm_get_end_time], 
      [acx_slurm_ok=yes
       AC_SUBST([SLURM_LIBS], "$slurm_lib")
       AC_DEFINE([CONFIG_SLURM], [], [Define if you want SLURM queue integration])
      ], 
      [
       if test "x$with_slurm" != xcheck; then
         AC_MSG_FAILURE([--with-slurm was given, but test for libslurm failed])
       fi
       LIBSLURM=""
       acx_slurm_ok=no
      ])
   AC_MSG_RESULT($acx_slurm_ok)
   LIBS="$save_LIBS"
fi


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

AX_BOOST_BASE([1.47], [], [AC_MSG_ERROR([boost libraries not found, aborting.])])
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

dnl Checks for library functions.
AC_FUNC_MEMCMP
AC_TYPE_SIGNAL
AC_CHECK_FUNCS(gettimeofday)
AC_HEADER_TIOCGWINSZ
AC_CHECK_HEADERS(sys/ioctl.h)
AC_FUNC_STRERROR_R
AC_C_BIGENDIAN

dnl CUDA runtime support

dnl if we're using cublas then force detection of cuda by default
if test "x$with_cuda" = x; then
   if test "x$with_cublas" != xno && test "x$with_cublas" != x; then
      with_cuda=yes
   fi
fi

if test "x$with_cuda" = xyes; then
   AX_CUDA_CONFIG

   dnl cublas support
   AX_CUBLAS([AC_CHECK_HEADERS(cublas_v2.h, [], [AC_MSG_ERROR([Cannot find cublas_v2.h header])])], 
   [ AS_IF([test "x$ax_cublas_ok" != xdisable], [AC_MSG_ERROR([Cannot find cuBLAS library])]) ])

   dnl cusovler support

   AX_CUSOLVER([AC_CHECK_HEADERS(cusolverDn.h, [], [AC_MSG_ERROR([Cannot find cusolverDn.h header])])], 
   [ AS_IF([test "x$ax_cusolver_ok" != xdisable], [AC_MSG_ERROR([Cannot find cuSOLVER library])]) ])
fi

dnl valgrind support
AC_CHECK_HEADERS(valgrind/memcheck.h)

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

dnl save the compiler flags to the config header
esccxxflags="${CXXFLAGS//\"/\\\"}"
AC_DEFINE_UNQUOTED([CONFIG_CXXFLAGS], ["${CXXFLAGS//\"/\\\"}"], [C++ compiler flags])

AC_OUTPUT(Makefile, echo timestamp > stamp-h)
