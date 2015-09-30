dnl -*- Autoconf -*-

AC_DEFUN([ACX_GMP],
 [
  AC_ARG_WITH(gmp,
    AC_HELP_STRING([--with-gmp=DIR],[top directory for GnuMP]))

  acx_gmp_dir=
  gmp_search_dir=
  if test x"$with_gmp" = xyes || test x"$with_gmp" = x; then
   gmp_search_dir="$HOME /usr/local /opt"
  elif test x"$with_gmp" != xno; then
   acx_gmp_dir="$with_gmp"
  fi

  AC_MSG_CHECKING([for GnuMP library])
  AC_LANG_PUSH(C++)
  acx_gmp=no
  acx_save_CPPFLAGS="$CPPFLAGS"
  acx_save_LDFLAGS="$LDFLAGS"
  if test x"$acx_gmp_dir" = x; then
   AC_TRY_LINK([
#include "gmp.h"
], [], acx_gmp=yes; gmp_search_dir= )
  fi
  for try_dir in $gmp_search_dir; do
   CPPFLAGS="$acx_save_CPPFLAGS -I$try_dir/include"
   LDFLAGS="$acx_save_LDFLAGS -L$try_dir/lib"
   AC_TRY_LINK([
#include "gmp.h"
], [], 
      acx_gmp_dir="$try_dir" ; break)
  done
  CPPFLAGS="$acx_save_CPPFLAGS"
  LDFLAGS="$acx_save_LDFLAGS"
  AC_LANG_POP(C++)

  if test x"$acx_gmp_dir" != x; then
   acx_gmp=yes
   CPPFLAGS=" $CPPFLAGS -I$acx_gmp_dir/include"
   LDFLAGS="$LDFLAGS -L$acx_gmp_dir/lib"
   AC_MSG_RESULT([[$acx_gmp_dir]])
  elif test x"$acx_gmp" = xyes; then
   AC_MSG_RESULT([yes])
  else
   AC_MSG_RESULT([no])
  fi

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
  if test x"$acx_gmp" = xyes; then
   ifelse([$1],,AC_DEFINE(HAVE_GMP,1,[Define if you have the GnuMP library.]),[$1])
        :
  else
	echo here
        acx_gmp=no
        $2
  fi
])dnl ACX_GMP

dnl
dnl check for 'restrict' keyword, in C98 but a non-standard extension in C++
dnl
AC_DEFUN(
 [ACX_CXX_RESTRICT],
 [AC_CACHE_CHECK(
   [whether the C++ compiler supports the restrict keyword],
   [acx_cv_cxx_restrict],
   acx_cv_cxx_restrict=unsupported
    AC_LANG_PUSH(C++)
    for acx_kw in restrict __restrict__ __restrict ; do
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[double* $acx_kw x; ]])], 
	               acx_cv_cxx_restrict="$acx_kw" ; break)
    done
    AC_LANG_POP(C++)
   
  )
  if test "$acx_cv_cxx_restrict" != unsupported; then
   acx_kw="$acx_cv_cxx_restrict"
   AC_DEFINE(HAVE_CXX_RESTRICT,,[Defined if the compiler supports the restrict keyword.])
  else
   acx_kw=""
  fi
  if test x"$acx_kw" != xrestrict; then
   AC_DEFINE_UNQUOTED(restrict,$acx_kw,[Make restrict keyword work])
  fi
 ]
)


dnl
dnl check for __attribute__((noreturn))
dnl
AC_DEFUN(
 [ACX_CXX_NORETURN],
 [AC_CACHE_CHECK(
  [whether the C++ compiler supports __attribute__((noreturn))],
  [acx_cv_cxx_noreturn],
   acx_cv_cxx_noreturn=unsupported
    AC_LANG_PUSH(C++)
    for acx_kw in "__attribute__((noreturn))" ; do
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[void foo() $acx_kw; ]])], 
	               acx_cv_cxx_noreturn="$acx_kw" ; break)
    done
    AC_LANG_POP(C++)
   
  )
  if test "$acx_cv_cxx_noreturn" != unsupported; then
   acx_kw="$acx_cv_cxx_noreturn"
   AC_DEFINE(HAVE_CXX_NORETURN,,[Defined if the compiler supports __attribute__((noreturn)).])
  else
   acx_kw=""
  fi
  AC_DEFINE_UNQUOTED(FUNC_NORETURN,$acx_kw,[((noreturn)) attribute])
 ]
)
  
dnl
dnl Check for boost headers
dnl
AC_DEFUN([ACX_HAVE_BOOST],
 [
  AC_ARG_WITH([boost], 
     AC_HELP_STRING([--with-boost=DIR],[Boost top directory]))

dnl boost-dir overrides the default search directories

  if test x"$with_boost" = xyes || test x"$with_boost" = x; then
   boost_search_dir=". $HOME/include $HOME/boost $HOME $prefix $prefix/boost $prefix /usr/local /usr/local/include /usr /usr/src"
  elif test x"$with_boost" = xno; then
   boost_search_dir=""
  else
   boost_search_dir=`echo "$with_boost" | sed 's,//*,/,g;s,/$,,'`
  fi

dnl search for boost include files

  AC_MSG_CHECKING([for Boost include files])
  AC_LANG_PUSH(C++)
  acx_boost_incdir=
  BOOST_CPPFLAGS=
  acx_save_CPPFLAGS="$CPPFLAGS"
  for try_incdir in $boost_search_dir; do
   CPPFLAGS="$acx_save_CPPFLAGS -I$try_incdir"
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <boost/config.hpp>
]], [])],
      acx_boost_incdir=$try_incdir ; break)
  done
  CPPFLAGS="$acx_save_CPPFLAGS"
  AC_LANG_POP(C++)

  if test x"$acx_boost_incdir" != x; then
   BOOST_CPPFLAGS="-I$acx_boost_incdir"
   AC_MSG_RESULT([[$BOOST_CPPFLAGS]])
  else
   AC_MSG_RESULT([no])
  fi

  CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
 ]
) dnl ACX_HAVE_BOOST


dnl
dnl AX_BOOST_PROGRAM_OPTIONS
dnl
AC_DEFUN([AX_BOOST_PROGRAM_OPTIONS],
 [
  AC_REQUIRE([ACX_HAVE_BOOST])
  AC_ARG_WITH([boost-program-options],AS_HELP_STRING([--with-boost-program-options],
   [specify the boost program-options library or suffix to use]),
   [
    if test "x$with_boost_program_options" != "xno"; then
      ax_program_options_lib=$with_boost_program_options
    else
      ax_program_options_lib=
    fi
   ])

  AC_MSG_CHECKING([whether the boost::program_options library is available])
  AC_LANG_PUSH(C++)
  save_LIBS="$LIBS"
  BOOST_PROGRAM_OPTIONS_LIB=
  ax_plib="$ax_program_options_lib"
  for ax_lib in "$ax_plib" "-l$ax_plib" "-lboost_program_options-$ax_plib" "-l$ax_plib-$CC" "-lboost_program_options" "-lboost_program_options-$CC"; do
    LIBS="$LIBS $ax_lib"
    AC_LINK_IFELSE([AC_LANG_SOURCE(
     [[
      #include <boost/program_options.hpp>
      namespace po = boost::program_options;
      int main(int ac, char** av){
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("compression", po::value<int>(), "set compression level")
        ;

        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm); }
      ]]
     )], [BOOST_PROGRAM_OPTIONS_LIB=$ax_lib ; break])
     LIBS="$save_LIBS"
  done
  AC_LANG_POP(C++)

  if test x"$BOOST_PROGRAM_OPTIONS_LIB" != x; then
   AC_MSG_RESULT([[$BOOST_PROGRAM_OPTIONS_LIB]])
  else
   AC_MSG_RESULT([no])
  fi
  AC_SUBST(BOOST_PROGRAM_OPTIONS_LIB)
 ]
) dnl AX_BOOST_PROGRAM_OPTIONS


dnl
dnl check convention for returning complex parameters from fortran functions
dnl
AC_DEFUN(
 [ACX_FORTRAN_COMPLEX_RETURN],
 [AC_REQUIRE([ACX_BLAS])
  AC_CACHE_CHECK(
  [convention for returning complex values from Fortran functions],
  [acx_cv_fortran_complex_return],
   acx_cv_fortran_complex_return=unsupported
    AC_LANG_PUSH(C++)
    CXXFLAGS_save="$CXXFLAGS"
    LIBS_save="$LIBS"
    LIBS="$BLAS_LIBS $FLIBS"
    CXXFLAGS="-I. -I$srcdir $CXXFLAGS"
    AC_RUN_IFELSE([AC_LANG_SOURCE([[
#define FORTRAN_COMPLEX_RETURN_FIRST_ARG
#include "common/blas1f.h"
#include <complex>
#include <vector>
typedef std::complex<double> complex;
int main()
{
   std::vector<complex> v1(10, complex(2.0, 2.0));
   std::vector<complex> v2(10, 1.0);
   complex res = BLAS::zdotu(10, &(*v1.begin()), 1, &(*v2.begin()), 1);
   return std::abs(res - complex(20.0, 20.0)) < 1E-10 ? 0 : 1;
}
    ]])], acx_cv_fortran_complex_return=pass_as_first_argument,
	[AC_RUN_IFELSE([AC_LANG_SOURCE([[
#define FORTRAN_COMPLEX_RETURN_IN_REGISTER
#include "common/blas1f.h"
#include <complex>
#include <vector>
typedef std::complex<double> complex;
int main()
{
   std::vector<complex> v1(10, complex(2.0, 2.0));
   std::vector<complex> v2(10, 1.0);
   complex res = BLAS::zdotu(10, &(*v1.begin()), 1, &(*v2.begin()), 1);
   return std::abs(res - complex(20.0, 20.0)) < 1E-10 ? 0 : 1;
}
    ]])], acx_cv_fortran_complex_return=return_in_register)])
    LIBS="$LIBS_save"
    CXXFLAGS="$CXXFLAGS_save"
    AC_LANG_POP(C++)
   )
   AS_IF([test "$acx_cv_fortran_complex_return" == pass_as_first_argument], [
     AC_DEFINE(FORTRAN_COMPLEX_RETURN_FIRST_ARG,,[Defined if the Fortran returns complex as first arg])
   ], [AS_IF([test "$acx_cv_fortran_complex_return" == return_in_register], [
     AC_DEFINE(FORTRAN_COMPLEX_RETURN_IN_REGISTER,,[Defined if the Fortran returns complex in registers])
      ]
   )])
 ]
) dnl ACX_FORTRAN_COMPLEX_RETURN


dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/acx_blas.html
dnl
AC_DEFUN([ACX_BLAS], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
acx_blas_ok=no

AC_ARG_WITH(blas,
	[AC_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
	yes | "") ;;
	no) acx_blas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
	*) BLAS_LIBS="-l$with_blas" ;;
esac

# Get fortran linker names of BLAS functions to check for.
AC_F77_FUNC(sgemm)
AC_F77_FUNC(dgemm)

acx_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test $acx_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
	AC_MSG_CHECKING([for $sgemm in $BLAS_LIBS])
	AC_TRY_LINK_FUNC($sgemm, [acx_blas_ok=yes], [BLAS_LIBS=""])
	AC_MSG_RESULT($acx_blas_ok)
	LIBS="$save_LIBS"
fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $acx_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS"
	AC_CHECK_FUNC($sgemm, [acx_blas_ok=yes])
	LIBS="$save_LIBS"
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_CHECK_LIB(f77blas, $sgemm,
		[AC_CHECK_LIB(cblas, cblas_dgemm,
			[acx_blas_ok=yes
			 BLAS_LIBS="-lcblas -lf77blas -latlas"],
			[], [-lf77blas -latlas])],
			[], [-latlas])])
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(dgemm, $dgemm,
		[AC_CHECK_LIB(sgemm, $sgemm,
			[acx_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
			[], [-lblas])],
			[], [-lblas])])
fi

# BLAS in Alpha CXML library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(cxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(dxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-ldxml"])
fi

# BLAS in Sun Performance library?
if test $acx_blas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, $sgemm,
        			[BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 acx_blas_ok=yes],[],[-lsunmath])])
	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(scs, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(complib.sgimath, $sgemm,
		     [acx_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(essl, $sgemm,
			[acx_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
			[], [-lblas $FLIBS])])
fi

# BLAS in AMD ACLM library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(acml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-lacml"])
fi

# Generic BLAS library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lblas"])
fi

AC_SUBST(BLAS_LIBS)

LIBS="$acx_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_blas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])
        :
else
        acx_blas_ok=no
        $2
fi
])dnl ACX_BLAS


dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/acx_lapack.html
dnl
AC_DEFUN([ACX_LAPACK], [
AC_REQUIRE([ACX_BLAS])
acx_lapack_ok=no

AC_ARG_WITH(lapack,
        [AC_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>])])
case $with_lapack in
        yes | "") ;;
        no) acx_lapack_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) LAPACK_LIBS="$with_lapack" ;;
        *) LAPACK_LIBS="-l$with_lapack" ;;
esac

# Get fortran linker name of LAPACK function to check for.
AC_F77_FUNC(cheev)

# We cannot use LAPACK if BLAS is not found
if test "x$acx_blas_ok" != xyes; then
        acx_lapack_ok=noblas
fi

# First, check LAPACK_LIBS environment variable
if test "x$LAPACK_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
        AC_MSG_CHECKING([for $cheev in $LAPACK_LIBS])
        AC_TRY_LINK_FUNC($cheev, [acx_lapack_ok=yes], [LAPACK_LIBS=""])
        AC_MSG_RESULT($acx_lapack_ok)
        LIBS="$save_LIBS"
        if test acx_lapack_ok = no; then
                LAPACK_LIBS=""
        fi
fi

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $acx_lapack_ok = no; then
        save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
        AC_CHECK_FUNC($cheev, [acx_lapack_ok=yes])
        LIBS="$save_LIBS"
fi

# Generic LAPACK library?
for lapack in lapack lapack_rs6k; do
        if test $acx_lapack_ok = no; then
                save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
                AC_CHECK_LIB($lapack, $cheev,
                    [acx_lapack_ok=yes; LAPACK_LIBS="-l$lapack"], [], [$FLIBS])
                LIBS="$save_LIBS"
        fi
done

AC_SUBST(LAPACK_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_lapack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_LAPACK,1,[Define if you have LAPACK library.]),[$1])
        :
else
        acx_lapack_ok=no
        $2
fi
])dnl ACX_LAPACK

dnl check for ARPACK
AC_DEFUN([ACX_ARPACK], 
[
        AC_REQUIRE([ACX_LAPACK])

        AC_ARG_WITH(arpack,
                [AS_HELP_STRING([--with-arpack=<lib>], 
                        [use ARPACK library, default=-larpack])],
                [with_arpack=$withval],
                [with_arpack=check])

        LIBARPACK=
        if test "x$with_arpack" = "xno"; then
                acx_want_arpack="no"
        elif test "x$with_arpack" = "xyes"; then
                acx_want_arpack="yes"
                LIBARPACK="-larpack"
        elif test "x$with_arpack" = "xcheck"; then
                acx_want_arpack="check"
                LIBARPACK="-larpack"
        else
                acx_want_arpack="yes"
                LIBARPACK="$with_arpack"
        fi
        
        if test "x$acx_want_arpack" != "xno"; then
                AC_MSG_CHECKING([for znaupd in $LIBARPACK])
                dnl get fortran name of the function we are interested in
                AC_F77_FUNC(znaupd)
                save_LIBS="$LIBS"; LIBS="$LIBARPACK $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
                acx_arpack_ok=no
                AC_TRY_LINK_FUNC([$znaupd], [acx_arpack_ok=yes], [LIBARPACK=""])
                LIBS="$save_LIBS"
                AC_MSG_RESULT($acx_arpack_ok)
                if test "x$acx_want_arpack" = "xyes" -a "x$acx_arpack_ok" != "xyes"; then
                        AC_MSG_FAILURE([--with-arpack was given, but the test for ARPACK failed])
                fi
                if test "x$acx_arpack_ok" = "xyes"; then
                        AC_DEFINE([HAVE_LIBARPACK], [1], [Define if you have the ARPACK library.])
                fi
        fi
        AC_SUBST(LIBARPACK)

]) dnl ACX_ARPACK

