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
#include <cstddef> // workaround for bug https://gcc.gnu.org/gcc-4.9/porting_to.html
#include "gmp.h"
], [], acx_gmp=yes; gmp_search_dir= )
  fi
  for try_dir in $gmp_search_dir; do
   CPPFLAGS="$acx_save_CPPFLAGS -I$try_dir/include"
   LDFLAGS="$acx_save_LDFLAGS -L$try_dir/lib"
   AC_TRY_LINK([
#include <cstddef> // workaround for bug https://gcc.gnu.org/gcc-4.9/porting_to.html
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
dnl check convention for returning complex parameters from fortran functions
dnl
AC_DEFUN(
 [ACX_FORTRAN_COMPLEX_RETURN],
 [AC_REQUIRE([AX_BLAS])
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


dnl check for ARPACK
AC_DEFUN([ACX_ARPACK], 
[
        AC_REQUIRE([AX_LAPACK])

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

