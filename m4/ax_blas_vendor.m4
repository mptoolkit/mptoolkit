# SYNOPSIS
#
#   AX_BLAS_VENDOR
#
# DESCRIPTION
#
#   Attempt to determine the vendor of the BLAS library

AC_DEFUN([AX_BLAS_VENDOR], [
AC_PREREQ([2.55])
AC_REQUIRE([AX_BLAS])

AC_MSG_CHECKING([BLAS vendor information])
save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
BLAS_VENDOR="generic"

# MKL has the MKL_Get_Version() function
AC_LINK_IFELSE([AC_LANG_CALL([], [MKL_Get_Version])], [BLAS_VENDOR="MKL"], [])

# OpenBLAS doesn't have a function to get the version, just a config string
AC_LINK_IFELSE([AC_LANG_CALL([], [openblas_get_config])], [BLAS_VENDOR="OpenBLAS"], [])

AC_MSG_RESULT([$BLAS_VENDOR])

AC_DEFINE_UNQUOTED([BLAS_VENDOR], [$BLAS_VENDOR], [The vendor of the BLAS libary])
AC_SUBST([BLAS_VENDOR])

])dnl AX_BLAS_VENDOR
