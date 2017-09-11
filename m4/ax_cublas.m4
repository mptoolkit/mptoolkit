# SYNOPSIS
#
#   AX_CUBLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])

AC_DEFUN([AX_CUBLAS], [
AC_REQUIRE([AX_CUDA])
ax_cublas_ok=no

AC_ARG_WITH(cublas,
        [AS_HELP_STRING([--with-cublas=<lib>], [use CUBLAS library <lib> (requires CUDA)])])

case $with_cublas in
        yes | "") ;;
        no) ax_cublas_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) CUBLAS_LIBS="$with_cublas" ;;
        *) CUBLAS_LIBS="-l$with_cublas" ;;
esac

# We cannot use CUBLAS if CUDA is not found
if test "x$ax_cuda_ok" = xdisable; then
        CUBLAS_LIBS=""
elif test "x$CUBLAS_LIBS" != x; then
# check CUBLAS_LIBS environment variable
        save_LIBS="$LIBS"; LIBS="$CUBLAS_LIBS $CUDA_LIBS $LIBS"
        AC_MSG_CHECKING([for cublasSetMatrix in $CUBLAS_LIBS])
        AC_TRY_LINK_FUNC(cublasSetMatrix, [ax_cublas_ok=yes], [CUBLAS_LIBS=""])
        AC_MSG_RESULT($ax_cublas_ok)
        LIBS="$save_LIBS"
        if test $ax_cublas_ok = no; then
                CUBLAS_LIBS=""
        fi
else
   # CUBLAS in standard location?
   for cublas in cublas; do
           if test $ax_cublas_ok = no; then
                   save_LIBS="$LIBS"; LIBS="$CUDA_LIBS $LIBS"
                   AC_CHECK_LIB($cublas, cublasSetMatrix,
                      [ax_cublas_ok=yes; CUBLAS_LIBS="-l$cublas"], [], [$FLIBS])
                   LIBS="$save_LIBS"
           fi
   done
fi

AC_SUBST(CUBLAS_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_cublas_ok" = xyes; then
        AC_DEFINE(HAVE_CUBLAS,[1],[Define if you have CUBLAS library.])
	AC_SUBST([HAVE_CUBLAS], [1])
        $1
else
        $2
fi
])dnl AX_CUBLAS
