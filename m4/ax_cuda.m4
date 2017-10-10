# SYNOPSIS
#
#   AX_CUDA([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#


AC_DEFUN([AX_CUDA], [
ax_cuda_ok=no

AC_ARG_WITH(cuda,
        [AS_HELP_STRING([--with-cuda=<lib>], [use CUDA library <lib>])])

case $with_cuda in
        yes) ;;
        no | "") ax_cuda_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) CUDA_LIBS="$with_cuda" ;;
        *) CUDA_LIBS="-l$with_cuda" ;;
esac

if test "x$ax_cuda_ok" = xdisable; then
   CUDA_LIBS=""
elif test "x$CUDA_LIBS" != x; then
# First, check CUDA_LIBS environment variable
        save_LIBS="$LIBS"; LIBS="$CUDA_LIBS $LIBS"
        AC_MSG_CHECKING([for cudaGetDeviceCount in $CUDA_LIBS])
        AC_TRY_LINK_FUNC(cudaGetDeviceCount, [ax_cuda_ok=yes], [CUDA_LIBS=""])
        AC_MSG_RESULT($ax_cuda_ok)
        LIBS="$save_LIBS"
        if test $ax_cuda_ok = no; then
                CUDA_LIBS=""
        fi
else
   # CUDA library in standard place?
   for cuda in cudart ; do
           if test $ax_cuda_ok = no; then
                   AC_CHECK_LIB($cuda, cudaGetDeviceCount,
                      [ax_cuda_ok=yes; CUDA_LIBS="-l$cuda"], [], [$LIBS])
           fi
   done
fi

AC_SUBST(CUDA_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:

if test x"$ax_cuda_ok" = xyes; then
        AC_DEFINE(HAVE_CUDA,[1],[Define if you have CUDA library.])
        AC_SUBST([HAVE_CUDA], [1])
        $1
else
        $2
        echo -n
fi
])dnl AX_CUDA
