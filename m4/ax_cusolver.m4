# SYNOPSIS
#
#   AX_CUSOLVER([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])

AC_DEFUN([AX_CUSOLVER], [
AC_REQUIRE([AX_CUDA])
AC_REQUIRE([AX_OPENMP])
ax_cusolver_ok=no

AC_ARG_WITH(cusolver,
        [AS_HELP_STRING([--with-cusolver=<lib>], [use CUSOLVER library <lib> (requires CUDA)])])

case $with_cusolver in
        yes | "") ;;
        no) ax_cusolver_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) CUSOLVER_LIBS="$with_cusolver" ;;
        *) CUSOLVER_LIBS="-l$with_cusolver" ;;
esac

# We cannot use CUSOLVER if CUDA is not found
if test "x$ax_cuda_ok" = xdisable; then
        CUSOLVER_LIBS=""
elif test "x$CUSOLVER_LIBS" != x; then
# check CUSOLVER_LIBS environment variable
        save_LIBS="$LIBS"; LIBS="$CUSOLVER_LIBS $CUDA_LIBS $LIBS"
        AC_MSG_CHECKING([for cusolverDnCreate in $CUSOLVER_LIBS])
        save_CXXFLAGS="$CXXFLAGS"; CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
        AC_TRY_LINK_FUNC(cusolverDnCreate, [ax_cusolver_ok=yes], [CUSOLVER_LIBS=""])
        AC_MSG_RESULT($ax_cusolver_ok)
        LIBS="$save_LIBS"
        CXXFLAGS="$save_CXXFLAGS"
        if test $ax_cusolver_ok = no; then
                CUSOLVER_LIBS=""
        fi
else
   # CUSOLVER in standard location?
   for cusolver in cusolver; do
           if test $ax_cusolver_ok = no; then
                   save_LIBS="$LIBS"; LIBS="$CUDA_LIBS $LIBS"
                   save_CXXFLAGS="$CXXFLAGS"; CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
                   AC_CHECK_LIB($cusolver, cusolverDnCreate,
                      [ax_cusolver_ok=yes; CUSOLVER_LIBS="-l$cusolver"], [], [$FLIBS])
                   LIBS="$save_LIBS"
                   CXXFLAGS="$save_CXXFLAGS"
           fi
   done
fi

AC_SUBST(CUSOLVER_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_cusolver_ok" = xyes; then
        AC_DEFINE(HAVE_CUSOLVER,[1],[Define if you have CUSOLVER library.])
	AC_SUBST([HAVE_CUSOLVER], [1])
        $1
else
        $2
fi
])dnl AX_CUSOLVER
