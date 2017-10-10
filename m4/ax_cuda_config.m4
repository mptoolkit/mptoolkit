# -*- Autoconf -*-


## ------------------------- ##
## Autoconf macros for CUDA. ##
## ------------------------- ##


# ----------------------------------------------------------------------
# CIT_CUDA_CONFIG
# ---------------
-------------------------------------------------------
# Determine the directory containing <cuda_runtime.h>
AC_DEFUN([AX_CUDA_CONFIG], [
ax_cuda_ok=no

AC_ARG_WITH(cuda,
        [AS_HELP_STRING([--with-cuda], [use CUDA])])

  # influential environment variables
  AC_ARG_VAR(NVCC, [NVIDIA nvcc compiler command])
  AC_ARG_VAR(NVFLAGS, [nvcc compiler flags])
  AC_ARG_VAR(CUDA_INC, [Location of CUDA include files])
  AC_ARG_VAR(CUDA_LIB, [Location of CUDA library libcudart])

  # tests NVCC variable
  AS_IF([test x"$NVCC" = x],[
    NVCC=nvcc
  ])

  # Check for compiler
  # checks if program in path
  AC_PATH_PROG(NVCC_PROG, $NVCC)
  if test -z "$NVCC_PROG" ; then
    AC_MSG_ERROR([cannot find '$NVCC' program, please check your PATH.])
  fi

  # Checks for compiling and linking
  AC_LANG_PUSH([C])
  AC_REQUIRE_CPP
  CFLAGS_save="$CFLAGS"
  LDFLAGS_save="$LDFLAGS"
  LIBS_save="$LIBS"

  # uses nvcc compiler
  CFLAGS="$CUDA_FLAGS"
  # NVFLAGS="$NVFLAGS $CXXFLAGS"
  if test "x$CUDA_INC" != "x"; then
    NVFLAGS="$NVFLAGS -I$CUDA_INC"
    CXXFLAGS="$CXXFLAGS -I$CUDA_INC"
  fi

  # Check for CUDA headers
  # runs test with nvcc
  AC_MSG_CHECKING([for cuda_runtime.h])
  ac_compile='$NVCC -c $NVFLAGS conftest.$ac_ext >&5'
  AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([[
#include <cuda.h>
#include <cuda_runtime.h>]],[[void* ptr = 0;
]])
  ], [
    AC_MSG_RESULT(yes)
  ], [
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([CUDA runtime header not found; try setting CUDA_INC.])
  ])

  # Check fo CUDA library
  if test "x$CUDA_LIB" != "x"; then
    CUDA_LDFLAGS="-L$CUDA_LIB"
    LDFLAGS="$CUDA_LDFLAGS $LDFLAGS"
  fi
  CUDA_LIBS="-lcudart"

  # runs compilation test with nvcc
  AC_MSG_CHECKING([nvcc compilation with cudaMalloc in -lcudart])
  ac_compile='$NVCC -c $NVFLAGS conftest.$ac_ext >&5'
  AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([[
#include <cuda.h>
#include <cuda_runtime.h>]],[[void* ptr = 0;cudaMalloc(&ptr, 1);
]])
  ], [
    AC_MSG_RESULT(yes)
  ], [
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([CUDA library function with nvcc compilation failed; try setting CUDA_INC.])
  ])

  # runs linking test with nvcc
  AC_MSG_CHECKING([nvcc linking with cudaMalloc in -lcudart])
  ac_link='$NVCC -o conftest$ac_exeext $NVFLAGS $CUDA_LDFLAGS conftest.$ac_ext $LIBS >&5'
  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>]],[[void* ptr = 0;cudaMalloc(&ptr, 1);
]])],
    [AC_MSG_RESULT(yes)],
    [AC_MSG_RESULT(no)
     AC_MSG_ERROR([CUDA library linking with nvcc failed; try setting CUDA_LIB.])
  ])

  # runs linking test with standard compiler
  AC_MSG_CHECKING([Trying to link a cuda program])

   # CUDA library in standard place?
   for cuda in cudart ; do
           if test $ax_cuda_ok = no; then
                   AC_CHECK_LIB($cuda, cudaGetDeviceCount,
                      [ax_cuda_ok=yes; CUDA_LIBS="-l$cuda"], [], [$LIBS])
           fi
   done

if test x"$ax_cuda_ok" = xyes; then
  AC_DEFINE(HAVE_CUDA,[1],[Define if you have CUDA library.])
  AC_SUBST([HAVE_CUDA], [1])
  AC_SUBST([NVCC])
  AC_SUBST([NVFLAGS])
  AC_SUBST([CUDA_LDFLAGS])
  AC_SUBST([CUDA_LIBS])
  $1
else
  
  $2
  echo -n
fi

  CFLAGS="$CFLAGS_save"
  LDFLAGS="$LDFLAGS_save"
  LIBS="$LIBS_save"
  AC_LANG_POP([C])

])dnl CUDA_CONFIG


dnl end of file
