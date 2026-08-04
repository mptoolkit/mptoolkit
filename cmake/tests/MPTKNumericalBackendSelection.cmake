cmake_minimum_required(VERSION 3.20)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/..")
include(MPTKNumericalBackend)

# Keep this regression independent of the numerical libraries installed on the
# test host.  ARPACK is required by MPToolkit, so model its observed GNU
# Fortran dependency and make the MKL discovery results deterministic.
function(mptk_detect_shared_library_dependencies library output_var)
  set(${output_var} "libgfortran.so.5" PARENT_SCOPE)
endfunction()

function(find_library output_var)
  if(output_var STREQUAL "_MPTK_MKL_GF_LP64_LIBRARY")
    set(${output_var} "/mock/libmkl_gf_lp64.so" PARENT_SCOPE)
  elseif(output_var STREQUAL "_MPTK_MKL_SEQUENTIAL_LIBRARY")
    set(${output_var} "/mock/libmkl_sequential.so" PARENT_SCOPE)
  elseif(output_var STREQUAL "_MPTK_MKL_CORE_LIBRARY")
    set(${output_var} "/mock/libmkl_core.so" PARENT_SCOPE)
  else()
    message(FATAL_ERROR "Unexpected find_library call for ${output_var}")
  endif()
endfunction()

function(assert_equal actual expected description)
  if(NOT "${actual}" STREQUAL "${expected}")
    message(FATAL_ERROR
      "${description}: expected '${expected}', got '${actual}'")
  endif()
endfunction()

set(MPTK_BLAS_VENDOR auto)
set(CMAKE_CXX_COMPILER_ID GNU)
set(ENV{BLA_VENDOR} "")
unset(BLA_VENDOR)

mptk_prepare_numerical_backend("/mock/libarpack.so")

assert_equal("${MPTK_DETECTED_ARPACK_FORTRAN_RUNTIME}" "gnu"
  "ARPACK runtime detection")
assert_equal("${MPTK_NUMERICAL_BACKEND_TRY_GNU_FORTRAN}" "TRUE"
  "Fortran language selection")
if(DEFINED CACHE{BLAS_LIBRARIES} OR DEFINED CACHE{LAPACK_LIBRARIES})
  message(FATAL_ERROR
    "Default BLA_VENDOR=All unexpectedly bypassed FindBLAS")
endif()

# An explicit non-MKL BLA_VENDOR must retain control of FindBLAS.
set(BLA_VENDOR OpenBLAS)

mptk_prepare_numerical_backend("/mock/libarpack.so")

if(DEFINED CACHE{BLAS_LIBRARIES} OR DEFINED CACHE{LAPACK_LIBRARIES})
  message(FATAL_ERROR
    "BLA_VENDOR=OpenBLAS unexpectedly preselected MKL libraries")
endif()
assert_equal("${MPTK_NUMERICAL_BACKEND_TRY_GNU_FORTRAN}" "FALSE"
  "OpenBLAS Fortran language selection")

# An explicit MKL request retains the compiler-independent GNU-interface
# fallback for systems where FindBLAS cannot infer the interface.
unset(BLA_VENDOR)
set(MPTK_BLAS_VENDOR MKL)

mptk_prepare_numerical_backend("/mock/libarpack.so")

set(_expected_mkl_libraries
  "/mock/libmkl_gf_lp64.so"
  "/mock/libmkl_sequential.so"
  "/mock/libmkl_core.so"
  -lm
  -ldl)
assert_equal("${BLAS_LIBRARIES}" "${_expected_mkl_libraries}"
  "explicit MKL BLAS selection")
assert_equal("${LAPACK_LIBRARIES}" "${_expected_mkl_libraries}"
  "explicit MKL LAPACK selection")
