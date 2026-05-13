#[[
Configure and validate the numerical backend ABI used by BLAS/LAPACK
and external Fortran libraries.

The CMake bootstrap currently supports the same LP64 integer ABI as the
existing source tree.  BLAS/LAPACK symbols are process-global on typical ELF
systems, so dependencies such as ARPACK must be compatible with the selected
provider and Fortran interface.
]]

function(mptk_prepare_numerical_backend)
  set(_arpack_library "")
  if(ARGC GREATER 0)
    set(_arpack_library "${ARGV0}")
  endif()

  set(MPTK_BLAS_INTEGER "lp64" CACHE STRING
      "BLAS/LAPACK integer ABI. The CMake bootstrap currently supports only lp64.")
  set_property(CACHE MPTK_BLAS_INTEGER PROPERTY STRINGS lp64)
  if(NOT MPTK_BLAS_INTEGER STREQUAL "lp64")
    message(FATAL_ERROR
      "MPTK_BLAS_INTEGER=${MPTK_BLAS_INTEGER} is not supported by this bootstrap. "
      "MPToolkit currently uses 32-bit Fortran integers; configure with MPTK_BLAS_INTEGER=lp64.")
  endif()
  set(BLA_SIZEOF_INTEGER 4 CACHE STRING "BLAS/LAPACK integer size in bytes")
  if(NOT BLA_SIZEOF_INTEGER EQUAL 4)
    message(FATAL_ERROR
      "BLA_SIZEOF_INTEGER=${BLA_SIZEOF_INTEGER} selects an ILP64-style BLAS/LAPACK ABI, "
      "but MPToolkit currently uses 32-bit Fortran integers. Use LP64 / BLA_SIZEOF_INTEGER=4.")
  endif()

  set(MPTK_MKL_INTERFACE "auto" CACHE STRING
      "MKL Fortran interface: auto, gnu, or intel")
  set_property(CACHE MPTK_MKL_INTERFACE PROPERTY STRINGS auto gnu intel)
  set(MPTK_MKL_INTERFACE_VALUES auto gnu intel)
  if(NOT MPTK_MKL_INTERFACE IN_LIST MPTK_MKL_INTERFACE_VALUES)
    message(FATAL_ERROR "MPTK_MKL_INTERFACE must be one of: auto, gnu, intel")
  endif()

  option(MPTK_VALIDATE_NUMERICAL_BACKEND
    "Validate BLAS/LAPACK/ARPACK ABI compatibility at configure time" ON)

  set(_arpack_fortran_runtime "")
  if(MPTK_VALIDATE_NUMERICAL_BACKEND AND _arpack_library)
    mptk_detect_shared_library_dependencies("${_arpack_library}" _arpack_deps)
    set(MPTK_DETECTED_ARPACK_DEPENDENCIES "${_arpack_deps}" CACHE INTERNAL
        "Shared-library dependencies detected for ARPACK" FORCE)
    mptk_detect_fortran_runtime(_arpack_deps _arpack_fortran_runtime)
    set(MPTK_DETECTED_ARPACK_FORTRAN_RUNTIME "${_arpack_fortran_runtime}" CACHE INTERNAL
        "Fortran runtime detected from ARPACK dependencies" FORCE)
  endif()

  set(_mptk_mkl_requested FALSE)
  if(MPTK_BLAS_VENDOR STREQUAL "MKL")
    set(_mptk_mkl_requested TRUE)
  endif()
  if(DEFINED BLA_VENDOR AND BLA_VENDOR MATCHES "^Intel")
    set(_mptk_mkl_requested TRUE)
  endif()
  if(DEFINED BLAS_LIBRARIES AND "${BLAS_LIBRARIES}" MATCHES "mkl")
    set(_mptk_mkl_requested TRUE)
  endif()

  set(_mptk_prefer_mkl_gnu FALSE)
  if(MPTK_MKL_INTERFACE STREQUAL "gnu")
    set(_mptk_prefer_mkl_gnu TRUE)
  elseif(MPTK_MKL_INTERFACE STREQUAL "intel")
    set(_mptk_prefer_mkl_gnu FALSE)
  elseif(_arpack_fortran_runtime STREQUAL "gnu")
    set(_mptk_prefer_mkl_gnu TRUE)
  elseif(_arpack_fortran_runtime STREQUAL "intel")
    set(_mptk_prefer_mkl_gnu FALSE)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    set(_mptk_prefer_mkl_gnu TRUE)
  endif()

  if(_mptk_mkl_requested AND _mptk_prefer_mkl_gnu
      AND NOT CMAKE_Fortran_COMPILER_LOADED)
    include(CheckLanguage)
    check_language(Fortran)
    if(CMAKE_Fortran_COMPILER)
      enable_language(Fortran)
      set(MPTK_DETECTED_FORTRAN_COMPILER_ID "${CMAKE_Fortran_COMPILER_ID}"
          CACHE INTERNAL "Detected Fortran compiler ID used for numerical dependency selection" FORCE)
      message(STATUS
        "Enabled Fortran (${CMAKE_Fortran_COMPILER_ID}) so FindBLAS can select the matching MKL interface")
    else()
      set(MPTK_DETECTED_FORTRAN_COMPILER_ID "None"
          CACHE INTERNAL "Detected Fortran compiler ID used for numerical dependency selection" FORCE)
    endif()
  elseif(CMAKE_Fortran_COMPILER_LOADED)
    set(MPTK_DETECTED_FORTRAN_COMPILER_ID "${CMAKE_Fortran_COMPILER_ID}"
        CACHE INTERNAL "Detected Fortran compiler ID used for numerical dependency selection" FORCE)
  endif()

  if(_mptk_mkl_requested AND _mptk_prefer_mkl_gnu
      AND NOT DEFINED CACHE{BLAS_LIBRARIES}
      AND NOT DEFINED CACHE{LAPACK_LIBRARIES})
    set(_mptk_mkl_hints "")
    if(DEFINED ENV{MKLROOT})
      file(TO_CMAKE_PATH "$ENV{MKLROOT}" _mptk_mklroot)
      list(APPEND _mptk_mkl_hints
        "${_mptk_mklroot}/lib/intel64"
        "${_mptk_mklroot}/lib")
    endif()

    find_library(_MPTK_MKL_GF_LP64_LIBRARY NAMES mkl_gf_lp64
      HINTS ${_mptk_mkl_hints})
    find_library(_MPTK_MKL_SEQUENTIAL_LIBRARY NAMES mkl_sequential
      HINTS ${_mptk_mkl_hints})
    find_library(_MPTK_MKL_CORE_LIBRARY NAMES mkl_core
      HINTS ${_mptk_mkl_hints})

    if(_MPTK_MKL_GF_LP64_LIBRARY AND _MPTK_MKL_SEQUENTIAL_LIBRARY AND _MPTK_MKL_CORE_LIBRARY)
      set(_mptk_mkl_gnu_libraries
        "${_MPTK_MKL_GF_LP64_LIBRARY}"
        "${_MPTK_MKL_SEQUENTIAL_LIBRARY}"
        "${_MPTK_MKL_CORE_LIBRARY}"
        -lm
        -ldl)
      set(BLAS_LIBRARIES "${_mptk_mkl_gnu_libraries}" CACHE STRING
          "BLAS libraries" FORCE)
      set(LAPACK_LIBRARIES "${_mptk_mkl_gnu_libraries}" CACHE STRING
          "LAPACK libraries" FORCE)
      set(MPTK_DETECTED_MKL_INTERFACE_DEFAULT "gnu" CACHE INTERNAL
          "Default MKL interface selected before FindBLAS" FORCE)
      message(STATUS
        "Defaulting MKL to GNU Fortran interface: ${_MPTK_MKL_GF_LP64_LIBRARY}")
    elseif(MPTK_MKL_INTERFACE STREQUAL "gnu")
      message(FATAL_ERROR
        "MPTK_MKL_INTERFACE=gnu was requested, but libmkl_gf_lp64, "
        "libmkl_sequential, and libmkl_core could not all be found. "
        "Source the MKL environment or set BLAS_LIBRARIES/LAPACK_LIBRARIES explicitly.")
    endif()
  endif()
endfunction()

function(mptk_detect_mkl_interface output_var)
  set(_interface "")
  foreach(_lib IN LISTS BLAS_LIBRARIES LAPACK_LIBRARIES)
    get_filename_component(_name "${_lib}" NAME)
    string(TOLOWER "${_name}" _name_lower)
    if(_name_lower MATCHES "mkl_gf_")
      set(_interface "gnu")
      break()
    elseif(_name_lower MATCHES "mkl_intel_")
      set(_interface "intel")
    elseif(_name_lower MATCHES "mkl_rt")
      set(_interface "runtime")
    endif()
  endforeach()

  set(${output_var} "${_interface}" PARENT_SCOPE)
  set(MPTK_DETECTED_MKL_INTERFACE "${_interface}" CACHE INTERNAL
      "Detected MKL Fortran interface from BLAS/LAPACK libraries" FORCE)
endfunction()

function(mptk_detect_shared_library_dependencies library output_var)
  set(_deps "")
  if(EXISTS "${library}" AND CMAKE_OBJDUMP AND NOT "${library}" MATCHES "\\.a$")
    execute_process(
      COMMAND "${CMAKE_OBJDUMP}" -p "${library}"
      RESULT_VARIABLE _objdump_result
      OUTPUT_VARIABLE _objdump_output
      ERROR_QUIET
    )
    if(_objdump_result EQUAL 0)
      string(REGEX MATCHALL "NEEDED[ \t]+[^ \t\r\n]+" _dep_lines "${_objdump_output}")
      foreach(_line IN LISTS _dep_lines)
        string(REGEX REPLACE ".*NEEDED[ \t]+([^ \t\r\n]+).*" "\\1" _dep "${_line}")
        list(APPEND _deps "${_dep}")
      endforeach()
      list(REMOVE_DUPLICATES _deps)
    endif()
  endif()

  set(${output_var} "${_deps}" PARENT_SCOPE)
endfunction()

function(mptk_detect_fortran_runtime deps_var output_var)
  set(_runtime "")
  foreach(_dep IN LISTS ${deps_var})
    string(TOLOWER "${_dep}" _dep_lower)
    if(_dep_lower MATCHES "libgfortran")
      set(_runtime "gnu")
      break()
    elseif(_dep_lower MATCHES "libifcore|libifport|libimf|libirc")
      set(_runtime "intel")
      break()
    endif()
  endforeach()

  set(${output_var} "${_runtime}" PARENT_SCOPE)
endfunction()

function(mptk_detect_libraries_fortran_runtime output_var)
  set(_runtime "")
  foreach(_lib IN LISTS ARGN)
    if(NOT _lib OR _lib MATCHES "^-")
      continue()
    endif()

    set(_candidate "")
    get_filename_component(_name "${_lib}" NAME)
    string(TOLOWER "${_name}" _name_lower)
    if(_name_lower MATCHES "mkl_gf_|libgfortran")
      set(_candidate "gnu")
    elseif(_name_lower MATCHES "mkl_intel_|libifcore|libifport|libimf|libirc")
      set(_candidate "intel")
    elseif(EXISTS "${_lib}")
      mptk_detect_shared_library_dependencies("${_lib}" _deps)
      mptk_detect_fortran_runtime(_deps _candidate)
    endif()

    if(_candidate)
      if(_runtime AND NOT _runtime STREQUAL _candidate)
        set(_runtime "mixed")
        break()
      endif()
      set(_runtime "${_candidate}")
    endif()
  endforeach()

  set(${output_var} "${_runtime}" PARENT_SCOPE)
endfunction()

function(mptk_validate_numerical_backend arpack_library)
  if(NOT MPTK_VALIDATE_NUMERICAL_BACKEND)
    return()
  endif()

  if(DEFINED BLA_SIZEOF_INTEGER AND NOT BLA_SIZEOF_INTEGER EQUAL 4)
    message(FATAL_ERROR
      "BLA_SIZEOF_INTEGER=${BLA_SIZEOF_INTEGER} selects an ILP64-style BLAS/LAPACK ABI, "
      "but MPToolkit currently uses 32-bit Fortran integers. Use LP64 / BLA_SIZEOF_INTEGER=4.")
  endif()

  foreach(_lib IN LISTS BLAS_LIBRARIES LAPACK_LIBRARIES)
    get_filename_component(_name "${_lib}" NAME)
    string(TOLOWER "${_name}" _name_lower)
    if(_name_lower MATCHES "ilp64")
      message(FATAL_ERROR
        "Detected ILP64 BLAS/LAPACK library '${_name}', but MPToolkit currently requires LP64.")
    endif()
  endforeach()

  mptk_detect_mkl_interface(_mptk_mkl_interface)
  mptk_detect_shared_library_dependencies("${arpack_library}" _arpack_deps)
  set(MPTK_DETECTED_ARPACK_DEPENDENCIES "${_arpack_deps}" CACHE INTERNAL
      "Shared-library dependencies detected for ARPACK" FORCE)
  mptk_detect_fortran_runtime(_arpack_deps _arpack_fortran_runtime)
  set(MPTK_DETECTED_ARPACK_FORTRAN_RUNTIME "${_arpack_fortran_runtime}" CACHE INTERNAL
      "Fortran runtime detected from ARPACK dependencies" FORCE)

  mptk_detect_libraries_fortran_runtime(_blas_lapack_fortran_runtime ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  set(MPTK_DETECTED_BLAS_LAPACK_FORTRAN_RUNTIME "${_blas_lapack_fortran_runtime}" CACHE INTERNAL
      "Fortran runtime detected from BLAS/LAPACK libraries" FORCE)

  if(_blas_lapack_fortran_runtime STREQUAL "mixed")
    message(FATAL_ERROR
      "Incompatible numerical backend: BLAS/LAPACK libraries appear to mix GNU and Intel "
      "Fortran ABIs (${BLAS_LIBRARIES};${LAPACK_LIBRARIES}). Use libraries built for one ABI.")
  endif()

  if(NOT MPTK_DETECTED_BLAS_VENDOR STREQUAL "MKL")
    if(_arpack_fortran_runtime AND _blas_lapack_fortran_runtime
        AND NOT _arpack_fortran_runtime STREQUAL _blas_lapack_fortran_runtime)
      message(FATAL_ERROR
        "Incompatible numerical backend: ARPACK appears to use the ${_arpack_fortran_runtime} "
        "Fortran ABI, but BLAS/LAPACK appears to use the ${_blas_lapack_fortran_runtime} "
        "Fortran ABI (${BLAS_LIBRARIES};${LAPACK_LIBRARIES}). Use matching BLAS/LAPACK "
        "and ARPACK builds.")
    endif()
    return()
  endif()

  if(MPTK_MKL_INTERFACE STREQUAL "gnu" AND NOT _mptk_mkl_interface STREQUAL "gnu")
    message(FATAL_ERROR
      "MPTK_MKL_INTERFACE=gnu was requested, but BLAS/LAPACK resolved to '${_mptk_mkl_interface}'. "
      "Use MKL libraries containing libmkl_gf_lp64, or adjust MPTK_MKL_INTERFACE.")
  endif()
  if(MPTK_MKL_INTERFACE STREQUAL "intel" AND NOT _mptk_mkl_interface STREQUAL "intel")
    message(FATAL_ERROR
      "MPTK_MKL_INTERFACE=intel was requested, but BLAS/LAPACK resolved to '${_mptk_mkl_interface}'. "
      "Use MKL libraries containing libmkl_intel_lp64, or adjust MPTK_MKL_INTERFACE.")
  endif()

  if(_mptk_mkl_interface STREQUAL "runtime")
    message(WARNING
      "MKL was detected through libmkl_rt; the CMake bootstrap cannot validate the active "
      "Fortran interface against ARPACK.")
    return()
  endif()

  if(_arpack_fortran_runtime STREQUAL "gnu" AND _mptk_mkl_interface STREQUAL "intel")
    message(FATAL_ERROR
      "Incompatible numerical backend: ARPACK depends on the GNU Fortran runtime, "
      "but MKL resolved to the Intel Fortran interface (${BLAS_LIBRARIES}). "
      "Use MKL's GNU interface library libmkl_gf_lp64, for example by enabling a GNU "
      "Fortran compiler for detection or by setting BLAS_LIBRARIES/LAPACK_LIBRARIES explicitly.")
  endif()

  if(_arpack_fortran_runtime STREQUAL "intel" AND _mptk_mkl_interface STREQUAL "gnu")
    message(FATAL_ERROR
      "Incompatible numerical backend: ARPACK appears to use the Intel Fortran runtime, "
      "but MKL resolved to the GNU Fortran interface (${BLAS_LIBRARIES}).")
  endif()
endfunction()
