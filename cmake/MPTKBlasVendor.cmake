#[[
Detect and normalize the BLAS vendor after find_package(BLAS).

The Autoconf build currently distinguishes MKL, OpenBLAS, and generic BLAS.
Keep the CMake bootstrap aligned with that public surface while making the
library-name matching less brittle than a whole-path substring search.
]]

function(_mptk_map_blas_vendor_name input_value output_var)
  if(input_value STREQUAL "")
    set(${output_var} "" PARENT_SCOPE)
    return()
  endif()

  # Prefer file basename over full path to avoid false positives from directory
  # names, then fall back to the full token for BLA_VENDOR-style values.
  get_filename_component(_candidate_file "${input_value}" NAME)
  string(TOLOWER "${_candidate_file}" _candidate_file_lower)
  string(TOLOWER "${input_value}" _candidate_path_lower)

  if(_candidate_file_lower MATCHES "mkl")
    set(${output_var} "MKL" PARENT_SCOPE)
    return()
  endif()
  if(_candidate_file_lower MATCHES "openblas")
    set(${output_var} "OpenBLAS" PARENT_SCOPE)
    return()
  endif()

  if(_candidate_path_lower MATCHES "mkl|intel")
    set(${output_var} "MKL" PARENT_SCOPE)
  elseif(_candidate_path_lower MATCHES "openblas")
    set(${output_var} "OpenBLAS" PARENT_SCOPE)
  elseif(_candidate_path_lower MATCHES "generic")
    set(${output_var} "generic" PARENT_SCOPE)
  else()
    set(${output_var} "" PARENT_SCOPE)
  endif()
endfunction()

function(mptk_detect_blas_vendor output_var)
  if(NOT BLAS_FOUND)
    set(_detected_vendor "generic")
  else()
    set(_detected_vendor "")
    foreach(_blas_lib IN LISTS BLAS_LIBRARIES)
      _mptk_map_blas_vendor_name("${_blas_lib}" _mapped_vendor)
      if(NOT _mapped_vendor STREQUAL "")
        set(_detected_vendor "${_mapped_vendor}")
        break()
      endif()
    endforeach()

    if(_detected_vendor STREQUAL "")
      if(DEFINED BLA_VENDOR AND NOT "${BLA_VENDOR}" STREQUAL "" AND NOT "${BLA_VENDOR}" STREQUAL "All")
        _mptk_map_blas_vendor_name("${BLA_VENDOR}" _mapped_vendor)
        if(NOT _mapped_vendor STREQUAL "")
          set(_detected_vendor "${_mapped_vendor}")
        endif()
      endif()
    endif()

    if(_detected_vendor STREQUAL "")
      set(_detected_vendor "generic")
    endif()
  endif()

  set(${output_var} "${_detected_vendor}" PARENT_SCOPE)
  set(MPTK_DETECTED_BLAS_VENDOR "${_detected_vendor}" CACHE INTERNAL "Detected BLAS vendor" FORCE)
  set(MPTK_DETECTED_BLAS_LIBRARIES "${BLAS_LIBRARIES}" CACHE INTERNAL "Detected BLAS libraries" FORCE)
  message(STATUS "Detected MPTK_DETECTED_BLAS_VENDOR: ${MPTK_DETECTED_BLAS_VENDOR}")
endfunction()
