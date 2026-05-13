cmake_minimum_required(VERSION 3.20)

# Keep the bootstrap Release default close to Autoconf and uni20 without
# clobbering user-provided or toolchain-provided release flags.
set(MPTK_DEFAULT_CXX_FLAGS_RELEASE "-DNDEBUG -O2")
if(NOT DEFINED CMAKE_CXX_FLAGS_RELEASE)
  set(CMAKE_CXX_FLAGS_RELEASE "${MPTK_DEFAULT_CXX_FLAGS_RELEASE}"
      CACHE STRING "Flags used by the CXX compiler during RELEASE builds.")
endif()

project(MPToolkit LANGUAGES C CXX)

if(CMAKE_CONFIGURATION_TYPES)
  message(FATAL_ERROR
    "MPToolkit's CMake bootstrap currently supports only single-config generators. "
    "Use Ninja (recommended) or Unix Makefiles with -DCMAKE_BUILD_TYPE=Release or Debug.")
endif()

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option(MPTK_WITH_OPENMP "Enable OpenMP support" OFF)
option(MPTK_WITH_PSTREAM "Enable persistent stream support" ON)
option(MPTK_BUILD_CONTRIB_MODELS "Build models/contrib generators" OFF)
option(MPTK_ENABLE_IPO "Enable IPO/LTO for Release and MinSizeRel builds when supported" ON)
if(CMAKE_GENERATOR MATCHES "Ninja")
  set(MPTK_DEFAULT_IPO_JOBS auto)
else()
  set(MPTK_DEFAULT_IPO_JOBS 1)
endif()
set(MPTK_IPO_JOBS "${MPTK_DEFAULT_IPO_JOBS}" CACHE STRING "GNU IPO/LTO worker count: positive integer, jobserver, or auto")
set_property(CACHE MPTK_IPO_JOBS PROPERTY STRINGS 1 2 4 8 16 jobserver auto)
set(MPTK_LINK_JOBS "1" CACHE STRING "Ninja-only maximum concurrent link jobs; empty uses the generator default")
set_property(CACHE MPTK_LINK_JOBS PROPERTY STRINGS "" 1 2 4 8 16)
set(MPTK_BLAS_VENDOR "auto" CACHE STRING "BLAS vendor for version reporting: auto, generic, OpenBLAS, or MKL")
set_property(CACHE MPTK_BLAS_VENDOR PROPERTY STRINGS auto generic OpenBLAS MKL)
set(MPTK_BLAS_VENDOR_VALUES auto generic OpenBLAS MKL)
if(NOT MPTK_BLAS_VENDOR IN_LIST MPTK_BLAS_VENDOR_VALUES)
  message(FATAL_ERROR "MPTK_BLAS_VENDOR must be one of: auto, generic, OpenBLAS, MKL")
endif()
if(NOT MPTK_LINK_JOBS STREQUAL "")
  if(NOT MPTK_LINK_JOBS MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR "MPTK_LINK_JOBS must be empty or a positive integer")
  endif()
  if(CMAKE_GENERATOR MATCHES "Ninja")
    set_property(GLOBAL APPEND PROPERTY JOB_POOLS "mptk_link_pool=${MPTK_LINK_JOBS}")
    if(NOT DEFINED CMAKE_JOB_POOL_LINK OR CMAKE_JOB_POOL_LINK STREQUAL "")
      set(CMAKE_JOB_POOL_LINK mptk_link_pool)
    endif()
  endif()
endif()

function(mptk_configure_gnu_ipo_jobs lang)
  if(NOT CMAKE_${lang}_COMPILER_ID STREQUAL "GNU")
    return()
  endif()
  if(NOT MPTK_IPO_JOBS MATCHES "^(auto|jobserver|[1-9][0-9]*)$")
    message(FATAL_ERROR "MPTK_IPO_JOBS must be a positive integer, jobserver, or auto")
  endif()
  if(DEFINED CMAKE_${lang}_COMPILE_OPTIONS_IPO)
    set(_ipo_options ${CMAKE_${lang}_COMPILE_OPTIONS_IPO})
    list(TRANSFORM _ipo_options REPLACE "^-flto(=.*)?$" "-flto=${MPTK_IPO_JOBS}")
    set(CMAKE_${lang}_COMPILE_OPTIONS_IPO ${_ipo_options} PARENT_SCOPE)
  endif()
endfunction()

if(MPTK_ENABLE_IPO)
  mptk_configure_gnu_ipo_jobs(C)
  mptk_configure_gnu_ipo_jobs(CXX)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT MPTK_IPO_SUPPORTED OUTPUT MPTK_IPO_ERROR LANGUAGES C CXX)
  if(MPTK_IPO_SUPPORTED)
    if(NOT DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION AND
       NOT DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE)
      set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    endif()
    if(NOT DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION AND
       NOT DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL)
      set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL ON)
    endif()
  else()
    message(STATUS "IPO/LTO not enabled: ${MPTK_IPO_ERROR}")
  endif()
endif()

set(MPTK_FORTRAN_COMPLEX_RETURN "auto" CACHE STRING "Fortran complex return convention: auto, register, or first_arg")
set_property(CACHE MPTK_FORTRAN_COMPLEX_RETURN PROPERTY STRINGS auto register first_arg)
set(MPTK_FORTRAN_COMPLEX_RETURN_VALUES auto register first_arg)
if(NOT MPTK_FORTRAN_COMPLEX_RETURN IN_LIST MPTK_FORTRAN_COMPLEX_RETURN_VALUES)
  message(FATAL_ERROR "MPTK_FORTRAN_COMPLEX_RETURN must be one of: auto, register, first_arg")
endif()

if(NOT MPTK_WITH_PSTREAM)
  message(FATAL_ERROR "MPTK_WITH_PSTREAM=OFF is not supported by the current source tree; PStream serialization is required")
endif()

include(CheckCXXSourceCompiles)
include(CheckCXXSourceRuns)
include(CheckIncludeFileCXX)
include(CheckSymbolExists)
include(CheckTypeSize)
include(MPTKBlasVendor)
include(MPTKNumericalBackend)
include(TestBigEndian)

find_library(ARPACK_LIBRARY NAMES arpack REQUIRED)
mptk_prepare_numerical_backend("${ARPACK_LIBRARY}")

find_package(Git QUIET)
if(NOT Boost_NO_BOOST_CMAKE)
  find_package(Boost 1.47 CONFIG QUIET COMPONENTS program_options)
endif()
if(NOT Boost_FOUND)
  if(POLICY CMP0167)
    cmake_policy(PUSH)
    cmake_policy(SET CMP0167 OLD)
  endif()
  find_package(Boost 1.47 REQUIRED MODULE COMPONENTS program_options)
  if(POLICY CMP0167)
    cmake_policy(POP)
  endif()
endif()
if(NOT TARGET Boost::program_options)
  message(FATAL_ERROR "Boost program_options target not found")
endif()
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_path(GMP_INCLUDE_DIR NAMES gmp.h REQUIRED)
find_library(GMP_LIBRARY NAMES gmp REQUIRED)

mptk_detect_blas_vendor(MPTK_DETECTED_BLAS_VENDOR)
mptk_validate_numerical_backend("${ARPACK_LIBRARY}")

if(NOT TARGET BLAS::BLAS)
  add_library(BLAS::BLAS INTERFACE IMPORTED)
  set_target_properties(BLAS::BLAS PROPERTIES INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}")
endif()

if(NOT TARGET LAPACK::LAPACK)
  add_library(LAPACK::LAPACK INTERFACE IMPORTED)
  set_target_properties(LAPACK::LAPACK PROPERTIES INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}")
endif()

add_library(ARPACK::ARPACK UNKNOWN IMPORTED)
set_target_properties(ARPACK::ARPACK PROPERTIES
  IMPORTED_LOCATION "${ARPACK_LIBRARY}"
  INTERFACE_LINK_LIBRARIES "LAPACK::LAPACK;BLAS::BLAS"
)

add_library(GMP::GMP UNKNOWN IMPORTED)
set_target_properties(GMP::GMP PROPERTIES
  IMPORTED_LOCATION "${GMP_LIBRARY}"
  INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}"
)

if(MPTK_WITH_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
  set(HAVE_OPENMP TRUE)
  set(MULTITHREAD TRUE)
  set(MPTK_OPENMP_CXXFLAGS "${OpenMP_CXX_FLAGS}")
else()
  set(MPTK_OPENMP_CXXFLAGS "")
endif()

check_include_file_cxx(fcntl.h HAVE_FCNTL_H)
check_include_file_cxx(inttypes.h HAVE_INTTYPES_H)
check_include_file_cxx(stdint.h HAVE_STDINT_H)
check_include_file_cxx(stdio.h HAVE_STDIO_H)
check_include_file_cxx(stdlib.h HAVE_STDLIB_H)
check_include_file_cxx(strings.h HAVE_STRINGS_H)
check_include_file_cxx(string.h HAVE_STRING_H)
check_include_file_cxx(sys/ioctl.h HAVE_SYS_IOCTL_H)
check_include_file_cxx(sys/stat.h HAVE_SYS_STAT_H)
check_include_file_cxx(sys/types.h HAVE_SYS_TYPES_H)
check_include_file_cxx(termios.h HAVE_TERMIOS_H)
check_include_file_cxx(unistd.h HAVE_UNISTD_H)

check_symbol_exists(gettimeofday "sys/time.h" HAVE_GETTIMEOFDAY)
check_symbol_exists(strerror_r "string.h" HAVE_STRERROR_R)
check_symbol_exists(TIOCGWINSZ "sys/ioctl.h" GWINSZ_IN_SYS_IOCTL)

if(HAVE_STRERROR_R)
  set(HAVE_DECL_STRERROR_R TRUE)
endif()

check_cxx_source_compiles("
  #include <string.h>
  int main()
  {
    char buffer[256];
    char* value = strerror_r(0, buffer, sizeof(buffer));
    (void)value;
    return 0;
  }
" STRERROR_R_CHAR_P)

check_type_size("short" SIZEOF_SHORT LANGUAGE CXX)
check_type_size("int" SIZEOF_INT LANGUAGE CXX)
check_type_size("long" SIZEOF_LONG LANGUAGE CXX)
check_type_size("long long" SIZEOF_LONG_LONG LANGUAGE CXX)
check_type_size("float" SIZEOF_FLOAT LANGUAGE CXX)
check_type_size("double" SIZEOF_DOUBLE LANGUAGE CXX)

test_big_endian(MPTK_WORDS_BIGENDIAN)
if(MPTK_WORDS_BIGENDIAN)
  set(WORDS_BIGENDIAN TRUE)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang|Intel")
  set(HAVE_PRETTYFUNC TRUE)
endif()
string(TOLOWER "${CMAKE_CXX_COMPILER_ID}" MPTK_CONFIG_COMPILER_VENDOR)

if(MPTK_WITH_PSTREAM)
  set(USE_PSTREAM TRUE)
endif()

if(MPTK_FORTRAN_COMPLEX_RETURN STREQUAL "auto")
  if(CMAKE_CROSSCOMPILING)
    message(FATAL_ERROR "MPTK_FORTRAN_COMPLEX_RETURN=auto cannot run the BLAS probe while cross compiling; set it to register or first_arg")
  endif()

  set(MPTK_REQUIRED_INCLUDES_SAVE "${CMAKE_REQUIRED_INCLUDES}")
  set(MPTK_REQUIRED_LIBRARIES_SAVE "${CMAKE_REQUIRED_LIBRARIES}")
  set(CMAKE_REQUIRED_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}")
  set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES})

  check_cxx_source_runs("
    #define FORTRAN_COMPLEX_RETURN_FIRST_ARG
    #include \"common/blas1f.h\"
    #include <complex>
    #include <vector>
    using complex = std::complex<double>;
    int main()
    {
      std::vector<complex> v1(10, complex(2.0, 2.0));
      std::vector<complex> v2(10, 1.0);
      complex res = BLAS::zdotu(10, v1.data(), 1, v2.data(), 1);
      return std::abs(res - complex(20.0, 20.0)) < 1E-10 ? 0 : 1;
    }
  " MPTK_FORTRAN_COMPLEX_RETURN_FIRST_ARG_WORKS)

  check_cxx_source_runs("
    #define FORTRAN_COMPLEX_RETURN_IN_REGISTER
    #include \"common/blas1f.h\"
    #include <complex>
    #include <vector>
    using complex = std::complex<double>;
    int main()
    {
      std::vector<complex> v1(10, complex(2.0, 2.0));
      std::vector<complex> v2(10, 1.0);
      complex res = BLAS::zdotu(10, v1.data(), 1, v2.data(), 1);
      return std::abs(res - complex(20.0, 20.0)) < 1E-10 ? 0 : 1;
    }
  " MPTK_FORTRAN_COMPLEX_RETURN_IN_REGISTER_WORKS)

  set(CMAKE_REQUIRED_INCLUDES "${MPTK_REQUIRED_INCLUDES_SAVE}")
  set(CMAKE_REQUIRED_LIBRARIES "${MPTK_REQUIRED_LIBRARIES_SAVE}")

  if(MPTK_FORTRAN_COMPLEX_RETURN_FIRST_ARG_WORKS)
    set(MPTK_FORTRAN_COMPLEX_RETURN_RESOLVED "first_arg")
  elseif(MPTK_FORTRAN_COMPLEX_RETURN_IN_REGISTER_WORKS)
    set(MPTK_FORTRAN_COMPLEX_RETURN_RESOLVED "register")
  else()
    message(FATAL_ERROR "Cannot determine how Fortran returns complex BLAS values; set MPTK_FORTRAN_COMPLEX_RETURN explicitly")
  endif()
else()
  set(MPTK_FORTRAN_COMPLEX_RETURN_RESOLVED "${MPTK_FORTRAN_COMPLEX_RETURN}")
endif()

if(MPTK_FORTRAN_COMPLEX_RETURN_RESOLVED STREQUAL "first_arg")
  set(FORTRAN_COMPLEX_RETURN_FIRST_ARG TRUE)
else()
  set(FORTRAN_COMPLEX_RETURN_IN_REGISTER TRUE)
endif()

if(MPTK_BLAS_VENDOR STREQUAL "auto")
  set(MPTK_CONFIG_BLAS_VENDOR "${MPTK_DETECTED_BLAS_VENDOR}")
else()
  set(MPTK_CONFIG_BLAS_VENDOR "${MPTK_BLAS_VENDOR}")
endif()

if(MPTK_CONFIG_BLAS_VENDOR STREQUAL "OpenBLAS")
  set(BLAS_VENDOR_OPENBLAS TRUE)
elseif(MPTK_CONFIG_BLAS_VENDOR STREQUAL "MKL")
  set(BLAS_VENDOR_MKL TRUE)
endif()

string(TOUPPER "${CMAKE_BUILD_TYPE}" MPTK_BUILD_TYPE_UPPER)
if(MPTK_BUILD_TYPE_UPPER STREQUAL "DEBUG")
  set(CONFIG_DEBUG TRUE)
else()
  set(CONFIG_NODEBUG TRUE)
endif()
set(MPTK_CONFIG_CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${MPTK_BUILD_TYPE_UPPER}}")
string(STRIP "${MPTK_CONFIG_CXXFLAGS}" MPTK_CONFIG_CXXFLAGS)

set(MPTK_GIT_VERSION_ESCAPED "")
if(Git_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse --symbolic-full-name --abbrev-ref HEAD
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE MPTK_GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" describe --dirty --always --tags
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE MPTK_GIT_TAG_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if(MPTK_GIT_BRANCH AND MPTK_GIT_TAG_VERSION)
    set(MPTK_GIT_VERSION "${MPTK_GIT_BRANCH}:${MPTK_GIT_TAG_VERSION}")
    string(REPLACE "\\" "\\\\" MPTK_GIT_VERSION_ESCAPED "${MPTK_GIT_VERSION}")
    string(REPLACE "\"" "\\\"" MPTK_GIT_VERSION_ESCAPED "${MPTK_GIT_VERSION_ESCAPED}")
  endif()
endif()

configure_file(cmake/config.h.in config.h @ONLY)

set(MPTK_CORE_SOURCES
  common/blas_vendor.cpp
  common/conflist.cpp
  common/halfint.cpp
  common/niftycounter.cpp
  common/openmp.cpp
  common/poolallocator.cpp
  common/proccontrol.cpp
  common/randutil.cpp
  common/sha256.cpp
  common/terminal.cpp
  common/unique.cpp
  interface/attributes.cpp
  interface/history.cpp
  interface/inittemp.cpp
  lattice/function.cpp
  lattice/infinite-parser.cpp
  lattice/infinitelattice.cpp
  lattice/latticesite.cpp
  lattice/siteoperator-parser.cpp
  lattice/siteoperator.cpp
  lattice/unitcell-parser.cpp
  lattice/unitcell.cpp
  lattice/unitcell_mpo.cpp
  linearalgebra/eigen.cpp
  linearalgebra/exponential.cpp
  mp/copyright.cpp
  mp-algorithms/transfer.cpp
  mpo/basic_finite_mpo.cpp
  mpo/basic_triangular_mpo.cpp
  mpo/e-optim.cpp
  mpo/f-optim.cpp
  mpo/generic_mpo.cpp
  mpo/infinite_mpo.cpp
  mpo/operator-diagonal.cpp
  mpo/operator_component.cpp
  mpo/operator_utilities.cpp
  mpo/product_mpo.cpp
  mps/density.cpp
  mps/packunpack.cpp
  mps/state_component.cpp
  parser/parser.cpp
  pheap/bufferalloc.cpp
  pheap/pagefile.cpp
  pheap/pagestream.cpp
  pheap/pheap.cpp
  pheap/pheapallocator.cpp
  pheap/pheapfsblock.cpp
  pheap/pheapfsv4.cpp
  pheap/pheapstream.cpp
  pheap/rawpagestream.cpp
  pstream/pfilestream.cpp
  pstream/pstream.cpp
  quantumnumbers/coupling.cpp
  quantumnumbers/null-quantumnumber.cpp
  quantumnumbers/quantumnumber.cpp
  quantumnumbers/su2.cpp
  quantumnumbers/symmetrybase.cpp
  quantumnumbers/symmetrylist.cpp
  quantumnumbers/u1.cpp
  quantumnumbers/z2.cpp
  quantumnumbers/zn.cpp
  tensor/basis.cpp
  tensor/regularize.cpp
  tensor/tensor_exponential.cpp
  tensor/tensor_eigen.cpp
  tensor/tensorproduct.cpp
  tensor/tensorsum.cpp
  wavefunction/canonicalwavefunction.cpp
  wavefunction/ea.cpp
  wavefunction/finitewavefunctionleft.cpp
  wavefunction/ibc.cpp
  wavefunction/infinitewavefunctionleft.cpp
  wavefunction/infinitewavefunctionright.cpp
  wavefunction/linearwavefunction.cpp
  wavefunction/mpwavefunction.cpp
  wavefunction/operator_actions.cpp
)

set(MPTK_ALGORITHM_SOURCES
  common/rangelist.cpp
  mp-algorithms/dmrg.cpp
  mp-algorithms/ef-matrix.cpp
  mp-algorithms/eigensolver.cpp
  mp-algorithms/excitation-ansatz.cpp
  mp-algorithms/expansion.cpp
  mp-algorithms/finitedmrg.cpp
  mp-algorithms/ibc-dmrg.cpp
  mp-algorithms/ibc-tdvp.cpp
  mp-algorithms/itdvp.cpp
  mp-algorithms/random_wavefunc.cpp
  mp-algorithms/stateslist.cpp
  mp-algorithms/tdvp.cpp
  mp-algorithms/tebd.cpp
  mp-algorithms/triangular_mpo_solver.cpp
  mp-algorithms/triangular_mpo_solver_cross.cpp
  mp-algorithms/triangular_mpo_solver_ea.cpp
  mp-algorithms/triangular_mpo_solver_helpers.cpp
  mp-algorithms/triangular_mpo_solver_simple.cpp
  mp-algorithms/write-matrix.cpp
  mpo/ef-optim.cpp
  parser/matrix-parser.cpp
  parser/number-parser.cpp
  tensor/wigner_eckart.cpp
  wavefunction/momentum_operations.cpp
)

add_library(mptoolkit_config INTERFACE)
target_compile_features(mptoolkit_config INTERFACE cxx_std_20)
target_compile_definitions(mptoolkit_config INTERFACE HAVE_CONFIG_H)
target_compile_definitions(mptoolkit_config INTERFACE _FILE_OFFSET_BITS=64)
if(MPTK_GIT_VERSION_ESCAPED)
  target_compile_definitions(mptoolkit_config INTERFACE "GIT_VERSION_STRING=\"${MPTK_GIT_VERSION_ESCAPED}\"")
endif()
target_include_directories(mptoolkit_config INTERFACE
  "${CMAKE_CURRENT_BINARY_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}"
)
target_link_libraries(mptoolkit_config INTERFACE
  Boost::program_options
  GMP::GMP
  LAPACK::LAPACK
  BLAS::BLAS
  ARPACK::ARPACK
)
if(MPTK_WITH_OPENMP)
  target_link_libraries(mptoolkit_config INTERFACE OpenMP::OpenMP_CXX)
endif()

set(MPTK_SYMMETRY_REGISTRAR_SOURCES
  quantumnumbers/null-quantumnumber.cpp
  quantumnumbers/su2.cpp
  quantumnumbers/u1.cpp
  quantumnumbers/z2.cpp
  quantumnumbers/zn.cpp
)
list(REMOVE_ITEM MPTK_CORE_SOURCES ${MPTK_SYMMETRY_REGISTRAR_SOURCES})

add_library(mptoolkit_core_objects OBJECT ${MPTK_CORE_SOURCES})
target_link_libraries(mptoolkit_core_objects PUBLIC mptoolkit_config)

add_library(mptoolkit_symmetry_registrar_objects OBJECT ${MPTK_SYMMETRY_REGISTRAR_SOURCES})
target_link_libraries(mptoolkit_symmetry_registrar_objects PUBLIC mptoolkit_config)

add_library(mptoolkit_core STATIC
  $<TARGET_OBJECTS:mptoolkit_core_objects>
  $<TARGET_OBJECTS:mptoolkit_symmetry_registrar_objects>
)
target_link_libraries(mptoolkit_core PUBLIC mptoolkit_config)

add_library(mptoolkit_algorithms_objects OBJECT ${MPTK_ALGORITHM_SOURCES})
target_link_libraries(mptoolkit_algorithms_objects PUBLIC mptoolkit_config)

add_library(mptoolkit_algorithms STATIC $<TARGET_OBJECTS:mptoolkit_algorithms_objects>)
target_link_libraries(mptoolkit_algorithms PUBLIC mptoolkit_core)

function(mptk_set_runtime_output target subdir)
  set(output_dir "${CMAKE_BINARY_DIR}/bin/${subdir}")
  set_target_properties("${target}" PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${output_dir}"
  )
  foreach(config Debug Release RelWithDebInfo MinSizeRel)
    string(TOUPPER "${config}" config_upper)
    set_target_properties("${target}" PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY_${config_upper} "${output_dir}/${config}"
    )
  endforeach()
endfunction()

function(mptk_add_executable target source subdir)
  add_executable("${target}" "${source}" ${ARGN})
  # Symmetry registration relies on static initializers in these translation
  # units. Link only those objects directly; the rest of core remains a normal
  # archive so unused object files can be discarded by the linker.
  target_sources("${target}" PRIVATE $<TARGET_OBJECTS:mptoolkit_symmetry_registrar_objects>)
  target_link_libraries("${target}" PRIVATE mptoolkit_algorithms)
  mptk_set_runtime_output("${target}" "${subdir}")
endfunction()

set(MPTK_TOOLS
  mp-idmrg
  mp-idmrg-s3e
  mp-ioverlap
  mp-info
  mp-history
  mp-reflect
  mp-conj
  mp-imoments
  mp-wigner-eckart
  mp-reorder-symmetry
  mp-ispectrum
  mp-lattice-info
  mp-show-operator
  mp-aux-matrix
  mp-aux-algebra
  mp-iapply
  mp-attr
  mp-idivide
  mp-ies
  mp-ies-cross
  mp-irepeat
  mp-norm
  mp-allcorrelation
  mp-irotate
  mp-icorrelation
  mp-itebd
  mp-coarsegrain
  mp-finegrain
  mp-dmrg
  mp-dmrg-3s
  mp-dmrg-2site
  mp-random
  mp-expectation
  mp-overlap
  mp-scale
  mp-matrix
  mp-tebd
  mp-apply
  mp-normalize
  mp-tdvp
  mp-itdvp
  mp-iexpectation-cross
  mp-imoments-cross
  mp-right-canonicalize
  mp-left-canonicalize
  mp-construct
  mp-change-lattice
  mp-finite-create
  mp-ibc-create
  mp-ibc-tdvp
  mp-ibc-apply
  mp-ibc-overlap
  mp-ibc-correlation
  mp-excitation-ansatz
  mp-ibc-wavepacket
  mp-ibc-splice
  mp-ibc-dmrg
  mp-ea-moments
  mp-ea-create
  mp-ea-extend
  mp-ea-change-k
  mp-ea-dmrg
)

set(MPTK_TOOL_TARGETS)
foreach(tool IN LISTS MPTK_TOOLS)
  set(extra_sources)
  if(tool STREQUAL "mp-ea-dmrg")
    list(APPEND extra_sources mp-algorithms/ea-dmrg.cpp)
  endif()
  if(tool MATCHES "^mp-(tdvp|itdvp|ibc-tdvp|ibc-dmrg)$")
    list(APPEND extra_sources mp-algorithms/tdvp-compositions.cpp)
  endif()
  mptk_add_executable("${tool}" "mp/${tool}.cpp" tools ${extra_sources})
  list(APPEND MPTK_TOOL_TARGETS "${tool}")
endforeach()
add_custom_target(tools DEPENDS ${MPTK_TOOL_TARGETS})

file(GLOB MPTK_DEFAULT_MODEL_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/models/*.cpp")
set(MPTK_MODEL_TARGETS)
foreach(model_source IN LISTS MPTK_DEFAULT_MODEL_SOURCES)
  get_filename_component(model_name "${model_source}" NAME_WE)
  mptk_add_executable("${model_name}" "${model_source}" models)
  list(APPEND MPTK_MODEL_TARGETS "${model_name}")
endforeach()
add_custom_target(models DEPENDS ${MPTK_MODEL_TARGETS})

if(MPTK_BUILD_CONTRIB_MODELS)
  file(GLOB MPTK_CONTRIB_MODEL_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/models/contrib/*.cpp")
  set(MPTK_CONTRIB_MODEL_TARGETS)
  foreach(model_source IN LISTS MPTK_CONTRIB_MODEL_SOURCES)
    get_filename_component(model_name "${model_source}" NAME_WE)
    mptk_add_executable("${model_name}" "${model_source}" models/contrib)
    list(APPEND MPTK_CONTRIB_MODEL_TARGETS "${model_name}")
  endforeach()
  add_custom_target(contrib-models DEPENDS ${MPTK_CONTRIB_MODEL_TARGETS})
  add_custom_target(all-models DEPENDS models contrib-models)
else()
  add_custom_target(all-models DEPENDS models)
endif()

add_custom_target(mptoolkit DEPENDS tools models)

enable_testing()

add_executable(testexponential linearalgebra/test/testexponential.cpp)
target_link_libraries(testexponential PRIVATE mptoolkit_core)
mptk_set_runtime_output(testexponential tests)
add_test(NAME linearalgebra.testexponential COMMAND testexponential)
set_tests_properties(linearalgebra.testexponential PROPERTIES TIMEOUT 30)

add_custom_target(test-integration
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/scripts/mptk-test"
          --bin-dir "${CMAKE_BINARY_DIR}"
          --config "$<CONFIG>"
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
  DEPENDS tools models
  USES_TERMINAL
)

install(TARGETS ${MPTK_TOOL_TARGETS}
  RUNTIME DESTINATION bin
  COMPONENT tools
)
install(TARGETS ${MPTK_MODEL_TARGETS}
  RUNTIME DESTINATION bin
  COMPONENT models
  EXCLUDE_FROM_ALL
)
if(MPTK_BUILD_CONTRIB_MODELS AND MPTK_CONTRIB_MODEL_TARGETS)
  install(TARGETS ${MPTK_CONTRIB_MODEL_TARGETS}
    RUNTIME DESTINATION bin
    COMPONENT contrib-models
    EXCLUDE_FROM_ALL
  )
endif()

add_custom_target(install-tools
  COMMAND "${CMAKE_COMMAND}" --install "${CMAKE_BINARY_DIR}" --component tools
  DEPENDS tools
  USES_TERMINAL
)
add_custom_target(install-models
  COMMAND "${CMAKE_COMMAND}" --install "${CMAKE_BINARY_DIR}" --component models
  DEPENDS models
  USES_TERMINAL
)
if(MPTK_BUILD_CONTRIB_MODELS)
  add_custom_target(install-contrib-models
    COMMAND "${CMAKE_COMMAND}" --install "${CMAKE_BINARY_DIR}" --component contrib-models
    DEPENDS contrib-models
    USES_TERMINAL
  )
  add_custom_target(install-all-models DEPENDS install-models install-contrib-models)
else()
  add_custom_target(install-all-models DEPENDS install-models)
endif()

function(mptk_print_configuration_summary)
  message(STATUS "")
  message(STATUS "===== MPToolkit Build Options =====")
  get_cmake_property(_mptk_cache_vars CACHE_VARIABLES)
  foreach(_var ${_mptk_cache_vars})
    if(_var MATCHES "^MPTK_" AND NOT _var MATCHES "^MPTK_DETECTED_")
      get_property(_type CACHE ${_var} PROPERTY TYPE)
      if(_type STREQUAL "INTERNAL")
        continue()
      endif()
      get_property(_help CACHE ${_var} PROPERTY HELPSTRING)
      message(STATUS "${_var} = ${${_var}}   # ${_help}")
    endif()
  endforeach()
  message(STATUS "")

  if(CMAKE_BUILD_TYPE)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _config_upper)
  else()
    set(_config_upper "RELEASE")
  endif()

  set(_well_known
    CMAKE_SYSTEM_VERSION
    CMAKE_SYSTEM_PROCESSOR
    CMAKE_BUILD_TYPE
    CMAKE_INSTALL_PREFIX
    BUILD_SHARED_LIBS
    BLA_VENDOR
    CMAKE_CXX_FLAGS
    CMAKE_CXX_FLAGS_${_config_upper}
    CMAKE_INTERPROCEDURAL_OPTIMIZATION_${_config_upper}
    CMAKE_JOB_POOL_LINK
  )

  message(STATUS "===== Common CMake Options =====")
  foreach(_var IN LISTS _well_known)
    if(DEFINED ${_var})
      if(DEFINED CACHE{${_var}})
        get_property(_help CACHE ${_var} PROPERTY HELPSTRING)
        message(STATUS "${_var} = ${${_var}}   # ${_help}")
      else()
        message(STATUS "${_var} = ${${_var}}")
      endif()
    else()
      message(STATUS "${_var} = # (not set)")
    endif()
  endforeach()
  message(STATUS "")

  message(STATUS "===== Detected Environment =====")
  foreach(_var ${_mptk_cache_vars})
    if(_var MATCHES "^MPTK_DETECTED_")
      get_property(_help CACHE ${_var} PROPERTY HELPSTRING)
      message(STATUS "${_var} = ${${_var}}   # ${_help}")
    endif()
  endforeach()
  message(STATUS "")
endfunction()

mptk_print_configuration_summary()
