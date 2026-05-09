cmake_minimum_required(VERSION 3.20)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

project(MPToolkit LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(MPTK_WITH_OPENMP "Enable OpenMP support" OFF)
option(MPTK_WITH_PSTREAM "Enable persistent stream support" ON)
option(MPTK_BUILD_CONTRIB_MODELS "Build models/contrib generators" OFF)
set(MPTK_BLAS_VENDOR "auto" CACHE STRING "BLAS vendor for version reporting: auto, generic, OpenBLAS, or MKL")
set_property(CACHE MPTK_BLAS_VENDOR PROPERTY STRINGS auto generic OpenBLAS MKL)
set(MPTK_FORTRAN_COMPLEX_RETURN "register" CACHE STRING "Fortran complex return convention: register or first_arg")
set_property(CACHE MPTK_FORTRAN_COMPLEX_RETURN PROPERTY STRINGS register first_arg)

include(CheckCXXSourceCompiles)
include(CheckIncludeFileCXX)
include(CheckSymbolExists)
include(CheckTypeSize)
include(TestBigEndian)

find_package(Boost 1.47 REQUIRED COMPONENTS program_options)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_library(ARPACK_LIBRARY NAMES arpack REQUIRED)
find_path(GMP_INCLUDE_DIR NAMES gmp.h REQUIRED)
find_library(GMP_LIBRARY NAMES gmp REQUIRED)

if(NOT TARGET BLAS::BLAS)
  add_library(BLAS::BLAS INTERFACE IMPORTED)
  set_target_properties(BLAS::BLAS PROPERTIES INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}")
endif()

if(NOT TARGET LAPACK::LAPACK)
  add_library(LAPACK::LAPACK INTERFACE IMPORTED)
  set_target_properties(LAPACK::LAPACK PROPERTIES INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}")
endif()

add_library(ARPACK::ARPACK UNKNOWN IMPORTED)
set_target_properties(ARPACK::ARPACK PROPERTIES IMPORTED_LOCATION "${ARPACK_LIBRARY}")

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

if(MPTK_WITH_PSTREAM)
  set(USE_PSTREAM TRUE)
endif()

if(MPTK_FORTRAN_COMPLEX_RETURN STREQUAL "first_arg")
  set(FORTRAN_COMPLEX_RETURN_FIRST_ARG TRUE)
else()
  set(FORTRAN_COMPLEX_RETURN_IN_REGISTER TRUE)
endif()

if(MPTK_BLAS_VENDOR STREQUAL "auto")
  string(TOLOWER "${BLAS_LIBRARIES}" MPTK_BLAS_LIBRARIES_LOWER)
  if(MPTK_BLAS_LIBRARIES_LOWER MATCHES "openblas")
    set(MPTK_CONFIG_BLAS_VENDOR "OpenBLAS")
  elseif(MPTK_BLAS_LIBRARIES_LOWER MATCHES "mkl")
    set(MPTK_CONFIG_BLAS_VENDOR "MKL")
  else()
    set(MPTK_CONFIG_BLAS_VENDOR "generic")
  endif()
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
  expokit/zgpadm.f
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

add_library(mptoolkit_core_objects OBJECT ${MPTK_CORE_SOURCES})
target_link_libraries(mptoolkit_core_objects PUBLIC mptoolkit_config)

add_library(mptoolkit_core STATIC $<TARGET_OBJECTS:mptoolkit_core_objects>)
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
  add_executable("${target}" "${source}"
    $<TARGET_OBJECTS:mptoolkit_core_objects>
    ${ARGN}
  )
  target_link_libraries("${target}" PRIVATE mptoolkit_algorithms)
  mptk_set_runtime_output("${target}" "${subdir}")
endfunction()

set(MPTK_TOOLS
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

add_executable(testexponential linearalgebra/test/testexponential.cpp
  $<TARGET_OBJECTS:mptoolkit_core_objects>
)
target_link_libraries(testexponential PRIVATE mptoolkit_core)
mptk_set_runtime_output(testexponential tests)
add_test(NAME linearalgebra.testexponential COMMAND testexponential)
set_tests_properties(linearalgebra.testexponential PROPERTIES TIMEOUT 30)

add_custom_target(test-integration
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/scripts/mptk-test" --bin-dir "${CMAKE_BINARY_DIR}"
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
  DEPENDS tools models
  USES_TERMINAL
)

install(TARGETS ${MPTK_TOOL_TARGETS} RUNTIME DESTINATION bin)
install(TARGETS ${MPTK_MODEL_TARGETS} RUNTIME DESTINATION bin)
