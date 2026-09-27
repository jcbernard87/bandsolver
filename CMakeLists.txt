cmake_minimum_required(VERSION 3.20)
project(bandsolver VERSION 0.1.0 LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

option(BANDSOLVER_BUILD_TESTS "Build Fortran/C++ tests (includes the legacy oracle)" ON)
option(BANDSOLVER_BUILD_EXAMPLES "Build the Fortran/C++ quickstart programs" ON)
option(BANDSOLVER_BUILD_PYTHON "Build the pybind11 extension bandsolver._core" OFF)
# Private, non-distributed copy of the archival BAND/MATINV used only as a test oracle
# (see docs/provenance.md). Legacy comparison tests are built only when it is present.
set(BANDSOLVER_LEGACY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../bandsolver-legacy-oracle"
    CACHE PATH "Directory holding the private archival-kernel test oracle")

if(BANDSOLVER_BUILD_TESTS)
  enable_testing()
endif()

add_subdirectory(fortran)
add_subdirectory(cpp)
if(BANDSOLVER_BUILD_PYTHON)
  add_subdirectory(python)
endif()
if(BANDSOLVER_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(BANDSOLVER_BUILD_TESTS)
  add_subdirectory(tests)
endif()
