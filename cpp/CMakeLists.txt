add_library(bandsolver_cpp src/band.cpp src/newton.cpp src/fd.cpp src/factor.cpp)
target_include_directories(bandsolver_cpp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(bandsolver_cpp PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic -Wconversion -Wno-sign-conversion>)
