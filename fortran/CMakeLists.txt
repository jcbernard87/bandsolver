add_library(bandsolver_f
  src/band_kernel.f90
  src/band_newton.f90
  src/band_fd.f90
  src/band_capi.f90)
target_include_directories(bandsolver_f PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(bandsolver_f PRIVATE
  $<$<Fortran_COMPILER_ID:GNU>:-std=f2008 -Wall -Wextra -fimplicit-none -Wno-compare-reals>
  # Intel's default "fast" FP model allows value-changing optimisations; a numerical library
  # needs IEEE-faithful arithmetic (and working NaN/Inf checks).
  $<$<AND:$<Fortran_COMPILER_ID:IntelLLVM,Intel>,$<PLATFORM_ID:Windows>>:/fp:precise>
  $<$<AND:$<Fortran_COMPILER_ID:IntelLLVM,Intel>,$<NOT:$<PLATFORM_ID:Windows>>>:-fp-model=precise>)
