#!/bin/sh
# Install a gfortran that matches the manylinux image's default gcc-toolset, so that
# C, C++ and Fortran come from the same GCC release.
set -eu
if command -v gfortran >/dev/null 2>&1; then gfortran --version | head -1; exit 0; fi
ver=$(gcc -dumpversion | cut -d. -f1)
dnf install -y "gcc-toolset-${ver}-gcc-gfortran" || dnf install -y gcc-gfortran
gfortran --version | head -1
