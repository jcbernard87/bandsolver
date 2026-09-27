# Changelog

## Unreleased

### Added
- Notebook **4 · Finite differences vs finite volumes**: an electrode/separator diffusion problem with a porosity/diffusivity jump, solved with nodal FD, FD with interface flux matching, and FV. It compares conservation (FV to round-off; nodal FD +8 % salt) and spatial order (FV 2, FD-matched 1, nodal FD 0), shows FV on a graded mesh (38–49× lower error at equal node count), shows FD and FV are identical for smooth single-material problems, and gives a decision guide.

## 0.1.1 — 2026-09-27

### Added
- **Finite-difference Jacobians** (Appendix C's "AUTOBAND" idea). `fd_jacobian`, `newton_fd` and `check_jacobian` are available in C++, Fortran, the C ABI and Python (both backends).
  - A block Jacobian, including the X/Y endpoint blocks, costs 3n + 1 residual evaluations regardless of the grid size.
  - `check_jacobian` reports the worst mismatch per block with its node, row and column.
- **Tutorial notebooks** (`notebooks/`), each with LaTeX derivations, implementation and validation against analytic solutions:
  - getting started;
  - coupled transient Nernst–Planck transport in a binary electrolyte;
  - the Newman–Tobias porous-electrode current distribution with Butler–Volmer kinetics.
- **Prebuilt wheels** for CPython 3.9–3.15 on Linux (x86_64, aarch64), macOS (arm64, x86_64) and Windows (x86_64), with the Fortran runtime bundled.
- **CI** on the same five platforms, including Windows with Intel ifx and MSVC.
- `CITATION.cff`, Zenodo metadata, and `bandsolver.__version__`.

### Fixed
- `band_solve` (Fortran) now validates its arguments before writing its output array.

## 0.1.0 — 2026-09-26

- First public release. It provides the BAND block solver with partial and legacy pivoting, and a Newton driver, in Fortran 2008 (with a C ABI) and C++17, with Python bindings for both.
- The documentation covers the formulation, the API, validation and provenance.
