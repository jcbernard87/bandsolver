# Validation report — bandsolver 0.1.0

**Date:** 2026-09-26.

**Platform:**
- macOS 14 on Apple Silicon (arm64).
- gfortran 14.2.0 (Homebrew GCC) and Apple clang 16.0.0.
- CMake 3.31.2, Release build (`-O3`).
- Python 3.13.3, numpy 2.5.3, pybind11 3.1.0.

**Reproduce:**
```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build   # the 2 legacy-comparison suites run only if the private oracle is present
pip install . --no-build-isolation
pytest
python examples/<name>.py
```

This report covers the generic linear kernel and the Newton driver. It does **not** validate any battery model from the archive. That work is Track 2 and remains gated by `../IMPLEMENTATION_PLAN.md`, Stages 0–4.

## 1. Linear kernel

| Check | Result |
|---|---|
| Backward error ‖KΔc−G‖∞ / (‖K‖∞‖Δc‖∞+‖G‖∞) on random block systems: n ∈ {1,2,5,12}, nj ∈ {3,4,50,500}, with and without X/Y | max 4.0e-16 (Fortran, partial pivot), 3.9e-16 (Fortran, legacy pivot), 5.0e-16 (C++, both modes). Pass threshold 1e-13. |
| Legacy pivot mode vs the frozen archival `BAND`/`MATINV` (private byte-exact extract, not distributed; see `provenance.md`) | **Bit-identical, 0 ulp**, on every case except the one below. |
| Partial pivot mode vs the archival kernel | ≤ 9.6e-16 relative. |
| C++ core vs Fortran library, same inputs | **Bit-identical, 0 ulp**, in both pivot modes. Tests allow 64 ε for other toolchains. |
| Both backends vs numpy dense `solve` of an independently assembled K (Python), n ≤ 6, nj ≤ 200 | ≤ 1e-12 relative; backward error < 1e-14. |
| Singular pivot block at the first, an interior, and the last node | `SINGULAR` status with the correct node, in both modes and all three interfaces (Fortran, C, C++). The Python layer raises `SingularBlockError(node)`. |
| Blocks that need an off-diagonal pivot (B = [[0,1],[1,0]]) | Solved in both modes, backward error < 1e-14. |
| NaN input; nj < 3; unknown pivot; null pointers | `NON_FINITE` / `INVALID_ARGUMENT`; Python raises `NonFiniteError` / `ValueError`. |

**Difference from the archival kernel (a fix, not a regression):** with `nj = 3` and both X and Y nonzero, the archival kernel omits node 1's `X` term when it eliminates `Y`. Its backward error is then 2.9e-3, which the test suite records. The library handles this case (backward error ≤ 5e-16).

**Second deliberate difference:** a pivot with |p| ≤ n·ε·max|block| is reported as singular. The archival kernel only stops on an exact zero, which misses rank-deficient blocks after rounding.

## 2. Newton driver

Test problem: -(D(c)c′)′ = f with D = 1 + c², c* = sin πx, and 101 nodes.

| Check | Result |
|---|---|
| Linear problem (D = 1) | Converges at iteration 2; the first solve reduces the residual from 9.9 to 2.6e-12. |
| Nonlinear step norms | 1.3, 2.8e-1, 4.8e-2, 1.2e-3, 6.9e-7, 2.4e-13: **quadratic**, 6 iterations. |
| C++ vs Fortran drivers on the same callback | Identical iterates and step histories (max diff 0). |
| Damping 0.5 | Converges in 35 iterations to within 1.4e-11 of the undamped solution. |
| One-step archival mode (`max_iter=1`, `require_convergence=false`) | Returns OK after one correction. |
| Error paths | Callback errors (Fortran `ierr`, C return code, C++ exception, Python exception) propagate. Not-converged and invalid options are reported. |

## 3. Worked examples (`examples/`, also run in `tests/python/test_examples.py`)

Observed order = log₂(error(h)/error(h/2)), using the max-norm error against the exact solution. Both backends gave identical tables.

**Nonlinear steady diffusion** (`nonlinear_bvp.py`). The expected spatial order is 2.

| nj | Newton iterations | max error | order |
|---|---|---|---|
| 21 | 6 | 2.750e-03 | – |
| 41 | 6 | 6.859e-04 | 2.003 |
| 81 | 6 | 1.714e-04 | 2.001 |
| 161 | 6 | 4.284e-05 | 2.000 |
| 321 | 6 | 1.071e-05 | 2.000 |

**Coupled 3-unknown differential-algebraic system** (`reaction_diffusion.py`). It has two diffusion equations, a nonlinear exchange term, and an algebraic row at every node. The Jacobian blocks are full 3×3.

| nj | Newton iterations | max error | order |
|---|---|---|---|
| 21 | 5 | 2.552e-03 | – |
| 41 | 5 | 6.373e-04 | 2.002 |
| 81 | 5 | 1.593e-04 | 2.000 |
| 161 | 5 | 3.982e-05 | 2.000 |
| 321 | 5 | 9.954e-06 | 2.000 |

Newton steps at nj = 81 were 1.0, 2.3e-1, 1.4e-3, 6.6e-9, 5.4e-15, which is quadratic.

**Transient diffusion with implicit Euler** (`transient_diffusion.py`). Setup: nj = 801, t = 0.1, and one Newton solve per step. The expected temporal order is 1.

| Δt | max error | order |
|---|---|---|
| 0.01 | 1.744e-02 | – |
| 0.005 | 8.893e-03 | 0.971 |
| 0.0025 | 4.492e-03 | 0.985 |
| 0.00125 | 2.258e-03 | 0.992 |

**X/Y endpoint blocks** (`xy_boundary.py`). The problem is -c″ + c = (π²+1)cos πx with c′(0) = c′(1) = 0. The second-order one-sided flux boundary needs the X and Y blocks.

| nj | error, 2nd-order closure (X/Y) | order | error, 1st-order closure (no X/Y) | order |
|---|---|---|---|---|
| 21 | 7.024e-04 | – | 1.134e-01 | – |
| 41 | 3.063e-04 | 1.20 | 5.687e-02 | 1.00 |
| 81 | 9.568e-05 | 1.68 | 2.847e-02 | 1.00 |
| 161 | 2.649e-05 | 1.85 | 1.424e-02 | 1.00 |
| 321 | 6.953e-06 | 1.93 | 7.124e-03 | 1.00 |
| 641 | 1.781e-06 | 1.97 | | |
| 1281 | 4.505e-07 | 1.98 | | |
| 2561 | 1.133e-07 | 1.99 | | |

The X/Y closure is pre-asymptotic on coarse grids and approaches order 2, whereas the first-order closure limits the whole solution to order 1.

## 4. Robustness builds and a rough timing check

**Hardened builds, all passing:**
- Fortran Debug build with `-fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=snan`: all 8 CTest suites pass.
- C++ tests and the C++ quickstart under AddressSanitizer and UndefinedBehaviorSanitizer (clang, `-O1`), linked against that checked Fortran library: no reports.

Across these differently optimized builds, the C++ and Fortran results differ by up to 1.6 ε rather than being bit-identical. This is expected, and it is why the tests use tolerances.

**Rough timing**, from Python on Apple Silicon (one core, Release build). These are indicative numbers, not a benchmark.

| n | nj | cpp | fortran |
|---|---|---|---|
| 1 | 100,000 | 4.0 ms | 8.7 ms |
| 5 | 20,000 | 11.6 ms | 13.8 ms |
| 10 | 2,000 | 3.9 ms | 6.5 ms |
| 30 | 2,000 | 54 ms | 85 ms |

The Fortran backend is probably slower because its C interface transposes the blocks and copies them on entry. That was not measured separately. Time scales linearly in nj and roughly as n³.

## 5. Not covered

- Other compilers and platforms (Linux, Intel, Windows) and other optimization levels. The bit-identical agreements above are specific to this toolchain; the tests use tolerances rather than exact equality.
- The Python extension links Homebrew's dynamic libgfortran, so a portable wheel has not been built or tested.
- No performance benchmark or comparison against LAPACK band, KLU, SUNDIALS, etc. That is Stage 7 of the plan.
- Very ill-conditioned or near-singular systems beyond the threshold tests: `min_rel_pivot` is reported, but no conditioning estimate is made.
- The archival battery models are not validated here (Track 2).
