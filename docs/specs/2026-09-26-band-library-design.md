# bandsolver — standalone Newman BAND library: design spec

**Date:** 2026-09-26 · **Status:** draft, approved to proceed under the Track 1 loop; open for user correction.
**Track:** 1 (generic solver library). Track 2 (battery-model audit, `../IMPLEMENTATION_PLAN.md` Stages 0–4) continues separately and does **not** gate this library.

## 1. Intent (what the user said vs. assumptions)

**Said:** port the BAND solvers used in the battery Fortran models (Newman, *Electrochemical Systems*, Appendix C) into a useful standalone library; publish in Fortran **and** C++, each accessible from Python; functionalize the PDE solver so other libraries can use it. Decision recorded 2026-09-26: split into Track 1 (this library) and Track 2 (model audit) — option A.

**Assumed (override any):**
- API level 2: linear block kernel + Newton driver with a user fill callback. Time stepping lives in examples, not the API (recommended default; user did not answer the level question before starting the loop).
- Correctness is established mathematically (residual/backward error, independent dense solves, manufactured solutions, convergence rates) and by comparison with a frozen copy of the archived kernel used as an oracle. No battery model is ported into the library.
- No external publication, push, package upload, or license choice without explicit user authorization (plan Stage 8).

**Success criteria**
1. Fortran and C++ implementations solve the Appendix C system (block tridiagonal plus endpoint `X`/`Y` blocks) with scaled backward error ≲ 1e‑13 on well-conditioned random systems, and agree with each other and with an independent dense solve.
2. `pivot="legacy"` reproduces the frozen archival `BAND`/`MATINV` result to within a few ulps on identical inputs; `pivot="partial"` (default) is documented as a deliberate, characterized change.
3. Singular/non-finite blocks return an error status with the node index — never print-and-continue.
4. Newton driver shows the expected convergence (one iteration for linear problems, quadratic near solution for nonlinear ones) and the expected spatial order on manufactured solutions.
5. Python users call either backend with numpy arrays through one package.

## 2. The mathematical contract (verified against archive source)

For `nj` nodes (`nj ≥ 3`) and `n` unknowns per node, solve for Δc ∈ ℝ^{n×nj}:

```
j = 1:        B_1 Δc_1 + D_1 Δc_2 + X Δc_3                 = G_1
1 < j < nj:   A_j Δc_{j-1} + B_j Δc_j + D_j Δc_{j+1}        = G_j
j = nj:       Y Δc_{nj-2} + A_nj Δc_{nj-1} + B_nj Δc_nj     = G_nj
```

This convention was derived from `BAND(J)` in `Research/Modelling/MnO2/Old Work/ZnMn02_v8.f95` lines 2041–2100 (file SHA-256 `8abba384…c191a`, matches `BAND_SOURCE_INVENTORY.md`) and matches Appendix C eqs C.10–C.15 per audit §8.18. `A_1`, `D_nj` are unused. For Newton use, `G = −F(c)` (negative residual) and the update is `c ← c + λ Δc`.

Algorithm (forward elimination, back substitution): at each node solve the small dense block system with multiple right-hand sides, storing `E_j` (n×n) and `e_j` (n) such that `Δc_j = E_j Δc_{j+1} + e_j` (plus the `X` correction at node 1). Cost O(nj·n³), storage O(nj·n²).

## 3. Architecture

```
bandsolver/
  fortran/src/     band_kernel.f90  (module bandsolver_kernel: block solve, legacy + partial pivot)
                   band_newton.f90  (module bandsolver_newton: Newton driver, fill callback)
                   band_capi.f90    (bind(C) API for C/C++/Python; C row-major layout)
  cpp/include/bandsolver/  band.hpp, newton.hpp, status.hpp
  cpp/src/         band.cpp, newton.cpp
  python/bandsolver/  __init__.py (solve, newton, backends), _core (pybind11 ext linking both cores)
  (../bandsolver-legacy-oracle/, private, outside the repo) frozen archival kernel + oracle wrapper; tests only
  tests/fortran, tests/cpp, tests/python
  examples/        transient diffusion (implicit Euler via newton), nonlinear BVP, X/Y boundary closure
  docs/            math formulation, API reference, validation report, provenance
```

Units and boundaries:
- **Kernel** (both languages): `solve(n, nj, A, B, D, G, X, Y, pivot) → Δc, status, diagnostics`. Pure, no global state; inputs are not modified (workspace is internal or caller-supplied). Diagnostics: failing node, minimum relative pivot seen.
- **Newton driver**: `newton(fill, c0, options) → c, result`. `fill(c) → A, B, D, G, X, Y` for the whole system (vectorizes naturally in numpy; node-by-node fill remains trivial for callers to write). Options: `rtol`, `atol`, `max_iter` (setting `max_iter=1` and `require_convergence=false` reproduces the archival one-correction-per-step usage), fixed damping `λ ∈ (0,1]`, `pivot`. Convergence: `max_ij |Δc_ij|/(atol + rtol|c_ij|) ≤ 1`. Result: converged flag, iterations, per-iteration update and residual (‖G‖∞) norms, status.
- **Layouts**: Fortran-native `A(n,n,nj)`, `G(n,nj)`. C++ and the C ABI use C row-major `A[nj][n][n]`, `G[nj][n]` (numpy shape `(nj,n,n)`/`(nj,n)`). The Fortran `bind(C)` layer transposes blocks on entry/exit (O(nj·n²), negligible vs O(nj·n³)).
- **Python**: one pybind11 extension `bandsolver._core` linking the C++ core and the Fortran library's C ABI; `bandsolver.solve(..., backend="cpp"|"fortran")`, `bandsolver.newton(fill, c0, ..., backend=...)`. Python callbacks cost is documented; native callbacks remain available from C++/Fortran.

## 4. Pivoting and errors

- `partial` (default): Gauss–Jordan/LU with row partial pivoting on each block.
- `legacy`: exact reproduction of the Appendix C / archival `MATINV` pivot heuristic (the `BMAX=1.1` ratio search), for oracle comparison and historical fidelity.
- Status codes (shared C enum): `0 OK`, `1 SINGULAR_BLOCK` (+node), `2 INVALID_ARGUMENT` (n<1, nj<3, bad pivot), `3 NOT_CONVERGED`, `4 NON_FINITE`, `5 CALLBACK_ERROR`. C++ offers exceptions on top of status (`bandsolver::Error`); Python raises `bandsolver.BandError` subclasses.
- Singularity: a block pivot with |p| ≤ n·ε·max|block| is reported as `SINGULAR_BLOCK` in both modes (the legacy kernel only stops on an exact zero). `min_rel_pivot` is reported as a conditioning diagnostic.
- Deliberate differences from legacy: no `PRINT` on singularity; status instead of continuing with garbage; inputs not overwritten; `nj ≥ 3` enforced; correct result for `nj = 3` with both `X` and `Y` nonzero (legacy drops node 1's `X` term there — see `docs/math.md`).

## 5. Testing and validation

TDD per component. Every test runs from `bandsolver/build` or pytest temp dirs — never in archive folders.
1. **Backward error:** random block systems (n ∈ {1,2,5,12}, nj ∈ {3,4,50,500}, with/without X/Y): ‖KΔc−G‖∞ / (‖K‖∞‖Δc‖∞ + ‖G‖∞) ≤ 1e‑13, K applied directly (independent of the solver).
2. **Legacy oracle:** frozen `ABDGXY`-free copy of archival `BAND`/`MATINV` wrapped in a test module (module globals preserved); `legacy` mode matches it to ≤ 4 ulp-scaled tolerance; `partial` matches within backward-error bound.
3. **Cross-backend and dense:** Fortran vs C++ vs numpy dense `solve` of the assembled K.
4. **Failure paths:** exact-singular block at first/interior/last node → `SINGULAR_BLOCK` with correct node; NaN input → `NON_FINITE`; zero pivot needing off-diagonal choice solves under both pivot modes.
5. **Newton:** linear diffusion converges in 1 iteration; nonlinear diffusion (D(c)=1+c²) with manufactured solution shows quadratic convergence and 2nd-order spatial accuracy; coupled 2-species reaction–diffusion; X/Y second-order boundary closure example.
6. **Python:** pytest for both backends, shapes/dtypes/contiguity validation, exceptions, callback errors propagating.

## 6. Build, tooling, provenance

- CMake ≥ 3.20 (languages C, CXX, Fortran), C++17, Fortran 2008; CTest for Fortran/C++ tests; scikit-build-core + pybind11 for the Python package; project venv at `bandsolver/.venv` (numpy, pytest, pybind11, scikit-build-core).
- Local git repo in `bandsolver/`; small commits per verified component; no remote/push without authorization.
- **Update 2026-09-26:** License is BSD-3-Clause, and the legacy extract is kept out of the repo and every release (`docs/provenance.md`). Originally: `legacy/` holds a byte-exact extract with source path, line range, and hashes. It is a test oracle only. Before any public release the user must decide licensing/attribution for the Appendix C–derived algorithm and whether `legacy/` ships (Stage 8 gate). No LICENSE file is chosen by the agent.
- Archive sources are read-only; the oracle extract is copied, never edited in place.

## 7. Out of scope (YAGNI)

Adaptive time stepping / DAE integration (examples only; SUNDIALS comparison is Stage 7), line search beyond fixed damping, sparse/general bandwidth beyond X/Y, AD, multithreading, battery physics, PyPI/conda packaging.
