# Track 1 loop — standalone BAND library

Self-paced `/loop`. Each iteration: read this file, take the **first unchecked task whose dependencies are done**, implement it test-first, run its verification, commit in `bandsolver/` (local only), tick the box, append a checkpoint, schedule the next wakeup. Stop the loop at a user gate or when all tasks are done.

Spec: [`docs/specs/2026-09-26-band-library-design.md`](docs/specs/2026-09-26-band-library-design.md). Project context: `../HANDOFF_PROGRESS.md`, `../AGENT_INSTRUCTIONS.md`.

## Rules
- Never edit, move, build, or run anything in archive folders (`~/Desktop/Research`, MnO2*, Blackbox Testing, etc.). Copy only; verify SHA-256.
- Build only in `bandsolver/build*`; Python env only `bandsolver/.venv`.
- No push, publish, upload, or license selection without the user's explicit OK.
- Report failures honestly in the checkpoint; do not tick a task whose verification did not pass.

## Tasks
- [x] T0 Spec, loop file, plan/handoff track-split note.
- [x] T1 Scaffold: `git init`, CMake skeleton (C/CXX/Fortran, CTest), `.venv` with numpy/pytest/pybind11/scikit-build-core, `legacy/` byte-exact extract of BAND/MATINV (v8 lines 1983–2100) + PROVENANCE.md with hashes + oracle wrapper module that compiles. Verify: cmake configure+build, oracle solves a 3-node system.
- [x] T2 Fortran kernel (`band_kernel.f90`): partial + legacy pivot, status codes. Tests: backward error sweep, legacy-oracle match, singular/NaN paths.
- [x] T3 Fortran Newton driver + `bind(C)` API (row-major). Tests: linear 1-iteration, nonlinear quadratic convergence, C-ABI round trip from a C test.
- [x] T4 C++ kernel (`band.hpp/.cpp`). Tests: backward error, match Fortran C-ABI, failure paths.
- [x] T5 C++ Newton driver. Tests mirror T3; cross-check vs Fortran.
- [x] T6 Python package (pybind11 `_core`, both backends). pytest: dense numpy comparison, both backends, errors, callbacks.
- [x] T7 Validation examples + report: manufactured nonlinear BVP (order 2), transient diffusion w/ implicit Euler, 2-species reaction–diffusion, X/Y closure. `docs/validation.md` with actual numbers.
- [x] T8 Docs: README, math formulation (Appendix C), API reference (Fortran/C++/C/Python), build instructions, provenance/licensing open questions.
- [x] T9 Final code review pass, fix findings, full test run, summary to user; stop loop.

## Checkpoints
- **2026-09-26 — T0.** User chose option A (Track 1 library decoupled from Track 2 audit gates) and asked for a loop. API level 2 assumed (user did not pick). Spec written; oracle source `ZnMn02_v8.f95` hash re-verified against inventory (`8abba384…c191a`). Next: T1.
- **2026-09-26 — T1.** Ran `git init` (local only). Added the CMake skeleton (C/CXX/Fortran, CTest) and `.venv` (numpy 2.5.3, pybind11 3.1.0, pytest 9.1.1, scikit-build-core 1.1.0). Made a byte-exact legacy extract (SHA `99eb0fb8…4ef8f`), reproduced both normalized inventory routine hashes, and wrote the oracle wrapper. `test_legacy_oracle` passes; its backward error is <1e-13 on random systems. **Finding:** the legacy kernel returns a wrong solution for nj=3 when both X and Y are nonzero (backward error 2.9e-3). It drops node 1's X term when substituting Y. This is characterized in the test and `legacy/PROVENANCE.md`; the library kernel will handle the case. Toolchain: gfortran 14.2.0, Apple clang 16, CMake 3.31.2, Python 3.13.3. Next: T2.
- **2026-09-26 — T2.** Added the `bandsolver_kernel` module with `band_solve` and `block_solve`, in partial and legacy pivot modes, with status codes. The test sweep covered n∈{1,2,5,12}, nj∈{3,4,50,500}, with and without X/Y, in both modes. Max backward error was 4.0e-16 (partial) and 3.9e-16 (legacy), including nj=3 with X and Y. **Legacy mode matches the frozen oracle exactly (0 ulp)** in a Release (-O3) gfortran 14.2 build. Partial mode differs from the oracle by at most 9.6e-16 relative. The singular-block tests, which report the failing node, pass at the first, interior and last nodes in both modes. Added a relative singular-pivot threshold, |p| ≤ n·ε·max|block|; the spec now documents it. The first test run showed that exact-zero detection misses rank-deficient blocks, which is the same weakness the legacy kernel has. The off-diagonal-pivot, NaN and invalid-argument tests pass, and inputs are unchanged after a solve. Next: T3.
- **2026-09-26 — T3.** Added the `bandsolver_newton` module. It uses an abstract `band_problem` type with a deferred `fill`, which is reentrant and uses no module globals. It has options and result types, including update, step and residual histories. Added `bandsolver_capi` with `bind(C)`, a row-major layout and a C callback with a `void*` context, plus the header `fortran/include/bandsolver_f.h`. Test results:
  - Linear BVP: converges at iteration 2, and the residual drops from 9.9 to 2.6e-12 after the first solve.
  - Nonlinear BVP with D=1+c², 101 nodes: converges in 6 iterations with quadratic step decay (4.8e-2 → 1.2e-3 → 6.9e-7 → 2.4e-13). Max error against sin(πx) is 1.1e-4.
  - The legacy-pivot Newton agrees with partial pivoting.
  - One-step archival mode returns OK; NOT_CONVERGED, callback-error and invalid-damping paths work.
  - C test: nonsymmetric row-major blocks with X/Y reach 1.4e-16 backward error; NULL X/Y, singular node 6 (1-based) and invalid nj are handled; C-callback Newton converges and callback errors propagate.

  All 4 CTest suites pass. Development note: gfortran rejects a procedure-pointer component for C function pointers, so the pointer is stored as `c_funptr` and converted at call time. Next: T4.
- **2026-09-26 — T4.** Added the C++17 core (`cpp/include/bandsolver/band.hpp`, `cpp/src/band.cpp`). It provides `SystemView` (non-owning, row-major, nullable X/Y), the owning `BlockSystem`, a non-throwing `solve()` that returns `SolveInfo` with a 0-based failing node, a throwing `solve(BlockSystem)` that raises `bandsolver::Error`, and the `Status`/`Pivot` enums. It builds warning-free under -Wall -Wextra -Wpedantic -Wconversion. Test results:
  - Sweep over the same n/nj set, both pivots, with X/Y: max backward error 5.0e-16.
  - **The C++ result equals the Fortran library result bit for bit (0 eps)** in both pivot modes (Apple clang 16 vs gfortran 14.2, Release). The test allows 64 eps for other toolchains.
  - Singular nodes 0, 3 and 6 are reported through both the status and the exception APIs; zero-diagonal pivots, NaN, nj<3, an unknown pivot, a null output and null X/Y are handled.

  All 5 CTest suites pass. Next: T5.
- **2026-09-26 — T5.** Added the C++ Newton driver (`newton.hpp/.cpp`): `NewtonOptions`, a `NewtonResult` with histories, and a `FillFunction` taking `(const double* c, BlockSystem&)`. An exception thrown by the fill function becomes `callback_error`, and `exception_ptr` keeps it so the Python binding can re-raise it. Test results:
  - It mirrors T3: the linear problem converges at iteration 2, and the nonlinear one converges quadratically in 6 iterations with the same step norms as Fortran and error 1.1e-4.
  - **The C++ and Fortran drivers give identical solutions and step histories** on the same C callback (max diff 0).
  - Legacy pivoting agrees; one-step mode, not_converged, the rethrowable callback exception, invalid damping and an empty fill are all handled.
  - Damped Newton (λ=0.5) converges in 35 iterations to within 1.4e-11 of the undamped solution.

  All 6 CTest suites pass. Next: T6 (Python).
- **2026-09-26 — T6.** Added the Python package `bandsolver`, built with scikit-build-core and pybind11 (`pyproject.toml`; `pip install . --no-build-isolation` in `.venv`). One extension, `_core`, links the C++ core and the Fortran C ABI.
  - API: `solve(A,B,D,G,X,Y,pivot,backend)` and `newton(fill,c0,...)` returning a `NewtonResult`, with exceptions `SingularBlockError` (0-based `.node`), `NonFiniteError`, and `NotConvergedError` (carries `.result`). Exceptions raised inside `fill` propagate unchanged on both backends; the Fortran path goes through a trampoline that catches the exception before it can unwind through Fortran frames.
  - Linear solves release the GIL.
  - pytest: **72 passed**. Results agree with the numpy dense solve (the dense matrix is assembled independently) to 1e-12 relative, with backward error <1e-14, across the n/nj/X-Y/pivot grid and both backends. The two backends agree with each other. Inputs are unchanged, and Fortran-ordered or strided inputs are accepted. Singular nodes, NaN input and argument validation are handled. The Newton tests show quadratic convergence on both backends with identical iterates; the linear problem with X/Y matches `solve`; callback exceptions, bad fill output, not-converged, one-step mode and a singular Jacobian are all handled.
  - Not yet checked: the wheel links the Homebrew libgfortran dynamically, so portability to machines without gfortran is untested (document in T8).

  Next: T7.
- **2026-09-26 — T7.** Added four examples: `nonlinear_bvp`, `reaction_diffusion` (n=3 with an algebraic row), `transient_diffusion` (implicit Euler through `newton`), and `xy_boundary`. They also run as `tests/python/test_examples.py`; pytest now has **82 passed**. Observed orders:
  - nonlinear BVP: 2.000
  - coupled DAE: 2.000, with quadratic Newton convergence
  - implicit Euler in time: 0.97 → 0.99
  - X/Y second-order Neumann closure: 1.97 → 1.99 as the grid refines, against 1.00 for the first-order closure.

  Both backends gave identical tables. Wrote `docs/validation.md` with all measured numbers and a "Not covered" section. Next: T8.
- **2026-09-26 — T8.** Docs: `README.md` (overview, build/install, three quickstarts, status and limitations, license not yet chosen), `docs/math.md` (linearization, BAND recurrence including the X/Y and nj=3 terms, pivoting, a table of differences from the archive), `docs/api.md` (Fortran/C/C++/Python conventions, status codes, fill contract, performance notes), and `docs/provenance.md` (origin, legacy oracle, five open decisions for the author: license, shipping legacy/, citation, package name, publication venue).
  - The Fortran and C++ quickstarts are compiled programs registered as CTest cases, and the README's Python snippet runs as a pytest case, so the documented code is executed.
  - Fixed the CMake `enable_testing()` order, which had kept the example tests from registering.
  - Checks: ctest 8/8 pass; pytest 83 passed.
  - The book edition and the Newman 1968 citation are marked for the author to verify.

  Next: T9.
- **2026-09-26 — T9 (final).**
  - Robustness builds: the Fortran `-fcheck=all` build with FPE traps and signalling-NaN initialization passes 8/8 CTest suites. The C++ tests and quickstart under ASan and UBSan are clean.
  - Reviewed the kernel, the Newton drivers, the C ABI and the bindings for index, lifetime, GIL and exception-unwinding issues; no defects found beyond those already fixed in T1–T8.
  - A rough timing check is recorded in `docs/validation.md` §4.
  - Final state: ctest 8/8 pass; pytest 83 passed.

  **Track 1 loop complete; loop stopped.** Open items for the user are the decisions in `docs/provenance.md`. Possible follow-ups: portable wheels (static libgfortran or delocate), Linux CI, a zero-copy Fortran path, and Stage 7 benchmarks.
- **2026-09-26 — Release decisions (user).** License: **BSD-3-Clause** (`LICENSE`, pyproject `license`). The legacy extract has been **removed from the repo and every release**. It now lives privately, hash-verified, in `../bandsolver-legacy-oracle/`, and CMake builds the two legacy comparison tests only when that folder is present. The oracle-free `test_kernel` was split from the new `tests/fortran/legacy/test_kernel_legacy.f90`.
  - Verified: with the oracle, ctest 9/9 pass and legacy mode is still 0 ulp; without it, 7/7 pass.
  - sdist and wheel built via `python -m build`, with the wheel built from the sdist. Neither contains legacy, oracle or log files. The wheel carries `License-Expression: BSD-3-Clause`, and pytest passes 83/83 against the installed wheel.
  - The stray tracked `build.log` was removed.
  - **Git history still contains `legacy/`** in earlier commits, so history must be rewritten or squashed before any push.
- **2026-09-26 — Fresh repository (user chose option 2).** The development history (12 commits, including `legacy/`) was moved to a private backup at `../bandsolver-dev-history.git`, which can be inspected with `git --git-dir=../bandsolver-dev-history.git log`. `bandsolver/` was re-initialized on branch `main` with one commit of the current tree: 45 files, and no `legacy/` anywhere in its history. A clean clone configures without the oracle (legacy tests skipped), passes ctest 7/7, installs with isolated `pip install .`, and passes pytest 83/83. Still not pushed or published.

## Feature loop: finite-difference Jacobians (branch `fd-jacobian`)
Spec: [`docs/specs/2026-09-27-fd-jacobian-design.md`](docs/specs/2026-09-27-fd-jacobian-design.md). Same rules as above. Work on the branch and merge via a PR once CI is green.
- [x] F1 C++ core: `fd_jacobian`, `fd_fill`, `newton_fd`, `check_jacobian` + tests (analytic match incl. nj=3 X/Y, eval counts, planted-error detection, convergence).
- [x] F2 Fortran core + C ABI: `band_residual_problem`, `band_fd_jacobian`, `band_newton_fd`, `band_check_jacobian`, `bandsolver_f_fd_jacobian`, `bandsolver_f_newton_fd` + tests; cross-check against C++.
- [x] F3 Python: `fd_jacobian`, `newton_fd`, `check_jacobian` on both backends + pytest; FD example; docs (api, math, README, validation).
- [x] F4 Open a PR, wait for green CI and wheels on all platforms, merge; update the handoff.

### Checkpoints
- **2026-09-27 — F0.** The author approved the design; spec written; branch `fd-jacobian` created. Next: F1.
- **2026-09-27 — F1.** Added `cpp/include/bandsolver/fd.hpp` and `cpp/src/fd.cpp`, containing `fd_jacobian` (period-3 colouring, 3n+1 evaluations), `fd_fill`, `newton_fd` and `check_jacobian`. `NewtonResult` gains `residual_evaluations`. Test results:
  - On a nonlinear n ∈ {1,3} problem with nonlinear X/Y terms, the finite-difference blocks match the analytic ones to ≤ 1.6e-8 relative for nj ∈ {3,4,5,10}, always with exactly 3n+1 evaluations.
  - `newton_fd` takes the same 9 iterations as analytic Newton and reaches the same solution within 3e-15. Its evaluation count is exactly (3n+1)·iterations.
  - **Design adjustment:** `check_jacobian` now measures each entry against its equation row's scale, so the floor is 1e-3·rowscale. The first test run showed that per-entry relative error flags ordinary finite-difference noise on tiny entries (3.6e-4 on a correct Jacobian). After the change, a correct fill scores 1.4e-5, a planted D error 1.6 (located at the right node, row and column), and a missing X entry about 1.
  - Exceptions propagate, and invalid options are rejected.
  - ctest passes 10/10, and the new tests are clean under ASan/UBSan.

  Next: F2.
- **2026-09-27 — F2.** Added the Fortran module `bandsolver_fd` (`fortran/src/band_fd.f90`): the abstract `band_residual_problem` type, `fd_options`, `band_fd_jacobian`, `band_newton_fd` (through an internal fill adapter) and `band_check_jacobian` (1-based locations, row-scaled metric). `newton_result` gains `residual_evaluations`. The C ABI gains `bandsolver_f_default_fd_options`, `bandsolver_f_fd_jacobian`, `bandsolver_f_newton_fd` and `bandsolver_f_check_jacobian`, with the new structs in `bandsolver_f.h`. Test results:
  - `test_fd` (Fortran): the finite-difference blocks match analytic for n ∈ {1,3}, nj ∈ {3,4,5,10} with 3n+1 evaluations; `band_newton_fd` matches analytic Newton to 4e-16 with (3n+1) evaluations per iteration; a correct fill scores 1.3e-6; the planted error is located at node 5, row 2, column 3 (1-based); the missing X entry is detected; errors and invalid options are handled.
  - `test_fd_cross`: **the C++ and Fortran results are bit-identical**, covering the Jacobian blocks, the Newton iterates and evaluation counts (9 iterations, 90 evaluations each), and the check report (same location after the 1-based offset, same score).
  - Shared test problem moved to `tests/cpp/fd_problem.hpp`.
  - ctest passes 12/12 in both Release and the Fortran `-fcheck=all` + FPE-trap build.

  Next: F3.
- **2026-09-27 — F3.** Python gains `fd_jacobian`, `newton_fd` (whose `NewtonResult.residual_evaluations` counts residual calls) and `check_jacobian`, which returns a `JacobianCheck` of `JacobianMismatch` values with `.max_error` and `.worst()`. All work on both backends, and residual and fill exceptions propagate unchanged, including through the Fortran path via a new trampoline.
  - Added `examples/fd_jacobian.py`: on the coupled n=3 DAE, `newton_fd` takes the same 5 iterations as analytic Newton, uses 50 residual calls independent of nj, matches within 3e-14, and has order 2.000. `check_jacobian` scores the correct fill 3.6e-9 and locates a planted sign bug at block B, row 1, column 0 with score 0.16. The score is modest because the entry is tiny next to the row's 2/h² scale, but it is well above the 1e-3 bug line.
  - Docs updated: math.md (new finite-difference section), api.md (all four interfaces), README (residual-only snippet), validation.md §5. `test_readme` now runs every README Python block.
  - Tests: pytest 103 passed; ctest 12/12.

  Next: F4.
- **2026-09-27 — F4 (in progress).** PR #2 opened. Every check passed except the Windows wheel job. There, `test_backends_identical` demanded bit-for-bit equality between the C++ core (MSVC) and the Fortran core (ifx); the solutions differed by ≤ 4.4e-16 (1–2 ulp) in 21 of 60 values. This is a test flaw, not a solver defect. The test is now `test_backends_agree`, with tolerances of 1e-13 on solutions and 1e-6 on the FD blocks (forward differences amplify an ulp in F to ~1e-8). validation.md now scopes the bit-identical claim to the reference toolchain. Pushed for re-run.
- **2026-09-27 — F4.** PR #2's re-run passed all 11 checks: CI on 5 platforms, wheels on 5 platforms, and the sdist. Merged into `main`. **Feature loop complete.**

## Benchmark loop: solver comparison, layers 1–2 (branch `benchmarks`)
Spec: [`docs/specs/2026-09-27-benchmarks-design.md`](docs/specs/2026-09-27-benchmarks-design.md). Work on the branch and merge via a PR. Benchmarks run locally; CI only smoke-tests the scripts.
- [x] B1 Harness: `benchmarks/` layout, pinned `requirements.txt`, `env.py`, native C++ timing harness + CMake option; smoke test.
- [x] B2 Layer 1: `linear.py` (BAND py/native, LAPACK band, SuperLU, dense), backward-error check, sweep → `results/linear.csv`, plots.
- [x] B3 Layer 2 models: bandsolver BE/BDF2, IDA DAE (band), SciPy BDF reduced ODE; verify all agree on the same discrete solution.
- [x] B4 Layer 2 sweeps: work-precision + mesh scaling → `results/transient*.csv`, plots.
- [x] B5 `docs/benchmarks.md`, notebook 05 (reads the CSVs), README link, CI smoke test; PR, CI, merge.

### Checkpoints
- **2026-09-27 — B0.** Author approved layers 1–2 (PyBaMM deferred). scipy 1.18.1 and scikit-sundae 1.1.3 installed in `.venv`. The IDA API (band solver, algebraic_idx, calc_initcond, nfev/njev) was checked. Spec written. Next: B1.
- **2026-09-27 — B1.** Added `benchmarks/`: `requirements.txt` (pins numpy 2.5.3, scipy 1.18.1, scikit-sundae 1.1.3, matplotlib 3.11.2), `env.py` (platform, CPU, compilers, versions, git commit to `results/env.json`), and `native_bench.cpp`. The native harness times the C++ and Fortran cores with no Python, taking the median of ≥5 repeats and ≥0.2 s per point; its CMake option `BANDSOLVER_BUILD_BENCHMARKS` is OFF by default, and a `--quick` smoke test is registered with CTest. Smoke run on an Apple M1 Pro: n=1, nj=25 takes 1.1 µs (C++) and 2.3 µs (Fortran); n=3, nj=50 takes 8.7 and 13.2 µs; backward error ≤ 1.5e-16. `benchmarks/results/` is excluded from the sdist. Next: B2.
- **2026-09-27 — B2.** `benchmarks/linear.py` (+ `plot.py`) ran the full sweep: n ∈ {1,3,5,10,20,30} × nj ∈ {25…2000} on an M1 Pro, 322 measurements in 80 s, all with backward error ≤ 8.7e-16. Results: `results/linear.csv`, `linear_time.png`, `linear_speedup.png`.
  - **BAND C++ (from Python) is fastest at every size.** It beats LAPACK banded (dgbsv, bandwidth widened to 3n−1 by X/Y) by 1.6–5.8× (median 2.4×) and SciPy SuperLU by 2.4–12.4× (median 4.8× natural ordering, 6.7× COLAMD). Dense only wins nowhere past about N = 50.
  - Factor storage at n=30, nj=2000: BAND 14.9 MB, LAPACK band 128.6 MB, SuperLU 65.5 MB.
  - Python-call overhead is about 5 µs fixed: 4.9× on n=1, nj=25, but ≤ 1.18× from nj ≈ 200.
  - Scaling is linear in nj (×2.05 per doubling); n 10→30 costs ×13.8 (below n³ = 27, a small-block efficiency effect).
  - **Finding:** the Fortran core is 1.5–2.2× slower than C++ *natively*, not just because of the C-ABI transposes. Its loops follow the legacy row-major order, which is cache-unfriendly in column-major Fortran. **Follow-up TODO:** optimise the partial-pivot path's loop order (legacy mode must keep its operation order).

  Next: B3.
- **2026-09-27 — B3.** Added `benchmarks/transient.py`. One `Model` class holds the notebook-2 finite-volume discretization (nj=81, t_end=5 s), and there are three stacks:
  - `run_bandsolver`: BE or BDF2 with fixed Δt, analytic blocks;
  - `run_ida`: IDA via scikit-sundae, band solver with lband=uband=3, `algebraic_idx`=φ, `calc_initcond='yp0'`;
  - `run_scipy`: BDF on the reduced ODE, φ eliminated exactly per face, tridiagonal `jac_sparsity`.

  **Consistency against an IDA rtol=1e-12 reference:** IDA at 1e-10 differs by 2.3e-8 mol/m³, SciPy at 1e-10 by 6.7e-8, bandsolver BDF2 at Δt=1e-4 by 1.3e-9, and BE at Δt=1e-4 by 1.3e-5 (its O(Δt) error). All stacks share one discrete solution. The DAE and the reduced ODE are equivalent, as designed.

  Early signal: for tight accuracy, fixed-step BE/BDF2 needs 50k steps (7 s), where adaptive IDA/SciPy need 365–710 steps (0.05 s). Next: B4 (work-precision and mesh sweeps).
- **2026-09-27 — B4.** `transient.py` now has the consistency check, work-precision (nj=81) and mesh sweeps (nj 41–1281); results are in `results/transient_wp.csv`, `transient_mesh.csv`, `transient_work_precision.png` and `transient_mesh.png`. The full run takes about 25 s.
  - **Per-step cost:** bandsolver full Newton runs 2.6–2.8 iterations/step; each iteration spends 51 µs in the Python fill and 10 µs in BAND at nj=81 (94/113 µs at nj=1281). IDA averages 1.3 residual evaluations/step and reuses its Jacobian (22 Jacobians in 143 steps).
  - **Added variant, BDF2 with 1 Newton iteration/step** (the archival linearized usage): same accuracy as full Newton on this mildly nonlinear problem, e.g. 7.37e-6 at Δt=0.01, 2.7× cheaper. Its per-step cost is 64–242 µs for nj 41–1281, below IDA and SciPy up to nj ≈ 1000 and tied with IDA at 1281 (242 vs 238 µs).
  - **Work-precision crossover** near an error of 1e-4 mol/m³ (1e-6 relative). For looser targets, linearized BDF2 is fastest: 8.1e-4 in 3.4 ms, vs SciPy 5.9e-4 in 8.3 ms and IDA 3.0e-3 in 7.6 ms. For tighter targets, adaptive order-5 BDF wins: IDA reaches 2.3e-8 in 51 ms, where BDF2 needs 343 ms for 7.3e-8.
  - Backward Euler is dominated everywhere (first order).
  - Conclusion: BAND is the faster kernel per step, and the remaining gap is time-integration strategy (adaptive step/order, Jacobian reuse), not linear algebra. Follow-up idea: adaptive BDF stepping, or use BAND as IDA's linear solver.

  Next: B5.
- **2026-09-27 — B5 (in progress).** Wrote `docs/benchmarks.md` (summary, methods, results tables, where the time goes, conclusions, fairness and limitations, reproduction commands) and `notebooks/05_benchmarks.ipynb` (reads the saved CSVs). Added README "Performance" and tutorial links, a CHANGELOG entry, and a CI benchmark smoke test (native harness via CTest; `linear.py --quick` and `transient.py --quick` on Linux and macOS). Corrected an overstatement before publishing: linearized BDF2 has the cheapest step *up to about 1000 nodes* and ties with IDA at 1281. Local tests: notebooks 5/5 and the README test pass. PR next.
- **2026-09-27 — B5 done.** On PR #7's first CI run, the Linux aarch64 benchmark smoke step failed: scikit-sundae 1.1.3 has no wheel there, and the source build needs SUNDIALS. Fixed so that CI installs scikit-sundae only as a binary; the `--quick` run skips the IDA cases without it, and the full run requires it and exits with a clear message. Both paths were tested locally in an environment without scikit-sundae. The rerun passed on all 5 platforms. Merged. **Benchmark loop complete.**

## Performance loop: fast Fortran kernel, Jacobian reuse, adaptive integrator
Spec: [`docs/specs/2026-09-27-performance-options-design.md`](docs/specs/2026-09-27-performance-options-design.md). There are three PRs (A, B, C), each with CI and a benchmark update. Every option can be switched on and off, and each is benchmarked against the current standard.
- [x] P1 (PR A) Fortran fast kernel + `kernel` option (Fortran, C ABI, Python); bit-identity tests; Layer 1 benchmark reference vs fast; PR, CI, merge.
- [x] P2 (PR B) C++ factor/solve (`Factorization`) + tests.
- [x] P3 (PR B) Fortran factor/solve + C ABI + Python `bs.factor`; tests incl. cross-check.
- [x] P4 (PR B) Newton `jacobian_reuse` option + residual-only callback (C++, Fortran, Python); tests; Layer 2 reuse rows; PR, CI, merge.
- [ ] P5 (PR C) C++ adaptive BDF1–2 DAE integrator with options; tests (orders, tolerance, DAE vs IDA, reuse, mask).
- [ ] P6 (PR C) Python binding `bs.integrate` + tests + example.
- [ ] P7 (PR C) Layer 2 benchmarks for all option combinations; docs/benchmarks.md before/after; notebook 5 update; PR, CI, merge.
- [ ] P8 (follow-up) Fortran port of the integrator.

### Checkpoints
- **2026-09-27 — P0.** The author approved the design. Spec written; branch `fortran-fast-kernel`. Next: P1.
- **2026-09-27 — P1.** Fortran fast kernel with the `kernel` option (Fortran, C ABI `bandsolver_f_solve_kernel` plus a trailing `kernel` field in the newton options, Python `kernel=`).
  - **Bit-identical to the reference loops:** 0 mismatches in 40 sweep cases in Fortran, plus the Python tests.
  - **First attempt:** the fast kernel was 5–22% *slower* at n ≤ 3. Fixes: the fast path falls back to the reference loop order below n = 4, and explicit-loop pivot search replaces maxval/maxloc temporaries (used in both paths; still bit-identical).
  - **Result:** fast is never slower. Reference/fast is 1.00 at n ≤ 3 and 1.1–1.48× for n ≥ 5 (n=20 ~1.0–1.08).
  - **Against C++:** the Fortran core called directly is at parity for n ≥ 5 (0.96–1.14×), versus 1.5–2.2× before. Through the C ABI it's 1.06–1.34×. At n ≤ 3 it's still 1.3–2.2× (a per-node fixed cost; follow-up).
  - Layer 1 re-run: 406 rows, backward error ≤ 8.7e-16. docs/benchmarks, api.md and the CHANGELOG are updated. ctest 13/13 (the `build` dir), pytest 114. **CI finding:** Windows ifx gave 8/40 bitwise mismatches (its default fast FP model vectorises the two loop forms differently). Fix: Intel builds use the precise FP model; the test requires exact equality on gfortran and ≤ 16ε elsewhere, printing the actual difference.
- **2026-09-27 — P1 merged** (PR #8, all 11 checks green after the ifx fix).
- **2026-09-27 — P2.** C++ `Factorization` / `factor()` (`cpp/include/bandsolver/factor.hpp`, `src/factor.cpp`): per-node LU with partial pivoting of B^_j, the E_j blocks, the effective A'_j (Y at the last node, B += Y X' when nj=3), and X'. The `solve` does a forward sweep and back substitution. Tests: factor+solve vs one-shot solve differ by ≤ 6.6e-16 relative over the n/nj/X-Y sweep incl. nj=3; repeated RHS work; singular node reported at factor time; non-finite and invalid inputs rejected. Clean under ASan/UBSan.
  - **Timing** (nj=1000, native): one factorization is 1.0–2.2× the cost of a one-shot solve (n=1: 66 vs 37 µs; n=10: 3.2 vs 1.8 ms; n=30: 59 vs 27 ms), and a re-solve is 2–12× cheaper (19 µs, 0.33 ms, 2.3 ms). Reuse breaks even after about 2–3 uses of one factorization.
- **2026-09-27 — P3.** Fortran `bandsolver_factor` (`band_factorization`, `band_factor`, `band_factor_solve`; column-oriented LU and sweeps), C ABI opaque handle (`bandsolver_f_factor` / `_factor_solve` / `_factor_free`), Python `bs.factor(...)` returning a `Factorization` with `.solve(G)` on both backends (GIL released during the solve).
  - **Tests:** Fortran factor+solve vs band_solve ≤ 5.97e-16 relative (sweep incl. nj=3, X/Y); repeated RHS work; singular node reported at factor time (1-based 6); invalid inputs rejected. Python: 23 tests on both backends (vs `solve`, cross-backend, singular 0-based node, shape and NaN errors).
  - Totals: ctest 15/15, pytest 137.
  - A gfortran `-Wdo-subscript` warning (index j−2 inside the j==nj branch) was removed by indexing with nj−2 explicitly.
- **2026-09-27 — P4.** The `jacobian_reuse` option is now in the C++ and Fortran Newton drivers (refresh when a factorization has been used reuse_max_iter times or step_k > contraction·step_{k−1}), with an optional residual-only callback: C++ `ResidualFunction`; Fortran overridable `band_problem%residual` with a default fallback; C `bandsolver_f_newton_ex`; Python `residual=`. `newton_fd` reuse iterations cost 1 residual evaluation. Results report jacobian_evaluations and factorizations.
  - **Tests:** C++ full Newton 9 it / 9 factorizations vs reuse 12 it / 4 factorizations, same root (7e-15); with the residual callback, fills == factorizations. The Fortran FD + reuse run took 8 it, 2 factorizations and 26 evaluations; the C++ and Fortran reuse drivers are identical through the C ABI (12/4/8, diff 0). Python has 10 new tests; pytest 147, ctest 16.
  - **Benchmark (honest):** within each BDF2 step, reuse saves only 12–18% (per step 140–546 µs vs 160–662 µs full Newton; the 1-iteration variant costs 64–242 µs). Each step needs about 3 Newton iterations to a very tight tolerance, and the Python residual dominates. Bigger gains need cross-step reuse and a matched Newton tolerance, which come with the integrator (PR C). docs/benchmarks, api.md and the CHANGELOG are updated. Next: PR B CI and merge, then P5.
