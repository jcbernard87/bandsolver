# Solver comparison (Layers 1–2): design

**Date:** 2026-09-27 · **Status:** approved by the author ("start with layers 1 and 2"). Layer 3 (PyBaMM) is deferred.

## Goal
Show where BAND is fast and accurate, and where other tools do better. The measurement is layered so that linear algebra, time integration and interpreter overhead are not conflated. Everything must be rerunnable from pinned versions, with raw results saved.

## Layer 1: linear solves
- **Systems:** random, well-conditioned Appendix C block systems, including the X/Y blocks (seeded).
  - nj ∈ {25, 50, 100, 200, 500, 1000, 2000}
  - n ∈ {1, 3, 5, 10, 20, 30}
- **Contenders:**
  - `bs.solve` with the cpp and fortran backends (Python call, including input conversion);
  - **native** BAND C++ and Fortran timings from a small C++ harness, with no Python;
  - LAPACK banded (`scipy.linalg.solve_banded`, dgbsv). The bandwidth is l = u = 3n−1 because of X/Y; the matrix is pre-assembled in band storage;
  - SciPy sparse SuperLU (`splu` factor + solve on CSC; natural and COLAMD orderings);
  - dense `numpy.linalg.solve`, only when nj·n ≤ 3000.
- **Timing:** only the factor+solve step; matrix assembly into each format is excluded and reported separately. Each timing is the median of repeats, with at least 0.2 s per point.
- **Metrics:** time, scaled backward error computed from the blocks, and storage in bytes (analytic for dense, band and BAND; nnz(L+U) for SuperLU).

## Layer 2: nonlinear transient (binary electrolyte, notebook 2)
- **Same finite-volume discretization throughout** (nj = 81 by default, plus a mesh-scaling sweep). The reference is the same discretization integrated with IDA at rtol = 1e-12 and atol ≈ 1e-12 × scale, so the errors measure *time-integration* error only.
- **Contenders:**
  - **bandsolver + implicit Euler** (fixed Δt, user stepping), and **bandsolver + BDF2** (fixed Δt, BE start). Stepping is user code, so the second-order variant is a fair "what a user would write".
  - **SUNDIALS IDA** via scikit-sundae 1.1.3: variable-order (1–5), adaptive BDF on the DAE (c, φ interleaved). It uses a band linear solver with lband = uband = 3, its internal FD Jacobian, `algebraic_idx` set to the φ entries, and consistent initial conditions from `calc_initcond='yp0'`. A variant supplies the analytic Jacobian (via bandsolver's fill blocks converted to band storage) if time permits.
  - **SciPy `solve_ivp(BDF)`**. SciPy has no DAE support, so φ is eliminated exactly at each face using current conservation (i = I at every face), which leaves an ODE in c with a tridiagonal Jacobian sparsity. The discrete solution is identical to the DAE one, which is verified.
- **Sweeps:** a tolerance sweep (IDA and SciPy rtol 1e-2…1e-10) and a Δt sweep for bandsolver give **work-precision** curves (max error at t_end against wall time). They also record steps, residual/RHS evaluations, Jacobian evaluations and Newton iterations. A mesh sweep (nj 41…1281) at fixed accuracy shows scaling.

## Deliverables
- `benchmarks/`: `requirements.txt` (pinned), `env.py` (records the platform, CPU and versions), `linear.py`, `native_bench.cpp` (+ CMake target, off by default), `transient.py`, `plot.py`, and `results/*.csv` + `results/env.json` + figures.
- `docs/benchmarks.md`: method, results, fairness notes and limitations.
- `notebooks/05_benchmarks.ipynb`: a summary that reads the saved CSVs and does not re-run the benchmarks.
- CI: a quick smoke test of the benchmark scripts at tiny sizes, so they don't rot.

## Fairness and limitations (to state in the results)
- There is a single machine (Apple Silicon); timings are indicative.
- Python callbacks are used for all Layer 2 contenders; native Layer 1 timings isolate the interpreter overhead.
- BAND is a direct block solver with no fill-in beyond the block structure; banded LAPACK stores the X/Y-widened band.
- IDA's error control is local, so its tolerance is not the global error, and the work-precision plot uses the measured global error.
