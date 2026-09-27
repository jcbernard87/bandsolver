# Performance options: fast Fortran kernel, Jacobian reuse, adaptive integrator

**Date:** 2026-09-27 · **Status:** approved by the author ("great lets do it"). Every improvement is an **option that can be turned on or off**, and each is benchmarked against the current standard method.

## A. Fast Fortran kernel (PR A)
- **Problem:** measured natively, the Fortran core is 1.5–2.2× slower than C++, because its innermost loops stride along the second (non-contiguous) index.
- **Change:** in the **partial-pivot** path, loop orders are rewritten so inner loops run down columns. Each matrix entry's accumulation order (over `l = 1..n`) is unchanged, so the results are expected to be **bit-identical** to the reference loops; a test enforces this. Legacy pivot mode always uses the reference loops.
- **Option:** `kernel = KERNEL_FAST (default) | KERNEL_REFERENCE`, applied as follows:
  - Fortran `band_solve(..., kernel=)` and `newton_options%kernel`;
  - the C ABI gains `bandsolver_f_solve_kernel(...)`, and `bandsolver_newton_options.kernel` is appended at the end of the struct. That is an ABI addition, noted in the CHANGELOG;
  - Python `solve(..., kernel="fast"|"reference")` and the same on `newton*`.

  It affects the Fortran backend only; the C++ loops are already contiguous.
- **Tests:** fast vs reference is bit-identical across the n/nj/X-Y sweep; results agree with C++ within tolerance; the legacy oracle still matches exactly.
- **Benchmark:** Layer 1 compares Fortran reference, Fortran fast and C++, both natively and from Python.

## B. Factor/solve split and Jacobian reuse (PR B)
- **Kernel:** `factor(system) -> Factorization` and `factor_solve(Factorization, G) -> dc`, in C++, Fortran, the C ABI and Python (`bs.factor(A,B,D,X,Y)` returning an object with `.solve(G)`).
  - The factor stores the per-node LU with partial pivoting of $\hat B_j$, the $E_j$ blocks, the effective $A_j$ (including the Y correction), $X'$ and Y.
  - One RHS then costs a forward sweep plus back substitution, O(nj·n²) instead of O(nj·n³).
  - The existing one-shot `solve` is unchanged (bit-for-bit).
- **Newton option:** `jacobian_reuse = False (default) | True`.
  - When True, the driver factors at the first iteration and reuses the factorization while updates contract (‖Δ_k‖ ≤ 0.5‖Δ_{k−1}‖). Otherwise, or after `reuse_max_iter` (default 5) iterations on one factorization, it refreshes.
  - An optional **residual-only callback** lets reuse iterations skip building the blocks: Python `newton(fill, c0, residual=...)`, C++ `NewtonOptions` + `ResidualFunction`, and the Fortran `band_problem` gets an optional `residual` type-bound procedure.
  - Result stats gain `jacobian_evaluations` and `factorizations`.
- **Tests:**
  - factor+solve matches `solve` to ≤ 1e-13 relative on the sweep, including X/Y and nj=3;
  - multiple RHS reuse works; singular blocks are reported at factor time;
  - reuse Newton converges to the same root, with fewer fills when the residual callback is given.

## C. Adaptive DAE integrator, order 1–2 (PR C; C++ and Python first, Fortran port follow-up)
- **Problem form:** $\mathbf F(t,\mathbf c,\dot{\mathbf c}) = 0$ on the BAND block structure, given as `residual(t, c, cdot) -> F`. The Jacobian is either:
  - `jacobian(t, c, cdot, alpha) -> (A,B,D,X,Y)` of $\partial F/\partial c + \alpha\,\partial F/\partial\dot c$, or
  - automatic, by finite-difference colouring of $c \mapsto F(t, c, \alpha c + \beta)$.
- **Method:** variable-step BDF of order 1–2 in fixed-leading-coefficient form, with predictor extrapolation and a local error estimate from the predictor–corrector difference.
  - The weighted RMS error norm uses rtol/atol and excludes variables flagged in an **algebraic mask**.
  - Step control: safety 0.9, growth clamped to [0.2, 5], and rejection with retry.
  - Order selection is 1 → 2 once enough history exists (max_order option).
  - Newton on each step uses the Jacobian reuse from PR B *across steps*; it refreshes on convergence failure or when α changes by more than 30%.
- **Options:**
  - `adaptive = True | False` (False = fixed Δt with the same method, which is the current standard);
  - `jacobian_reuse = True | False`;
  - `max_order = 1 | 2`;
  - `rtol`, `atol`, `dt0`, `dt_min`, `dt_max`, `newton_rtol` (default 0.1 × rtol).
- **Output:** the solution at requested times, plus stats (steps, rejected steps, residual evaluations, Jacobian evaluations, factorizations, Newton iterations).
- **Tests:**
  - the scalar linear ODE $\dot y = -\lambda y$ shows orders 1 and 2 with a fixed step;
  - the adaptive run meets its tolerance on $\dot y=-\lambda y$ and on a stiff 2-component problem;
  - the binary-electrolyte DAE matches the IDA reference;
  - with reuse, fewer Jacobian evaluations give the same answer within tolerance;
  - the algebraic mask is respected.
- **Benchmark:** Layer 2 work-precision with every combination of adaptive {on, off} × reuse {on, off} × order {1, 2}, against the current standard (fixed-step BDF2 with full Newton, and its linearized variant), with IDA and SciPy BDF as external references. Results go into `docs/benchmarks.md` (before/after) and notebook 5.

## Out of scope
- Orders 3–5; they are revisited after the benchmarks.
- BAND as an IDA linear solver.
- Event detection.
