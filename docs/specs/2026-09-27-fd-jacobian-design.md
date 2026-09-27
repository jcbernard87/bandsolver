# Finite-difference Jacobians — design

**Date:** 2026-09-27 · **Status:** approved by the author.

## Goal
Users should only have to supply the residual `F(c)` for the whole grid (shape `(nj, n)`). The library then builds the Appendix C blocks `A, B, D, X, Y` and sets `G = −F`. A companion `check_jacobian` verifies a hand-written `fill` against finite differences. Appendix C describes both ideas: the AUTOBAND program, and checking the programmed equations.

## Method
- **Stencil assumption.** `F_j` depends on `c_{j−1}, c_j, c_{j+1}` only. The end rows also reach one node further: `F_0` on `c_2` (the X block) and `F_{nj−1}` on `c_{nj−3}` (the Y block). All of these lie within three consecutive nodes.
- **Colouring.** Perturb unknown `k` at every node with `j ≡ r (mod 3)` simultaneously, for r = 0, 1, 2 and k = 1..n. Each residual row then sees at most one perturbed node per colour, so its derivative is unambiguous.
- **Cost.** 3n evaluations plus the base evaluation, i.e. 3n + 1 per Jacobian, independent of nj.
- **Step.** Forward differences with `h = rel_step · max(|c_jk|, typical) · sign(c_jk)`, recomputed as `(c+h) − c` so the division uses the exact step taken. Defaults: `rel_step = √ε ≈ 1.49e-8` and `typical = 1`. Both are options.
- **Extraction.** From the perturbation at node m, colour r and unknown k:
  - `B_m[:,k]` comes from row m;
  - `D_{m−1}[:,k]` from row m−1, when m−1 ≥ 0;
  - `A_{m+1}[:,k]` from row m+1, when m+1 ≤ nj−1;
  - `X[:,k]` from row 0 when m = 2;
  - `Y[:,k]` from row nj−1 when m = nj−3.
  - With nj = 3, node 2's perturbation feeds both `D_1` and `X` from rows 1 and 0, and node 0's feeds both `A_1` and `Y`, from rows 1 and 2.
- **Implementation.** The finite-difference code is an adapter that turns a residual into a `fill` function, so `newton_fd` reuses the existing Newton driver unchanged. It is implemented natively in both cores.
- **Jacobian check.** `check_jacobian(fill, c)` takes `F = −G` from the user's own fill, differentiates it by finite differences, and reports the largest scaled mismatch `|J_user − J_fd| / max(|J_user|, |J_fd|, 1e-3·rowscale)` per block type, with its node, row and column. `rowscale` is the largest entry of that equation row. This was adjusted during F1: a pure per-entry relative error flagged forward-difference noise on tiny entries.

## API
- **C++:**
  - `ResidualFunction = std::function<void(const double* c, double* F)>`
  - `fd_fill(residual, FdOptions) -> FillFunction`
  - `fd_jacobian(n, nj, residual, c, BlockSystem&, FdOptions)`
  - `newton_fd(n, nj, residual, c, NewtonOptions, FdOptions)`
  - `check_jacobian(n, nj, fill, c, FdOptions) -> JacobianCheck`
  - `FdOptions`: `rel_step` and `typical`. The result exposes `residual_evaluations`.
- **Fortran:**
  - an abstract type `band_residual_problem` with a deferred `residual(self, n, nj, c, F, ierr)`;
  - `band_fd_jacobian`, `band_newton_fd` and `band_check_jacobian`.
- **C ABI:** `bandsolver_f_fd_jacobian` and `bandsolver_f_newton_fd`.
- **Python:**
  - `fd_jacobian(residual, c)`
  - `newton_fd(residual, c0, ..., rel_step, typical, backend)`
  - `check_jacobian(fill, c, backend) -> JacobianCheck`

## Tests
- The finite-difference blocks match the analytic Jacobians of the nonlinear BVP, the n = 3 coupled DAE and the X/Y problem to ≲ 1e-6 relative, including nj = 3 with X and Y.
- `newton_fd` converges to the analytic-Jacobian solutions within 1e-10.
- Residual evaluations number exactly (3n+1) per iteration.
- `check_jacobian` passes for correct fills and locates a planted error by block, node, row and column.
- The C++ and Fortran backends agree.
- Callback errors propagate.

## Out of scope
Central differences, automatic differentiation, sparsity inside blocks, and verification of the stencil assumption (documented as a user responsibility; `check_jacobian` against a trusted Jacobian catches violations).
