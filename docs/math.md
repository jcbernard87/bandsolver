# Mathematical formulation

Reference: J. Newman and K. E. Thomas-Alyea, *Electrochemical Systems*, Appendix C ("Numerical Solution of Coupled, Ordinary Differential Equations"), §§C.1–C.4. The correspondence between that listing and the author's archival battery-model code was checked in the author's private source audit, which is not included here.

## Problem class

The library solves systems where each node carries `n` unknowns `c_j ∈ ℝⁿ`, `j = 0 … nj−1`. The residual at node `j` depends only on the neighbouring nodes:

```
F_j(c_{j-1}, c_j, c_{j+1}) = 0,
```

The first and last nodes may also reach one node further, to `c₂` and `c_{nj−3}` respectively. This covers one-dimensional finite-difference and finite-volume discretizations of coupled ODE/PDE systems, including:
- mixed differential and algebraic rows, such as electroneutrality or potential equations;
- interface conditions;
- the equations of each implicit time step.

## Newton linearization

Let `c` be the current iterate and `Δc` the correction. Linearizing `F(c + Δc) = 0` gives the block system

```
A_j Δc_{j-1} + B_j Δc_j + D_j Δc_{j+1} = G_j,
A_j = ∂F_j/∂c_{j-1},   B_j = ∂F_j/∂c_j,   D_j = ∂F_j/∂c_{j+1},   G_j = −F_j(c),
```

Two endpoint blocks extend this: `X = ∂F_0/∂c_2` and `Y = ∂F_{nj−1}/∂c_{nj−3}`. The update is `c ← c + λΔc`, where λ ∈ (0, 1] is the damping (default 1). Because `G` is the negative residual, it vanishes at the solution. Appendix C recommends this "change variable" form because the right-hand side then directly measures how far the equations are from being satisfied.

**Convergence test** (after each update): `max_ij |Δc_ij| / (atol + rtol·|c_ij|) ≤ 1`, with defaults `rtol = 1e-10` and `atol = 1e-12`. The library also records `max|G|` at each linearization. Near a solution with a correct Jacobian and λ = 1, convergence is quadratic (see `validation.md`). Convergence from a poor initial guess is not guaranteed, as Appendix C also cautions. Damping, a better initial guess, or continuation are the user's tools there.

Setting `max_iter = 1` and `require_convergence = false` performs one linearized correction per call. The archival battery models use it this way, as a semi-implicit time step (audit §8.10). The result is not a converged Newton solve.

## Elimination (the BAND recurrence)

Forward sweep: seek `Δc_j = E_j Δc_{j+1} + e_j`.

- **Node 0.** Solve `B₀ [S_D | S_X | s_G] = [D₀ | X | G₀]`. Then `E₀ = −S_D`, `e₀ = s_G`, and `X′ = −S_X`. This gives `Δc₀ = E₀Δc₁ + e₀ + X′Δc₂`.
- **Node 1.** Substituting `Δc₀` adds `A₁X′` to the upper block: `D₁ ← D₁ + A₁X′`.
- **Last node (Y).** Substitute `Δc_{nj−3} = E_{nj−3}Δc_{nj−2} + e_{nj−3}`. This gives `A ← A + Y E_{nj−3}` and `G ← G − Y e_{nj−3}`. If `nj = 3`, node 0's `X′Δc₂` term falls on the last node itself, so also `B ← B + Y X′`.
- **Interior and last nodes.** Form `B̂_j = B_j + A_j E_{j−1}`, then solve `B̂_j [S | s] = [D_j | G_j − A_j e_{j−1}]`, giving `E_j = −S` and `e_j = s`.

Back substitution: `Δc_{nj−1} = e_{nj−1}`, then `Δc_j = E_jΔc_{j+1} + e_j` going down to `j = 0`, and finally `Δc₀ += X′Δc₂`.

Each node's block system is a dense `n×n` solve with `n+1` right-hand sides; node 0 has `2n+1`. The total cost is O(nj·n³) flops and O(nj·n²) storage for the `E_j`. This is block LU without pivoting across nodes. It is stable when the `B̂_j` are well conditioned, which is typical for diffusion-dominated discretizations. If a `B̂_j` is singular, the method stops at that node rather than pivoting across nodes. `min_rel_pivot` is reported so that ill conditioning can be seen.

## Finite-difference Jacobians

Appendix C suggests an "AUTOBAND" program that computes the coefficients `A, B, D` by numerical derivatives of the governing equations. The library does this in `fd_jacobian` and `newton_fd`.

**Colouring.** `F_j` depends only on nodes `j−1, j, j+1`, plus the endpoint reach to `c₂` and `c_{nj−3}`. So every equation row involves at most three consecutive nodes. Perturbing unknown `k` at all nodes with `j ≡ r (mod 3)` at once therefore changes each row through at most one perturbed node. For one colour `r` and unknown `k`, the change in rows `m−1, m, m+1` gives column `k` of `D_{m−1}`, `B_m` and `A_{m+1}`. The endpoint rows give `X` (row 0, when m = 2) and `Y` (row nj−1, when m = nj−3). One full Jacobian costs **3n + 1** residual evaluations, independent of `nj`.

**Step size.** Forward differences use `h = rel_step · max(|c|, typical) · sign(c)`, with `rel_step = √ε` by default. The divisor is the step actually representable, `(c + h) − c`. Appendix C notes the trade-off: too small a step amplifies round-off, and too large a step gives a poor derivative. The expected error is about √ε relative to the magnitude of each equation row. For badly scaled unknowns, set `typical` to their characteristic size.

**Requirement.** The residual must respect the stencil. If `F_j` depends on `c_{j±2}`, the colouring silently mixes derivatives.

**Checking a hand-written Jacobian.** `check_jacobian(fill, c)` differentiates the fill's own `G` (`F = −G`) and compares each entry with the user's:

```
score = |J_user − J_fd| / max(|J_user|, |J_fd|, 10⁻³ · rowscale)
```

Here `rowscale` is the largest Jacobian entry of that equation row. Forward-difference noise is about √ε·rowscale, so correct Jacobians score about 1e-9 to 1e-5. Mistakes in entries that matter at the 10⁻³ row level score well above 10⁻³.

## Block solve and pivoting

- **`partial` (default).** Gauss–Jordan elimination with row partial pivoting on each block.
- **`legacy`.** The Appendix C `MATINV` pivot rule. For each unused row, it finds the largest and second-largest magnitudes among the unused columns. It pivots on the row whose second/first ratio is smallest, and moves that row into the pivot column's position. Operation order follows the archival source, so results agree bit for bit with it on the test platform. This mode exists for historical comparison, not because it is more accurate.

In both modes, a block is **singular** if a chosen pivot satisfies `|p| ≤ n·ε·max|block|`.

## Differences from the archival kernel

| Behaviour | Archival `BAND`/`MATINV` | bandsolver |
|---|---|---|
| Singular block | Prints `DETERM=0 AT J=` and continues with an undefined result | Returns status `SINGULAR` with the node index |
| Rank-deficient block after rounding | Not detected (exact-zero test only) | Detected by the relative pivot threshold |
| `nj = 3` with nonzero X and Y | Drops node 0's `X` term when eliminating `Y`, giving a wrong answer (backward error ≈ 3e-3) | Correct |
| Inputs | Overwrites module-global A, B, D, G, X | Inputs unchanged; no global state |
| NaN/Inf | Propagates silently | Status `NON_FINITE` |
| Nonlinear iteration | The caller loops; the archival models apply one correction per time step | Optional Newton driver with convergence test, damping, and histories |
