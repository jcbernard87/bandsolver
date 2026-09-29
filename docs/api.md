# API reference

## Conventions shared by all interfaces

| | Fortran | C / C++ / Python |
|---|---|---|
| Blocks `A, B, D` | `A(n,n,nj)`; `A(i,k,j)` couples equation `i` of node `j` to unknown `k` of node `j−1` | row-major `[nj][n][n]`; numpy shape `(nj, n, n)` |
| `G`, `c`, `dc` | `G(n,nj)` | `[nj][n]`; numpy `(nj, n)` (same memory layout as Fortran) |
| `X`, `Y` | `X(n,n)`, optional | `[n][n]`, nullable / `None` |
| Node index in errors | 1-based | C ABI 1-based; C++ and Python 0-based |
| Unused | `A(:,:,1)`, `D(:,:,nj)` | `A[0]`, `D[nj−1]` |
| Requirement | `n ≥ 1`, `nj ≥ 3` | same |

**Status codes**, identical everywhere: `0 OK`, `1 SINGULAR`, `2 INVALID_ARGUMENT`, `3 NOT_CONVERGED`, `4 NON_FINITE`, `5 CALLBACK_ERROR`.

**Pivot modes:** `0 partial` (default), `1 legacy`. See [math.md](math.md).

**Singular-block rule** (`solve` only, added in 0.1.2): `0 relative` (default; a pivot ≤ n·ε·max|block| is singular) or `1 exact` (only an exactly zero pivot, as in the archival `MATINV`). It is available in Fortran as `band_solve(..., singular=SINGULAR_EXACT)`, in C as `bandsolver_f_solve_ex`, in C++ as `solve(view, dc, pivot, Singular::exact)`, and in Python as `singular="exact"`. `factor`, `newton` and `integrate` always use the relative rule.

**Fortran kernel** (Fortran backend only): `0 fast` (default; column-major loops) or `1 reference` (archival loop order). With gfortran the two give bit-identical results; other compilers agree to rounding. `reference` exists for comparison and benchmarking. The option is available in Fortran as `band_solve(..., kernel=)` and `newton_options%kernel`, in C as `bandsolver_f_solve_kernel` and `bandsolver_newton_options.kernel`, and in Python as `kernel="fast"|"reference"`.

**Fill contract (Newton):** at state `c`, set `G = −F(c)` and the Jacobian blocks. Output arrays arrive zeroed, so only nonzero entries need to be written.

## Fortran (`fortran/src`)

```fortran
use bandsolver_kernel   ! band_solve, block_solve, BAND_* status and PIVOT_* constants
call band_solve(n, nj, A, B, D, G, dc, status [, X=X] [, Y=Y] [, pivot=PIVOT_PARTIAL] &
                [, fail_node=node] [, min_rel_pivot=r] [, kernel=KERNEL_FAST] [, singular=SINGULAR_RELATIVE])

use bandsolver_newton   ! band_problem, newton_options, newton_result, band_newton
type, extends(band_problem) :: my_problem
    ! your data
contains
    procedure :: fill      ! subroutine fill(self, n, nj, c, A, B, D, G, X, Y, ierr)
end type
type(newton_options) :: opts       ! rtol=1e-10, atol=1e-12, damping=1, max_iter=50,
                                   ! pivot=PIVOT_PARTIAL, require_convergence=.true.
type(newton_result) :: res         ! status, iterations, fail_node, converged,
                                   ! update_norm(:), step_norm(:), residual_norm(:)
call band_newton(problem, n, nj, c, opts, res)     ! c is updated in place
```

A nonzero `ierr` from `fill` stops the iteration with `BAND_CALLBACK_ERROR`. There is no module-level mutable state, so separate `band_problem` instances may be solved concurrently.

## C (`fortran/include/bandsolver_f.h`, backed by the Fortran library)

```c
int bandsolver_f_solve(int n, int nj, const double *A, const double *B, const double *D,
                       const double *G, const double *X, const double *Y, int pivot,
                       double *dc, int *fail_node, double *min_rel_pivot);
/* every option: pivot, kernel (enum bandsolver_kernel), singular (enum bandsolver_singular) */
int bandsolver_f_solve_ex(int n, int nj, const double *A, const double *B, const double *D,
                          const double *G, const double *X, const double *Y, int pivot, int kernel,
                          int singular, double *dc, int *fail_node, double *min_rel_pivot);

typedef int (*bandsolver_fill_fn)(int n, int nj, const double *c, double *A, double *B,
                                  double *D, double *G, double *X, double *Y, void *ctx);
void bandsolver_f_default_options(bandsolver_newton_options *opts);
int bandsolver_f_newton(int n, int nj, bandsolver_fill_fn fill, void *ctx, double *c,
                        const bandsolver_newton_options *opts, bandsolver_newton_result *res,
                        double *update_history, double *step_history, double *residual_history);
```

History arrays (length ≥ `max_iter`) may be `NULL`. The fill callback returns nonzero to abort. Link against `bandsolver_f` and the Fortran runtime; CMake adds the runtime automatically.

## C++ (`cpp/include/bandsolver`)

```cpp
#include <bandsolver/band.hpp>
bandsolver::SystemView v{n, nj, A, B, D, G, X /*nullable*/, Y /*nullable*/};
bandsolver::SolveInfo info = bandsolver::solve(v, dc, bandsolver::Pivot::partial);
// or solve(v, dc, bandsolver::Pivot::legacy, bandsolver::Singular::exact)
//   info.status, info.fail_node (0-based or -1), info.min_rel_pivot

bandsolver::BlockSystem sys(n, nj);           // owning, zeroed; sys.A(j,i,k), sys.G(j,i), sys.X(i,k) ...
std::vector<double> dc = bandsolver::solve(sys);   // throws bandsolver::Error{status(), node()}

#include <bandsolver/newton.hpp>
bandsolver::NewtonOptions opts;               // rtol, atol, damping, max_iter, pivot, require_convergence
bandsolver::NewtonResult r = bandsolver::newton(n, nj,
    [&](const double* c, bandsolver::BlockSystem& s) { /* fill s, G = -F(c) */ }, c, opts);
//   r.status, r.iterations, r.converged, r.fail_node, r.update_norm, r.step_norm,
//   r.residual_norm, r.callback_exception (rethrow with std::rethrow_exception)
```

`solve(SystemView, …)` and `newton` report errors through their return values. Only `std::bad_alloc` can propagate out of them.

## Factor once, solve many

- **C++:** `bandsolver::factor(SystemView) -> Factorization` (`status()`, `fail_node()`, `solve(G, dc)`).
- **Fortran:** `band_factor(n, nj, A, B, D, f [, X, Y])` and `band_factor_solve(f, G, dc, status)`.
- **C:** `bandsolver_f_factor` / `_factor_solve` / `_factor_free`, using an opaque handle.
- **Python:** `bs.factor(A, B, D, X, Y, backend=) -> Factorization` with `.solve(G)`.

A factorization costs 1–2× a one-shot solve; each re-solve is O(nj·n²).

## Jacobian reuse (modified Newton)

The option is off by default.
- **Fortran:** `newton_options%jacobian_reuse` (with `reuse_max_iter` and `reuse_contraction`). Override `band_problem%residual` to supply a residual-only evaluation.
- **C++:** `NewtonOptions::jacobian_reuse` etc., plus an optional `ResidualFunction` argument to `newton()`.
- **C:** fields appended to `bandsolver_newton_options`, and `bandsolver_f_newton_ex(..., residual, ...)`.
- **Python:** `newton(..., jacobian_reuse=True, residual=f)` and `newton_fd(..., jacobian_reuse=True)`.

In every interface, results gain `jacobian_evaluations` and `factorizations`.

## Adaptive DAE integrator (C++ core, Python)

This solves $F(t, c, \dot c) = 0$ with variable-step BDF of order 1–2 (`<bandsolver/integrate.hpp>`, `bs.integrate`). A Fortran port is planned.

```python
r = bs.integrate(residual, c0, t_out, *, jacobian=None, algebraic=None, cdot0=None, t0=0.0,
                 adaptive=True, max_order=2, rtol=1e-6, atol=1e-8, dt=None, dt0=None, dt_min=0, dt_max=None,
                 jacobian_reuse=True, reuse_alpha_change=0.3, max_newton_iter=4, newton_tol=0.33,
                 max_steps=1_000_000)
#   residual(t, c, cdot) -> F ; jacobian(t, c, cdot, alpha) -> (A, B, D[, X, Y]) of dF/dc + alpha dF/dcdot
#   r.t, r.y (len(t), nj, n), r.stats {steps, rejected_error, rejected_newton, newton_iterations,
#                                      jacobian_evaluations, factorizations, residual_evaluations}
```

**Options:**
- `adaptive=False` gives fixed-step BDF with step `dt`.
- `jacobian_reuse=False` refactors at every Newton iteration.
- If `jacobian` is omitted, the Jacobian is built by finite-difference colouring.
- `algebraic` flags the entries (for example potentials) that are excluded from error control.

**Errors:** an early stop raises `IntegrationError` with `.result`.

In C++, the equivalent is `integrate(n, nj, residual, jacobian, t0, c0, cdot0, t_out, algebraic, IntegratorOptions)`, which returns an `IntegrationResult`.

## Finite-difference Jacobians (all interfaces)

These build `A, B, D, X, Y` and `G = −F` from a residual `F(c)` using 3n + 1 residual evaluations; see [math.md](math.md#finite-difference-jacobians). Options are `rel_step` (default √ε) and `typical` (default 1).

**C++** (`<bandsolver/fd.hpp>`):
```cpp
using ResidualFunction = std::function<void(const double* c, double* F)>;   // [nj][n]
long evals = bandsolver::fd_jacobian(n, nj, residual, c, sys /*BlockSystem&*/, fd_opts);
bandsolver::NewtonResult r = bandsolver::newton_fd(n, nj, residual, c, newton_opts, fd_opts);
//   r.residual_evaluations == (3n+1) * r.iterations
bandsolver::JacobianCheck chk = bandsolver::check_jacobian(n, nj, fill, c, fd_opts);
//   chk.A/B/D/X/Y: {error, node, row, col, user, fd} (0-based); chk.max_error()
bandsolver::FillFunction f = bandsolver::fd_fill(n, nj, residual, fd_opts, &eval_counter);
```

**Fortran** (`use bandsolver_fd`):
```fortran
type, extends(band_residual_problem) :: my_residual
contains
    procedure :: residual   ! subroutine residual(self, n, nj, c, F, ierr)
end type
call band_fd_jacobian(prob, n, nj, c, A, B, D, G, X, Y, status [, opts=fd_options(...)] [, evaluations=k])
call band_newton_fd(prob, n, nj, c, newton_opts, res [, fd_opts])   ! res%residual_evaluations
call band_check_jacobian(fill_problem, n, nj, c, check, status [, opts])  ! 1-based locations; max_error(check)
```

**C:**
```c
bandsolver_f_default_fd_options(&fd);
bandsolver_f_fd_jacobian(n, nj, residual_fn, ctx, c, &fd, A, B, D, G, X, Y, &evals);
bandsolver_f_newton_fd(n, nj, residual_fn, ctx, c, &opts, &fd, &res, NULL, NULL, NULL, &evals);
bandsolver_f_check_jacobian(n, nj, fill_fn, ctx, c, &fd, &check);           // 1-based locations
```
Here `residual_fn` has the signature `int (*)(int n, int nj, const double* c, double* F, void* ctx)`, and `fd` may be `NULL`.

## Python (`bandsolver`)

```python
bandsolver.solve(A, B, D, G, X=None, Y=None, *, pivot="partial", backend="cpp", kernel="fast",
                 singular="relative") -> ndarray (nj, n)

bandsolver.newton(fill, c0, *, rtol=1e-10, atol=1e-12, damping=1.0, max_iter=50,
                  pivot="partial", require_convergence=True, backend="cpp") -> NewtonResult
#   fill(c) -> (A, B, D, G) or (A, B, D, G, X, Y); c is a copy with shape (nj, n)
#   NewtonResult: c, converged, iterations, status, update_norm, step_norm, residual_norm

bandsolver.fd_jacobian(residual, c, *, rel_step=sqrt(eps), typical=1.0, backend="cpp")
    -> (A, B, D, G, X, Y)
bandsolver.newton_fd(residual, c0, *, <newton options>, rel_step, typical, backend) -> NewtonResult
#   residual(c) -> F with shape (nj, n); NewtonResult.residual_evaluations
bandsolver.check_jacobian(fill, c, *, rel_step, typical, backend) -> JacobianCheck
#   .A/.B/.D/.X/.Y: JacobianMismatch(error, node, row, col, user, fd), 0-based
#   .max_error, .worst() -> (block_name, JacobianMismatch)

bandsolver.BACKENDS == ("cpp", "fortran")
```

**Errors.**
- `SingularBlockError` has a 0-based `.node`.
- `NonFiniteError` signals NaN or Inf values.
- `NotConvergedError` carries the last iterate in `.result`.
- All three derive from `BandError`, which has `.status`.
- Invalid arguments raise `ValueError`.
- An exception raised inside `fill` propagates unchanged on both backends.

**Performance notes.**
- Inputs are converted to contiguous float64 arrays; that is a copy only if they are not already contiguous float64.
- `solve` releases the GIL during the native solve.
- `newton` calls `fill` once per iteration while holding the GIL, and each call copies the returned arrays into native storage. A vectorized numpy `fill` keeps that overhead small compared with the O(nj·n³) solve for moderate `n`.
- For the tightest loops, write `fill` in C++ or Fortran.
