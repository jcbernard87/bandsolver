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

**Fill contract (Newton):** at state `c`, set `G = −F(c)` and the Jacobian blocks. Output arrays arrive zeroed, so only nonzero entries need to be written.

## Fortran (`fortran/src`)

```fortran
use bandsolver_kernel   ! band_solve, block_solve, BAND_* status and PIVOT_* constants
call band_solve(n, nj, A, B, D, G, dc, status [, X=X] [, Y=Y] [, pivot=PIVOT_PARTIAL] &
                [, fail_node=node] [, min_rel_pivot=r])

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

## Python (`bandsolver`)

```python
bandsolver.solve(A, B, D, G, X=None, Y=None, *, pivot="partial", backend="cpp") -> ndarray (nj, n)

bandsolver.newton(fill, c0, *, rtol=1e-10, atol=1e-12, damping=1.0, max_iter=50,
                  pivot="partial", require_convergence=True, backend="cpp") -> NewtonResult
#   fill(c) -> (A, B, D, G) or (A, B, D, G, X, Y); c is a copy with shape (nj, n)
#   NewtonResult: c, converged, iterations, status, update_norm, step_norm, residual_norm

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
