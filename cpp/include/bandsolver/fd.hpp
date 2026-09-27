// Finite-difference Jacobians for the Appendix C block system (the "AUTOBAND" idea).
//
// The user supplies only the residual F(c) (c and F are [nj][n], row-major). Because F_j
// depends only on c_{j-1}, c_j, c_{j+1} (plus c_2 for F_0 and c_{nj-3} for F_{nj-1}),
// perturbing unknown k at every third node at once separates all derivatives, so one
// Jacobian costs 3n + 1 residual evaluations regardless of nj. Forward differences with
// step rel_step * max(|c|, typical). The residual MUST respect this stencil; if F_j
// depends on c_{j±2}, derivatives are silently mixed (check_jacobian against a trusted
// Jacobian will reveal it).
#ifndef BANDSOLVER_FD_HPP
#define BANDSOLVER_FD_HPP

#include <functional>

#include "bandsolver/band.hpp"
#include "bandsolver/newton.hpp"

namespace bandsolver {

struct FdOptions {
    double rel_step = 1.4901161193847656e-08;  // sqrt(machine epsilon)
    double typical = 1.0;                      // step floor scale for unknowns near zero
};

// Evaluate F(c) into F (n*nj values). May throw to signal failure.
using ResidualFunction = std::function<void(const double* c, double* F)>;

// Fill sys with the finite-difference blocks and G = -F(c). Returns the number of residual
// evaluations (3n + 1). Exceptions from the residual propagate; invalid sizes/options throw
// bandsolver::Error(invalid_argument).
long fd_jacobian(int n, int nj, const ResidualFunction& residual, const double* c, BlockSystem& sys,
                 const FdOptions& opts = {});

// A FillFunction that builds its blocks by finite differences of `residual`. If `evaluations`
// is non-null, it is incremented by the number of residual evaluations performed.
FillFunction fd_fill(int n, int nj, ResidualFunction residual, FdOptions opts = {}, long* evaluations = nullptr);

// Newton iteration using finite-difference Jacobians; sets NewtonResult::residual_evaluations.
NewtonResult newton_fd(int n, int nj, const ResidualFunction& residual, double* c, const NewtonOptions& opts = {},
                       const FdOptions& fd = {});

// Largest discrepancy found in one block type:
//   error = |user - fd| / max(|user|, |fd|, 1e-3 * rowscale),
// where rowscale is the largest entry (any block) of the same equation row. Entries far
// smaller than their row's dominant entry are thus judged at 1e-3 of the row scale, where
// forward-difference noise (~sqrt(eps) relative to the row) cannot masquerade as an error.
// node/row/col are 0-based (-1 if the block has no entries); X is reported at node 0 and
// Y at node nj-1.
struct JacobianMismatch {
    double error = 0.0;
    int node = -1;
    int row = -1;
    int col = -1;
    double user = 0.0;
    double fd = 0.0;
};

struct JacobianCheck {
    JacobianMismatch A, B, D, X, Y;
    double max_error() const noexcept;
};

// Compare a hand-written fill's blocks against finite differences of its own G (F = -G) at c.
// A correct Jacobian typically scores 1e-8 to 1e-5; scores above ~1e-3 indicate a bug.
JacobianCheck check_jacobian(int n, int nj, const FillFunction& fill, const double* c, const FdOptions& opts = {});

}  // namespace bandsolver

#endif
