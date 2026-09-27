// Newton iteration on the Appendix C block system — C++17 core.
//
// `fill(c, sys)` evaluates, at state c ([nj][n], row-major), the Jacobian blocks and
// G = -F(c) into `sys`, which arrives zeroed. Each iteration solves K dc = G and updates
// c <- c + damping*dc; it has converged when max |dc|/(atol + rtol|c|) <= 1 after the
// update. max_iter = 1 with require_convergence = false reproduces the archival usage of
// one linearized correction per call. An exception thrown by `fill` stops the iteration
// with Status::callback_error; it is kept in NewtonResult::callback_exception.
#ifndef BANDSOLVER_NEWTON_HPP
#define BANDSOLVER_NEWTON_HPP

#include <exception>
#include <functional>
#include <vector>

#include "bandsolver/band.hpp"

namespace bandsolver {

struct NewtonOptions {
    double rtol = 1e-10;
    double atol = 1e-12;
    double damping = 1.0;  // in (0, 1]
    int max_iter = 50;
    Pivot pivot = Pivot::partial;
    bool require_convergence = true;
    // Jacobian reuse (modified Newton): factor once, keep solving with that factorization
    // while the updates contract (step_k <= reuse_contraction * step_{k-1}) and it has been
    // used fewer than reuse_max_iter times; otherwise refresh. Off by default.
    bool jacobian_reuse = false;
    int reuse_max_iter = 5;
    double reuse_contraction = 0.5;
};

struct NewtonResult {
    Status status = Status::ok;
    int iterations = 0;
    bool converged = false;
    int fail_node = -1;                  // 0-based, when status == singular
    std::vector<double> update_norm;     // scaled update norm per iteration
    std::vector<double> step_norm;       // max |dc| per iteration
    std::vector<double> residual_norm;   // max |G| at the start of each iteration
    std::exception_ptr callback_exception;
    long residual_evaluations = 0;       // residual-only callback calls (reuse) / newton_fd count
    int jacobian_evaluations = 0;        // fill calls (blocks + G)
    int factorizations = 0;              // block factorizations (full Newton: one per iteration)
};

using FillFunction = std::function<void(const double* c, BlockSystem& sys)>;

// Evaluate F(c) into F (n*nj values). May throw to signal failure.
using ResidualFunction = std::function<void(const double* c, double* F)>;

// Iterate in place on c (n*nj values). Errors are reported through the result. With
// opts.jacobian_reuse, an optional residual-only callback lets iterations that reuse the
// factorization skip building the blocks (otherwise fill is called and its blocks ignored).
NewtonResult newton(int n, int nj, const FillFunction& fill, double* c, const NewtonOptions& opts = {},
                    const ResidualFunction& residual = {});

}  // namespace bandsolver

#endif
