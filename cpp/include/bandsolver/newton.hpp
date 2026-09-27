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
    long residual_evaluations = 0;       // set by newton_fd only
};

using FillFunction = std::function<void(const double* c, BlockSystem& sys)>;

// Iterate in place on c (n*nj values). Errors are reported through the result.
NewtonResult newton(int n, int nj, const FillFunction& fill, double* c, const NewtonOptions& opts = {});

}  // namespace bandsolver

#endif
