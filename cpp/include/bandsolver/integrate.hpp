// Adaptive BDF (order 1-2) integrator for DAEs F(t, c, cdot) = 0 on the BAND block structure.
//
// Each step solves the implicit BDF equation with Newton's method using the BAND
// factor/solve split. Options (every one can be switched):
//   adaptive        true: local error control with step rejection/growth; false: fixed dt
//   max_order       1 (backward Euler) or 2 (variable-step BDF2)
//   jacobian_reuse  keep one factorization across Newton iterations and steps while it
//                   still converges and alpha = d(cdot)/dc changes by < reuse_alpha_change
// The Jacobian dF/dc + alpha dF/dcdot comes from a user callback or, if none is given, from
// finite differences (fd_jacobian colouring, 3n+1 residual evaluations).
//
// Error control uses a weighted RMS norm, weights 1/(atol + rtol|c|), over the differential
// variables (entries flagged in `algebraic` are excluded; they are still solved for). The
// local error is estimated from the corrector-predictor difference (Milne), and steps are
// shortened to land exactly on each requested output time.
#ifndef BANDSOLVER_INTEGRATE_HPP
#define BANDSOLVER_INTEGRATE_HPP

#include <exception>
#include <functional>
#include <string>
#include <vector>

#include "bandsolver/band.hpp"
#include "bandsolver/fd.hpp"

namespace bandsolver {

// F(t, c, cdot) -> F, all n*nj values, row-major [nj][n].
using DaeResidual = std::function<void(double t, const double* c, const double* cdot, double* F)>;
// Fill sys.A/B/D/X/Y with dF/dc + alpha * dF/dcdot at (t, c, cdot). sys arrives zeroed; G is ignored.
using DaeJacobian = std::function<void(double t, const double* c, const double* cdot, double alpha, BlockSystem& sys)>;

struct IntegratorOptions {
    bool adaptive = true;
    int max_order = 2;              // 1 or 2
    double rtol = 1e-6;
    double atol = 1e-8;
    double dt = 0.0;                // fixed step (adaptive = false); required > 0 there
    double dt0 = 0.0;               // initial step (adaptive); 0 = automatic
    double dt_min = 0.0;
    double dt_max = 0.0;            // 0 = unbounded
    bool jacobian_reuse = true;
    double reuse_alpha_change = 0.3;
    int max_newton_iter = 4;
    double newton_tol = 0.33;       // on the weighted norm, with a convergence-rate estimate
    int max_steps = 1000000;
    FdOptions fd;                   // used when no Jacobian callback is given
};

struct IntegratorStats {
    int steps = 0;                  // accepted steps
    int rejected_error = 0;
    int rejected_newton = 0;
    int newton_iterations = 0;
    int jacobian_evaluations = 0;   // user-Jacobian calls or FD Jacobians
    int factorizations = 0;
    long residual_evaluations = 0;  // all residual calls, including those inside FD Jacobians
};

struct IntegrationResult {
    Status status = Status::ok;
    std::vector<double> t;                   // the output times reached
    std::vector<std::vector<double>> y;      // solution at each output time (n*nj values)
    IntegratorStats stats;
    double t_reached = 0.0;
    std::string message;
    std::exception_ptr callback_exception;
};

// Integrate from t0 (state c0; cdot0 optional, used only to predict the first step) and
// return the solution at each time in t_out (strictly increasing, > t0). `algebraic` is empty
// or has n*nj flags (nonzero = algebraic, excluded from error control).
IntegrationResult integrate(int n, int nj, const DaeResidual& residual, const DaeJacobian& jacobian, double t0,
                            const double* c0, const double* cdot0, const std::vector<double>& t_out,
                            const std::vector<char>& algebraic, const IntegratorOptions& opts = {});

}  // namespace bandsolver

#endif
