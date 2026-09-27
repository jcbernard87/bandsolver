// F1: finite-difference Jacobians (C++ core). A nonlinear n=3 test problem with full
// neighbour coupling and nonlinear X/Y endpoint terms has an exact analytic Jacobian.
#include <bandsolver/fd.hpp>

#include "fd_problem.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <stdexcept>
#include <vector>

using namespace bandsolver;

static int failures = 0;
static void check(bool cond, const char* msg) {
    std::printf("%s %s\n", cond ? "PASS" : "FAIL", msg);
    if (!cond) ++failures;
}

// Max |fd - exact| / max|exact| over all blocks (relative to the Jacobian's scale).
static double jac_diff(BlockSystem& a, BlockSystem& e) {
    double d = 0, s = 0;
    const std::vector<double>* ea[] = {&e.A(), &e.B(), &e.D(), &e.X(), &e.Y()};
    const std::vector<double>* aa[] = {&a.A(), &a.B(), &a.D(), &a.X(), &a.Y()};
    for (int b = 0; b < 5; ++b)
        for (std::size_t i = 0; i < ea[b]->size(); ++i) {
            // A[0] and D[nj-1] are unused and not produced by either side.
            d = std::max(d, std::abs((*aa[b])[i] - (*ea[b])[i]));
            s = std::max(s, std::abs((*ea[b])[i]));
        }
    for (std::size_t i = 0; i < e.G().size(); ++i) d = std::max(d, std::abs(a.G()[i] - e.G()[i]) / s);
    return d / s;
}

int main() {
    for (int n : {1, 3})
        for (int nj : {3, 4, 5, 10}) {
            Problem p(n, nj, 17u + n * 100 + nj);
            auto c = start(n, nj);
            BlockSystem exact(n, nj), fd(n, nj);
            p.fill(c.data(), exact);
            long evals = fd_jacobian(n, nj, [&](const double* x, double* F) { p.residual(x, F); }, c.data(), fd);
            double d = jac_diff(fd, exact);
            char msg[128];
            std::snprintf(msg, sizeof msg, "fd_jacobian n=%d nj=%d matches analytic (rel %.1e), %ld evals", n, nj, d, evals);
            check(d < 1e-6 && evals == 3 * n + 1, msg);
        }

    {  // newton_fd reaches the analytic-Jacobian solution; evaluation count is exact.
        const int n = 3, nj = 40;
        Problem p(n, nj, 5);
        auto ca = start(n, nj), cf = ca;
        NewtonResult ra = newton(n, nj, [&](const double* c, BlockSystem& s) { p.fill(c, s); }, ca.data());
        NewtonResult rf = newton_fd(n, nj, [&](const double* c, double* F) { p.residual(c, F); }, cf.data());
        double d = 0;
        for (int i = 0; i < n * nj; ++i) d = std::max(d, std::abs(ca[i] - cf[i]));
        std::printf("analytic: %d iterations; fd: %d iterations, %ld residual evaluations; max diff %.1e\n",
                    ra.iterations, rf.iterations, rf.residual_evaluations, d);
        check(ra.converged && rf.converged && d < 1e-10, "newton_fd converges to the analytic-Jacobian solution");
        check(rf.residual_evaluations == static_cast<long>(3 * n + 1) * rf.iterations, "residual evaluations = (3n+1) per iteration");
        check(rf.iterations <= ra.iterations + 2, "newton_fd needs at most 2 more iterations than analytic");
    }

    {  // check_jacobian: clean on a correct fill, pinpoints a planted error.
        const int n = 3, nj = 10;
        Problem p(n, nj, 9);
        auto c = start(n, nj);
        JacobianCheck ok = check_jacobian(n, nj, [&](const double* x, BlockSystem& s) { p.fill(x, s); }, c.data());
        std::printf("correct fill: max error %.1e\n", ok.max_error());
        check(ok.max_error() < 1e-4, "check_jacobian passes a correct fill");
        JacobianCheck bad = check_jacobian(
            n, nj,
            [&](const double* x, BlockSystem& s) {
                p.fill(x, s);
                s.D(4, 1, 2) += 0.5;
            },
            c.data());
        std::printf("planted D error: block D node %d row %d col %d error %.2e (A %.1e, B %.1e, X %.1e, Y %.1e)\n",
                    bad.D.node, bad.D.row, bad.D.col, bad.D.error, bad.A.error, bad.B.error, bad.X.error, bad.Y.error);
        check(bad.D.node == 4 && bad.D.row == 1 && bad.D.col == 2 && bad.D.error > 1e-2, "planted error located");
        check(std::max({bad.A.error, bad.B.error, bad.X.error, bad.Y.error}) < 1e-4, "other blocks stay clean");
        JacobianCheck badx = check_jacobian(
            n, nj,
            [&](const double* x, BlockSystem& s) {
                p.fill(x, s);
                s.X(2, 0) = 0;  // forgot the X term
            },
            c.data());
        check(badx.X.row == 2 && badx.X.col == 0 && badx.X.error > 1e-2, "missing X entry detected");
    }

    {  // Errors propagate: residual exception inside newton_fd, and directly from fd_jacobian.
        std::vector<double> c(9, 0.1);
        NewtonResult r = newton_fd(3, 3, [](const double*, double*) { throw std::runtime_error("bad residual"); }, c.data());
        check(r.status == Status::callback_error && r.callback_exception, "newton_fd captures residual exceptions");
        BlockSystem s(3, 3);
        bool threw = false;
        try {
            fd_jacobian(3, 3, [](const double*, double*) { throw std::runtime_error("x"); }, c.data(), s);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw, "fd_jacobian propagates residual exceptions");
        FdOptions bad_opts;
        bad_opts.rel_step = 0;
        threw = false;
        try {
            fd_jacobian(3, 3, [](const double*, double* F) { std::fill(F, F + 9, 0.0); }, c.data(), s, bad_opts);
        } catch (const Error& e) {
            threw = e.status() == Status::invalid_argument;
        }
        check(threw, "rel_step <= 0 rejected");
    }
    return failures ? 1 : 0;
}
