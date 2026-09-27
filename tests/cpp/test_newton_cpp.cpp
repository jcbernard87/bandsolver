// T5: C++ Newton driver — mirrors the Fortran T3 tests and cross-checks the Fortran
// driver (via its C ABI) on the same C-callback problem.
#include <bandsolver/newton.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

#include "bandsolver_f.h"

using namespace bandsolver;

static int failures = 0;
static void check(bool cond, const char* msg) {
    std::printf("%s %s\n", cond ? "PASS" : "FAIL", msg);
    if (!cond) ++failures;
}

static const double pi = std::acos(-1.0);

// -(D(c) c')' = f on [0,1], c(0)=c(1)=0, manufactured c* = sin(pi x). D = 1 + nl*c^2.
struct Bvp { double nl; };
static int fill_bvp(int, int nj, const double* c, double* A, double* B, double* D, double* G,
                    double*, double*, void* ctx) {
    const double nl = static_cast<Bvp*>(ctx)->nl, h = 1.0 / (nj - 1);
    auto Dc = [nl](double v) { return 1 + nl * v * v; };
    auto Dd = [nl](double v) { return 2 * nl * v; };
    B[0] = 1; G[0] = -c[0];
    B[nj - 1] = 1; G[nj - 1] = -c[nj - 1];
    for (int j = 1; j < nj - 1; ++j) {
        const double x = j * h, s = std::sin(pi * x), sx = pi * std::cos(pi * x), sxx = -pi * pi * s;
        const double f = -(Dd(s) * sx * sx + Dc(s) * sxx);
        const double cm = (c[j] + c[j - 1]) / 2, cp = (c[j] + c[j + 1]) / 2;
        const double Dminus = Dc(cm), Dplus = Dc(cp);
        const double F = -(Dplus * (c[j + 1] - c[j]) - Dminus * (c[j] - c[j - 1])) / (h * h) - f;
        G[j] = -F;
        A[j] = (0.5 * Dd(cm) * (c[j] - c[j - 1]) - Dminus) / (h * h);
        D[j] = -(Dplus + 0.5 * Dd(cp) * (c[j + 1] - c[j])) / (h * h);
        B[j] = (Dplus + Dminus - 0.5 * Dd(cp) * (c[j + 1] - c[j]) + 0.5 * Dd(cm) * (c[j] - c[j - 1])) / (h * h);
    }
    return 0;
}

static FillFunction wrap(Bvp& p) {
    return [&p](const double* c, BlockSystem& s) {
        fill_bvp(s.n(), s.nj(), c, s.A().data(), s.B().data(), s.D().data(), s.G().data(), s.X().data(),
                 s.Y().data(), &p);
    };
}

int main() {
    const int nj = 101;
    Bvp lin{0.0}, nl{1.0};

    std::vector<double> c(nj, 0.0);
    NewtonResult r = newton(1, nj, wrap(lin), c.data());
    check(r.status == Status::ok && r.converged && r.iterations == 2, "linear converges at iteration 2");
    check(r.residual_norm[1] < 1e-10 * r.residual_norm[0], "linear residual eliminated by first solve");

    std::fill(c.begin(), c.end(), 0.0);
    r = newton(1, nj, wrap(nl), c.data());
    std::printf("nonlinear: iterations %d, step norms:", r.iterations);
    for (double s : r.step_norm) std::printf(" %.2e", s);
    std::printf("\n");
    check(r.status == Status::ok && r.converged && r.iterations <= 8, "nonlinear converges in <= 8 iterations");
    bool quad = true;
    for (std::size_t k = 0; k + 1 < r.step_norm.size(); ++k)
        if (r.step_norm[k] < 0.1 && r.step_norm[k + 1] > 1e-13) quad = quad && r.step_norm[k + 1] <= 10 * r.step_norm[k] * r.step_norm[k];
    check(quad, "quadratic convergence: s(k+1) <= 10 s(k)^2 once s(k) < 0.1");
    double err = 0;
    for (int j = 0; j < nj; ++j) err = std::max(err, std::abs(c[j] - std::sin(pi * j / (nj - 1.0))));
    std::printf("max error vs sin(pi x): %.2e\n", err);
    check(err < 1e-3, "nonlinear solution within O(h^2) of exact");

    // Cross-check against the Fortran driver on the identical callback.
    std::vector<double> cf(nj, 0.0), hu(50), hs(50), hr(50);
    bandsolver_newton_options fo; bandsolver_newton_result fr;
    bandsolver_f_default_options(&fo);
    int st = bandsolver_f_newton(1, nj, fill_bvp, &nl, cf.data(), &fo, &fr, hu.data(), hs.data(), hr.data());
    double d = 0;
    for (int j = 0; j < nj; ++j) d = std::max(d, std::abs(c[j] - cf[j]));
    bool hist_same = fr.iterations == r.iterations;
    for (int k = 0; hist_same && k < fr.iterations; ++k) hist_same = std::abs(hs[k] - r.step_norm[k]) <= 1e-12 * hs[0];
    std::printf("C++ vs Fortran Newton: max |c diff| %.2e, iterations %d vs %d\n", d, r.iterations, fr.iterations);
    check(st == 0 && d < 1e-13 && hist_same, "C++ and Fortran Newton agree (solution and step history)");

    NewtonOptions o;
    o.pivot = Pivot::legacy;
    std::vector<double> cl(nj, 0.0);
    r = newton(1, nj, wrap(nl), cl.data(), o);
    d = 0;
    for (int j = 0; j < nj; ++j) d = std::max(d, std::abs(c[j] - cl[j]));
    check(r.converged && d < 1e-12, "legacy pivot Newton agrees");

    o = NewtonOptions{};
    o.max_iter = 1; o.require_convergence = false;
    std::fill(c.begin(), c.end(), 0.0);
    r = newton(1, nj, wrap(nl), c.data(), o);
    check(r.status == Status::ok && r.iterations == 1 && !r.converged, "one-step legacy usage returns ok");
    o.require_convergence = true;
    std::fill(c.begin(), c.end(), 0.0);
    r = newton(1, nj, wrap(nl), c.data(), o);
    check(r.status == Status::not_converged, "max_iter exhausted -> not_converged");

    r = newton(1, nj, [](const double*, BlockSystem&) { throw std::runtime_error("boom"); }, c.data());
    bool rethrown = false;
    try {
        if (r.callback_exception) std::rethrow_exception(r.callback_exception);
    } catch (const std::runtime_error& e) {
        rethrown = std::string(e.what()) == "boom";
    }
    check(r.status == Status::callback_error && rethrown, "callback exception captured and rethrowable");

    o = NewtonOptions{};
    o.damping = 0;
    check(newton(1, nj, wrap(nl), c.data(), o).status == Status::invalid_argument, "damping = 0 rejected");
    check(newton(1, nj, FillFunction{}, c.data()).status == Status::invalid_argument, "empty fill rejected");

    // Damped Newton still converges (linearly) to the same solution.
    o = NewtonOptions{};
    o.damping = 0.5; o.max_iter = 200;
    std::fill(cl.begin(), cl.end(), 0.0);
    r = newton(1, nj, wrap(nl), cl.data(), o);
    std::fill(c.begin(), c.end(), 0.0);
    newton(1, nj, wrap(nl), c.data());
    d = 0;
    for (int j = 0; j < nj; ++j) d = std::max(d, std::abs(c[j] - cl[j]));
    std::printf("damped (0.5): iterations %d, diff to undamped %.2e\n", r.iterations, d);
    check(r.converged && d < 1e-9, "damped Newton converges to the same solution");

    return failures ? 1 : 0;
}
