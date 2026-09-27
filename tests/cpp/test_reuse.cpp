// P4: Newton with Jacobian reuse (modified Newton) reaches the same root as full Newton with
// fewer factorizations; with a residual-only callback it also makes fewer fill calls.
#include <bandsolver/newton.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "fd_problem.hpp"

using namespace bandsolver;

static int failures = 0;
static void check(bool cond, const char* msg) {
    std::printf("%s %s\n", cond ? "PASS" : "FAIL", msg);
    if (!cond) ++failures;
}

int main() {
    const int n = 3, nj = 40;
    Problem p(n, nj, 5);
    auto start = [&] {
        std::vector<double> c(n * nj);
        for (int i = 0; i < n * nj; ++i) c[i] = 0.3 + 0.01 * i;
        return c;
    };
    long fills = 0, residuals = 0;
    FillFunction fill = [&](const double* x, BlockSystem& s) { ++fills; p.fill(x, s); };
    ResidualFunction residual = [&](const double* x, double* F) { ++residuals; p.residual(x, F); };

    auto cf = start();
    NewtonResult full = newton(n, nj, fill, cf.data());
    std::printf("full Newton: %d iterations, %d factorizations, %d fill calls\n", full.iterations,
                full.factorizations, full.jacobian_evaluations);
    check(full.converged && full.factorizations == full.iterations, "full Newton factors every iteration");

    NewtonOptions o;
    o.jacobian_reuse = true;
    auto cr = start();
    fills = residuals = 0;
    NewtonResult reuse = newton(n, nj, fill, cr.data(), o);
    double d = 0;
    for (int i = 0; i < n * nj; ++i) d = std::max(d, std::abs(cr[i] - cf[i]));
    std::printf("reuse (fill only): %d iterations, %d factorizations, %d fill calls, max diff %.1e\n",
                reuse.iterations, reuse.factorizations, reuse.jacobian_evaluations, d);
    check(reuse.converged && d < 1e-9, "reuse converges to the same root");
    check(reuse.factorizations < reuse.iterations, "reuse needs fewer factorizations than iterations");

    auto cr2 = start();
    fills = residuals = 0;
    NewtonResult reuse2 = newton(n, nj, fill, cr2.data(), o, residual);
    d = 0;
    for (int i = 0; i < n * nj; ++i) d = std::max(d, std::abs(cr2[i] - cf[i]));
    std::printf("reuse (+ residual callback): %d iterations, %d factorizations, %ld fills, %ld residual calls, diff %.1e\n",
                reuse2.iterations, reuse2.factorizations, fills, residuals, d);
    check(reuse2.converged && d < 1e-9, "reuse with residual callback converges to the same root");
    check(fills == reuse2.factorizations && residuals == reuse2.iterations - reuse2.factorizations &&
              reuse2.residual_evaluations == residuals,
          "fill is called only to refactor; other iterations use the residual callback");

    NewtonOptions o1 = o;
    o1.reuse_max_iter = 1;
    auto c1 = start();
    NewtonResult r1 = newton(n, nj, fill, c1.data(), o1);
    check(r1.factorizations == r1.iterations, "reuse_max_iter = 1 refactors every iteration (full Newton)");

    NewtonOptions ostrict = o;
    ostrict.reuse_contraction = 1e-12;          // any step that doesn't collapse forces a refresh
    auto c2 = start();
    NewtonResult r2 = newton(n, nj, fill, c2.data(), ostrict);
    check(r2.converged && r2.factorizations >= reuse.factorizations, "a strict contraction threshold refreshes more often");

    NewtonOptions bad = o;
    bad.reuse_max_iter = 0;
    auto c3 = start();
    check(newton(n, nj, fill, c3.data(), bad).status == Status::invalid_argument, "reuse_max_iter < 1 rejected");
    return failures ? 1 : 0;
}
