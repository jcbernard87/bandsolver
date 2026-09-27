// P5: adaptive BDF(1-2) DAE integrator. Exact-solution tests for orders, tolerance control,
// stiffness, algebraic variables, coupled PDE blocks, Jacobian reuse and FD Jacobians.
#include <bandsolver/integrate.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace bandsolver;

static int failures = 0;
static void check(bool cond, const char* msg) {
    std::printf("%s %s\n", cond ? "PASS" : "FAIL", msg);
    if (!cond) ++failures;
}
static const double pi = std::acos(-1.0);

// Decoupled linear ODEs y_j' = -lam_j y_j on nj = 3 nodes, n = 1: F = cdot + lam y.
static const double lam[3] = {1.0, 3.0, 10.0};
static void lin_res(double, const double* y, const double* yd, double* F) {
    for (int j = 0; j < 3; ++j) F[j] = yd[j] + lam[j] * y[j];
}
static void lin_jac(double, const double*, const double*, double alpha, BlockSystem& s) {
    for (int j = 0; j < 3; ++j) s.B(j, 0, 0) = lam[j] + alpha;
}
static double lin_err(const IntegrationResult& r, double t) {
    double e = 0;
    for (int j = 0; j < 3; ++j) e = std::max(e, std::abs(r.y.back()[j] - std::exp(-lam[j] * t)));
    return e;
}

int main() {
    std::vector<double> y0 = {1.0, 1.0, 1.0};
    // 1. Fixed-step orders.
    for (int order : {1, 2}) {
        double prev = 0, ratio = 0;
        for (double dt : {0.02, 0.01, 0.005}) {
            IntegratorOptions o;
            o.adaptive = false; o.dt = dt; o.max_order = order;
            auto r = integrate(1, 3, lin_res, lin_jac, 0.0, y0.data(), nullptr, {1.0}, {}, o);
            double e = lin_err(r, 1.0);
            if (prev > 0) ratio = prev / e;
            prev = e;
            if (r.status != Status::ok) ratio = -1;
        }
        char msg[96];
        std::snprintf(msg, sizeof msg, "fixed-step BDF%d converges at order %d (error ratio per halving %.2f)", order,
                      order, ratio);
        check(std::abs(ratio - (order == 1 ? 2.0 : 4.0)) < 0.35, msg);
    }
    // 2. Adaptive: error follows the tolerance; fewer steps than fixed step.
    double e_prev = 1;
    bool monotone = true;
    for (double tol : {1e-3, 1e-5, 1e-7}) {
        IntegratorOptions o;
        o.rtol = tol; o.atol = tol;
        auto r = integrate(1, 3, lin_res, lin_jac, 0.0, y0.data(), nullptr, {1.0}, {}, o);
        double e = lin_err(r, 1.0);
        std::printf("adaptive rtol=%g: error %.2e, %d steps, %d rejected, %d factorizations\n", tol, e, r.stats.steps,
                    r.stats.rejected_error + r.stats.rejected_newton, r.stats.factorizations);
        // Local error control at order 2: global error exceeds tol by a factor that grows as tol
        // shrinks. Calibration: SUNDIALS IDA with max_order=2 on this problem gives 1x, 7x, 29x.
        monotone = monotone && r.status == Status::ok && e < e_prev && e < 200 * tol;
        e_prev = e;
    }
    check(monotone, "adaptive error decreases with tolerance (< 200x tol, as for order-2 IDA)");

    // 3. Stiff: lambda = 1e4 alongside 1; BDF must not need tiny steps after the transient.
    {
        auto res = [](double, const double* y, const double* yd, double* F) {
            F[0] = yd[0] + 1e4 * (y[0] - std::cos(0.0 + 0 * y[0]));   // relaxes fast to 1
            F[1] = yd[1] + y[1];
            F[2] = yd[2] + y[2];
        };
        auto jac = [](double, const double*, const double*, double a, BlockSystem& s) {
            s.B(0, 0, 0) = 1e4 + a; s.B(1, 0, 0) = 1 + a; s.B(2, 0, 0) = 1 + a;
        };
        std::vector<double> z0 = {0.0, 1.0, 1.0};
        IntegratorOptions o;
        o.rtol = 1e-4; o.atol = 1e-6;
        auto r = integrate(1, 3, res, jac, 0.0, z0.data(), nullptr, {10.0}, {}, o);
        std::printf("stiff: %d steps, y0(10)=%.6f, y1(10)=%.3e (exact %.3e)\n", r.stats.steps, r.y.back()[0],
                    r.y.back()[1], std::exp(-10.0));
        check(r.status == Status::ok && r.stats.steps < 300 && std::abs(r.y.back()[0] - 1) < 1e-4 &&
                  std::abs(r.y.back()[1] - std::exp(-10.0)) < 1e-4,
              "stiff problem: few steps, accurate");
    }
    // 4. Index-1 DAE with an algebraic variable: y' = -y + z, 0 = z - y^2; y(t) = 1/(1 + (1/y0 - 1) e^t).
    {
        auto res = [](double, const double* u, const double* ud, double* F) {
            for (int j = 0; j < 3; ++j) {
                const double y = u[2 * j], z = u[2 * j + 1];
                F[2 * j] = ud[2 * j] + y - z;
                F[2 * j + 1] = z - y * y;
            }
        };
        auto jac = [](double, const double* u, const double*, double a, BlockSystem& s) {
            for (int j = 0; j < 3; ++j) {
                s.B(j, 0, 0) = 1 + a; s.B(j, 0, 1) = -1;
                s.B(j, 1, 0) = -2 * u[2 * j]; s.B(j, 1, 1) = 1;
            }
        };
        std::vector<double> u0 = {0.5, 0.0, 0.5, 0.0, 0.5, 0.0};     // z inconsistent on purpose
        std::vector<char> alg = {0, 1, 0, 1, 0, 1};
        IntegratorOptions o;
        o.rtol = 1e-6; o.atol = 1e-8;
        auto r = integrate(2, 3, res, jac, 0.0, u0.data(), nullptr, {0.5, 2.0}, alg, o);
        const double ye = 1 / (1 + std::exp(2.0));
        std::printf("DAE: y(2)=%.8f exact %.8f, z(2)=%.8f (= y^2 %.8f), %d steps\n", r.y[1][0], ye, r.y[1][1],
                    ye * ye, r.stats.steps);
        check(r.status == Status::ok && r.t.size() == 2 && std::abs(r.t[0] - 0.5) < 1e-14,
              "output at the requested times");
        check(std::abs(r.y[1][0] - ye) < 2e-5 && std::abs(r.y[1][1] - ye * ye) < 2e-5,
              "index-1 DAE with inconsistent algebraic start: accurate");
    }
    // 5. Heat equation (coupled A/D blocks): u_t = u_xx, u = sin(pi x) e^{-pi^2 t}, Dirichlet ends.
    {
        const int nj = 41;
        const double h = 1.0 / (nj - 1);
        auto res = [&](double, const double* u, const double* ud, double* F) {
            F[0] = u[0]; F[nj - 1] = u[nj - 1];
            for (int j = 1; j < nj - 1; ++j) F[j] = ud[j] - (u[j + 1] - 2 * u[j] + u[j - 1]) / (h * h);
        };
        auto jac = [&](double, const double*, const double*, double a, BlockSystem& s) {
            s.B(0, 0, 0) = 1; s.B(nj - 1, 0, 0) = 1;
            for (int j = 1; j < nj - 1; ++j) {
                s.A(j, 0, 0) = s.D(j, 0, 0) = -1 / (h * h);
                s.B(j, 0, 0) = 2 / (h * h) + a;
            }
        };
        std::vector<double> u0(nj);
        for (int j = 0; j < nj; ++j) u0[j] = std::sin(pi * j * h);
        std::vector<char> alg(nj, 0);
        alg[0] = alg[nj - 1] = 1;
        IntegratorOptions o;
        o.rtol = 1e-6; o.atol = 1e-9;
        auto ra = integrate(1, nj, res, jac, 0.0, u0.data(), nullptr, {0.1}, alg, o);
        auto rf = integrate(1, nj, res, DaeJacobian{}, 0.0, u0.data(), nullptr, {0.1}, alg, o);
        IntegratorOptions onr = o;
        onr.jacobian_reuse = false;
        auto rn = integrate(1, nj, res, jac, 0.0, u0.data(), nullptr, {0.1}, alg, onr);
        // semi-discrete exact solution: eigenvalue of the 3-point Laplacian
        const double mu = 4 / (h * h) * std::pow(std::sin(pi * h / 2), 2);
        double e = 0, dfd = 0, dnr = 0;
        for (int j = 0; j < nj; ++j) {
            e = std::max(e, std::abs(ra.y.back()[j] - std::sin(pi * j * h) * std::exp(-mu * 0.1)));
            dfd = std::max(dfd, std::abs(ra.y.back()[j] - rf.y.back()[j]));
            dnr = std::max(dnr, std::abs(ra.y.back()[j] - rn.y.back()[j]));
        }
        std::printf("heat: error vs semi-discrete exact %.2e; %d steps; factorizations reuse %d vs no-reuse %d; "
                    "FD-Jacobian diff %.1e (%ld residual evals)\n",
                    e, ra.stats.steps, ra.stats.factorizations, rn.stats.factorizations, dfd,
                    rf.stats.residual_evaluations);
        check(ra.status == Status::ok && e < 1e-4, "heat equation: adaptive solution accurate");
        check(ra.stats.factorizations < rn.stats.factorizations && dnr < 1e-4,
              "Jacobian reuse: fewer factorizations, same answer within tolerance");
        check(rf.status == Status::ok && dfd < 1e-5, "finite-difference Jacobian path agrees");
    }
    // 6. Errors.
    {
        IntegratorOptions o;
        o.max_order = 3;
        auto r = integrate(1, 3, lin_res, lin_jac, 0.0, y0.data(), nullptr, {1.0}, {}, o);
        check(r.status == Status::invalid_argument, "max_order > 2 rejected");
        IntegratorOptions f;
        f.adaptive = false;          // dt missing
        r = integrate(1, 3, lin_res, lin_jac, 0.0, y0.data(), nullptr, {1.0}, {}, f);
        check(r.status == Status::invalid_argument, "fixed-step mode requires dt > 0");
        auto bad = [](double, const double*, const double*, double*) { throw std::runtime_error("boom"); };
        r = integrate(1, 3, bad, lin_jac, 0.0, y0.data(), nullptr, {1.0}, {}, IntegratorOptions{});
        check(r.status == Status::callback_error && r.callback_exception, "residual exception captured");
    }
    return failures ? 1 : 0;
}
