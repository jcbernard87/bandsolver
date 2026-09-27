// F2: the Fortran finite-difference routines (via the C ABI) against the C++ core on the
// same residual/fill callbacks: identical Jacobians, Newton iterates, and check reports.
#include <bandsolver/fd.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "bandsolver_f.h"
#include "fd_problem.hpp"

using namespace bandsolver;

static int failures = 0;
static void check(bool cond, const char* msg) {
    std::printf("%s %s\n", cond ? "PASS" : "FAIL", msg);
    if (!cond) ++failures;
}

struct Ctx {
    const Problem* p;
    bool plant = false;
    bool fail = false;
};

extern "C" int c_residual(int, int, const double* c, double* F, void* ctx) {
    auto* x = static_cast<Ctx*>(ctx);
    if (x->fail) return 1;
    x->p->residual(c, F);
    return 0;
}

extern "C" int c_fill(int n, int nj, const double* c, double* A, double* B, double* D, double* G, double* X,
                      double* Y, void* ctx) {
    auto* x = static_cast<Ctx*>(ctx);
    BlockSystem s(n, nj);
    x->p->fill(c, s);
    if (x->plant) s.D(4, 1, 2) += 0.5;
    std::copy(s.A().begin(), s.A().end(), A);
    std::copy(s.B().begin(), s.B().end(), B);
    std::copy(s.D().begin(), s.D().end(), D);
    std::copy(s.G().begin(), s.G().end(), G);
    std::copy(s.X().begin(), s.X().end(), X);
    std::copy(s.Y().begin(), s.Y().end(), Y);
    return 0;
}

int main() {
    double worst = 0;
    bool counts = true;
    for (int n : {1, 3})
        for (int nj : {3, 4, 10}) {
            Problem p(n, nj, 77u + n + nj);
            Ctx ctx{&p};
            std::vector<double> c(n * nj);
            for (int i = 0; i < n * nj; ++i) c[i] = 0.3 + 0.01 * i;
            BlockSystem s(n, nj);
            long ec = fd_jacobian(n, nj, [&](const double* x, double* F) { p.residual(x, F); }, c.data(), s);
            std::vector<double> A(n * n * nj), B(A.size()), D(A.size()), G(n * nj), X(n * n), Y(n * n);
            long ef = 0;
            int st = bandsolver_f_fd_jacobian(n, nj, c_residual, &ctx, c.data(), nullptr, A.data(), B.data(),
                                              D.data(), G.data(), X.data(), Y.data(), &ef);
            counts = counts && st == 0 && ec == ef && ef == 3 * n + 1;
            auto diff = [&](const std::vector<double>& a, const std::vector<double>& b) {
                for (std::size_t i = 0; i < a.size(); ++i) worst = std::max(worst, std::abs(a[i] - b[i]));
            };
            diff(s.A(), A); diff(s.B(), B); diff(s.D(), D); diff(s.G(), G); diff(s.X(), X); diff(s.Y(), Y);
        }
    std::printf("max |C++ - Fortran| over all FD blocks: %.2e\n", worst);
    check(counts, "Fortran fd_jacobian returns OK with 3n+1 evaluations, like C++");
    check(worst <= 1e-12, "C++ and Fortran finite-difference Jacobians agree");

    {
        const int n = 3, nj = 40;
        Problem p(n, nj, 5);
        Ctx ctx{&p};
        std::vector<double> cc(n * nj), cf(n * nj);
        for (int i = 0; i < n * nj; ++i) cc[i] = cf[i] = 0.3 + 0.01 * i;
        NewtonResult rc = newton_fd(n, nj, [&](const double* x, double* F) { p.residual(x, F); }, cc.data());
        bandsolver_newton_options o;
        bandsolver_newton_result rf;
        bandsolver_f_default_options(&o);
        long ef = 0;
        int st = bandsolver_f_newton_fd(n, nj, c_residual, &ctx, cf.data(), &o, nullptr, &rf, nullptr, nullptr,
                                        nullptr, &ef);
        double d = 0;
        for (int i = 0; i < n * nj; ++i) d = std::max(d, std::abs(cc[i] - cf[i]));
        std::printf("newton_fd: C++ %d it / %ld evals, Fortran %d it / %ld evals, max diff %.2e\n", rc.iterations,
                    rc.residual_evaluations, rf.iterations, ef, d);
        check(st == 0 && rf.converged && rc.iterations == rf.iterations && rc.residual_evaluations == ef && d < 1e-13,
              "C++ and Fortran newton_fd agree");
        ctx.fail = true;
        st = bandsolver_f_newton_fd(n, nj, c_residual, &ctx, cf.data(), &o, nullptr, &rf, nullptr, nullptr, nullptr,
                                    nullptr);
        check(st == BANDSOLVER_CALLBACK_ERROR, "C residual error propagates (Fortran newton_fd)");
    }

    {  // Jacobian reuse through the C ABI (Fortran) vs the C++ driver, with a residual callback
        const int n = 3, nj = 40;
        Problem p(n, nj, 5);
        Ctx ctx{&p};
        std::vector<double> cc(n * nj), cf(n * nj);
        for (int i = 0; i < n * nj; ++i) cc[i] = cf[i] = 0.3 + 0.01 * i;
        NewtonOptions o;
        o.jacobian_reuse = true;
        NewtonResult rc = newton(n, nj, [&](const double* x, BlockSystem& s) { p.fill(x, s); }, cc.data(), o,
                                 [&](const double* x, double* F) { p.residual(x, F); });
        bandsolver_newton_options fo;
        bandsolver_newton_result fr;
        bandsolver_f_default_options(&fo);
        fo.jacobian_reuse = 1;
        int st = bandsolver_f_newton_ex(n, nj, c_fill, c_residual, &ctx, cf.data(), &fo, &fr, nullptr, nullptr, nullptr);
        double d = 0;
        for (int i = 0; i < n * nj; ++i) d = std::max(d, std::abs(cc[i] - cf[i]));
        std::printf("reuse: C++ %d it / %d factorizations / %ld residual calls; Fortran %d / %d / %d; max diff %.1e\n",
                    rc.iterations, rc.factorizations, rc.residual_evaluations, fr.iterations, fr.factorizations,
                    fr.residual_evaluations, d);
        check(st == 0 && fr.converged && rc.iterations == fr.iterations && rc.factorizations == fr.factorizations &&
                  rc.residual_evaluations == fr.residual_evaluations && d < 1e-12,
              "C++ and Fortran Jacobian-reuse drivers agree (iterations, factorizations, residual calls, root)");
        fo.jacobian_reuse = 0;
        std::vector<double> cg(n * nj);
        for (int i = 0; i < n * nj; ++i) cg[i] = 0.3 + 0.01 * i;
        st = bandsolver_f_newton_ex(n, nj, c_fill, c_residual, &ctx, cg.data(), &fo, &fr, nullptr, nullptr, nullptr);
        check(st == 0 && fr.factorizations == fr.iterations && fr.residual_evaluations == 0,
              "reuse off: full Newton through newton_ex (residual callback unused)");
    }

    {
        const int n = 3, nj = 10;
        Problem p(n, nj, 9);
        Ctx ctx{&p};
        ctx.plant = true;
        std::vector<double> c(n * nj);
        for (int i = 0; i < n * nj; ++i) c[i] = 0.3 + 0.01 * i;
        bandsolver_jacobian_check fchk;
        int st = bandsolver_f_check_jacobian(n, nj, c_fill, &ctx, c.data(), nullptr, &fchk);
        JacobianCheck cchk = check_jacobian(n, nj, [&](const double* x, BlockSystem& s) {
            p.fill(x, s);
            s.D(4, 1, 2) += 0.5;
        }, c.data());
        std::printf("check_jacobian D: Fortran node/row/col %d/%d/%d err %.3e; C++ %d/%d/%d err %.3e\n", fchk.D.node,
                    fchk.D.row, fchk.D.col, fchk.D.error, cchk.D.node, cchk.D.row, cchk.D.col, cchk.D.error);
        check(st == 0 && fchk.D.node - 1 == cchk.D.node && fchk.D.row - 1 == cchk.D.row &&
                  fchk.D.col - 1 == cchk.D.col && std::abs(fchk.D.error - cchk.D.error) <= 1e-12 * cchk.D.error,
              "C++ and Fortran check_jacobian report the same mismatch (Fortran 1-based)");
        check(std::abs(fchk.B.error - cchk.B.error) <= 1e-9, "clean blocks score identically");
    }
    return failures ? 1 : 0;
}
