// T4: C++ kernel — backward error, agreement with the Fortran library, failure paths.
#include <bandsolver/band.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "bandsolver_f.h"

using namespace bandsolver;

static int failures = 0;
static void check(bool cond, const char* msg) {
    std::printf("%s %s\n", cond ? "PASS" : "FAIL", msg);
    if (!cond) ++failures;
}

static std::mt19937_64 rng(4242);
static double urand() { return std::uniform_real_distribution<double>(-1.0, 1.0)(rng); }

static BlockSystem random_system(int n, int nj, bool xy) {
    BlockSystem s(n, nj);
    for (auto* v : {&s.A(), &s.B(), &s.D(), &s.G()})
        for (auto& x : *v) x = urand();
    if (xy)
        for (auto* v : {&s.X(), &s.Y()})
            for (auto& x : *v) x = urand();
    for (int j = 0; j < nj; ++j)
        for (int i = 0; i < n; ++i) s.B(j, i, i) += 3.0 * n + 3.0;
    return s;
}

// ||K dc - G||_inf / (||K||_inf ||dc||_inf + ||G||_inf), applying K from the blocks.
static double backward_error(BlockSystem& s, const std::vector<double>& dc) {
    const int n = s.n(), nj = s.nj();
    double err = 0, knorm = 0, xnorm = 0, gnorm = 0;
    for (int j = 0; j < nj; ++j)
        for (int i = 0; i < n; ++i) {
            double r = -s.G(j, i), rs = 0;
            for (int k = 0; k < n; ++k) {
                r += s.B(j, i, k) * dc[j * n + k]; rs += std::abs(s.B(j, i, k));
                if (j > 0) { r += s.A(j, i, k) * dc[(j - 1) * n + k]; rs += std::abs(s.A(j, i, k)); }
                if (j < nj - 1) { r += s.D(j, i, k) * dc[(j + 1) * n + k]; rs += std::abs(s.D(j, i, k)); }
                if (j == 0) { r += s.X(i, k) * dc[2 * n + k]; rs += std::abs(s.X(i, k)); }
                if (j == nj - 1) { r += s.Y(i, k) * dc[(nj - 3) * n + k]; rs += std::abs(s.Y(i, k)); }
            }
            err = std::max(err, std::abs(r)); knorm = std::max(knorm, rs);
            xnorm = std::max(xnorm, std::abs(dc[j * n + i])); gnorm = std::max(gnorm, std::abs(s.G(j, i)));
        }
    return err / (knorm * xnorm + gnorm);
}

int main() {
    double berr = 0, fdiff = 0;
    bool all_ok = true;
    for (int n : {1, 2, 5, 12})
        for (int nj : {3, 4, 50, 500})
            for (bool xy : {false, true}) {
                BlockSystem s = random_system(n, nj, xy);
                for (Pivot p : {Pivot::partial, Pivot::legacy}) {
                    std::vector<double> dc(n * nj), df(n * nj);
                    SolveInfo info = solve(s.view(), dc.data(), p);
                    all_ok = all_ok && info.status == Status::ok;
                    berr = std::max(berr, backward_error(s, dc));
                    int fnode; double mrp;
                    int st = bandsolver_f_solve(n, nj, s.A().data(), s.B().data(), s.D().data(), s.G().data(),
                                                s.X().data(), s.Y().data(), static_cast<int>(p), df.data(), &fnode, &mrp);
                    all_ok = all_ok && st == 0;
                    double scale = 0, d = 0;
                    for (int i = 0; i < n * nj; ++i) {
                        scale = std::max(scale, std::abs(df[i]));
                        d = std::max(d, std::abs(dc[i] - df[i]));
                    }
                    fdiff = std::max(fdiff, d / (scale * 2.220446049250313e-16));
                }
            }
    std::printf("max backward error %.2e; max C++ vs Fortran diff %.1f eps*max|dc|\n", berr, fdiff);
    check(all_ok, "all sweep solves return ok (C++ and Fortran)");
    check(berr < 1e-13, "backward error sweep < 1e-13 (n up to 12, nj up to 500, with X/Y, both pivots)");
    check(fdiff <= 64, "C++ agrees with Fortran library within 64 eps (both pivot modes)");

    for (Pivot p : {Pivot::partial, Pivot::legacy})
        for (int node : {0, 3, 6}) {
            BlockSystem s = random_system(2, 7, false);
            for (int k = 0; k < 2; ++k) {
                s.A(node, 0, k) = s.A(node, 1, k) = 0;
                s.B(node, 1, k) = 2 * s.B(node, 0, k);
            }
            std::vector<double> dc(14);
            SolveInfo info = solve(s.view(), dc.data(), p);
            check(info.status == Status::singular && info.fail_node == node, "singular block reported at 0-based node");
            try {
                solve(s, p);
                check(false, "throwing solve throws on singular block");
            } catch (const Error& e) {
                check(e.status() == Status::singular && e.node() == node, "throwing solve throws Error with node");
            }
        }

    {
        BlockSystem s = random_system(2, 6, true);
        for (auto* v : {&s.A(), &s.D(), &s.X(), &s.Y()})
            for (auto& x : *v) x *= 0.1;
        std::fill(s.B().begin(), s.B().end(), 0.0);
        for (int j = 0; j < 6; ++j) s.B(j, 0, 1) = s.B(j, 1, 0) = 1;
        for (Pivot p : {Pivot::partial, Pivot::legacy}) {
            std::vector<double> dc(12);
            SolveInfo info = solve(s.view(), dc.data(), p);
            check(info.status == Status::ok && backward_error(s, dc) < 1e-14, "zero-diagonal blocks solve");
        }
    }
    // Nearly singular node-0 block (pivot 2^-52 relative), decoupled so every elimination
    // product is exact: the relative rule stops, the exact rule divides by the tiny pivot.
    for (Pivot p : {Pivot::partial, Pivot::legacy}) {
        BlockSystem s = random_system(2, 6, false);
        std::fill(s.X().begin(), s.X().end(), 0.0);
        std::fill(s.Y().begin(), s.Y().end(), 0.0);
        for (int i = 0; i < 2; ++i)
            for (int k = 0; k < 2; ++k) s.D(0, i, k) = s.A(1, i, k) = 0;
        s.B(0, 0, 0) = 1; s.B(0, 0, 1) = 0.5; s.B(0, 1, 0) = 2; s.B(0, 1, 1) = 1 + 2.220446049250313e-16;
        std::vector<double> dc(12);
        SolveInfo info = solve(s.view(), dc.data(), p);
        check(info.status == Status::singular && info.fail_node == 0, "nearly singular block: relative rule reports it");
        info = solve(s.view(), dc.data(), p, Singular::exact);
        check(info.status == Status::ok && backward_error(s, dc) < 1e-13, "nearly singular block: exact rule solves");
        std::vector<double> df(12);
        int fnode; double mrp;
        int st = bandsolver_f_solve_ex(2, 6, s.A().data(), s.B().data(), s.D().data(), s.G().data(), nullptr, nullptr,
                                       static_cast<int>(p), 0, 1, df.data(), &fnode, &mrp);
        check(st == 0 && df == dc, "nearly singular block: C++ and Fortran exact rule bit-identical");
        s.B(0, 1, 1) = 1;
        check(solve(s.view(), dc.data(), p, Singular::exact).status == Status::singular,
              "exactly singular block: exact rule reports it");
        check(solve(s.view(), dc.data(), p, static_cast<Singular>(5)).status == Status::invalid_argument,
              "unknown singular rule rejected");
    }

    {
        BlockSystem s = random_system(2, 5, false);
        s.G(1, 0) = std::nan("");
        std::vector<double> dc(10);
        check(solve(s.view(), dc.data()).status == Status::non_finite, "NaN input reported non-finite");
        SystemView v = s.view();
        v.nj = 2;
        check(solve(v, dc.data()).status == Status::invalid_argument, "nj < 3 rejected");
        check(solve(s.view(), dc.data(), static_cast<Pivot>(9)).status == Status::invalid_argument, "unknown pivot rejected");
        check(solve(s.view(), nullptr).status == Status::invalid_argument, "null output rejected");
        s.G(1, 0) = 0.5;
        v = s.view(); v.X = nullptr; v.Y = nullptr;
        check(solve(v, dc.data()).status == Status::ok, "null X/Y accepted");
        try {
            BlockSystem bad(0, 5);
            check(false, "BlockSystem rejects n < 1");
        } catch (const Error& e) {
            check(e.status() == Status::invalid_argument, "BlockSystem rejects n < 1");
        }
    }
    return failures ? 1 : 0;
}
