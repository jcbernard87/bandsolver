// P2: factor once, solve many. Factorization::solve must match the one-shot solve() on the
// same systems (incl. X/Y and nj=3), support repeated right-hand sides, and report failures.
#include <bandsolver/band.hpp>
#include <bandsolver/factor.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using namespace bandsolver;

static int failures = 0;
static void check(bool cond, const char* msg) {
    std::printf("%s %s\n", cond ? "PASS" : "FAIL", msg);
    if (!cond) ++failures;
}

static std::mt19937_64 rng(2718);
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

static double rel_diff(const std::vector<double>& a, const std::vector<double>& b) {
    double d = 0, s = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        d = std::max(d, std::abs(a[i] - b[i]));
        s = std::max(s, std::abs(b[i]));
    }
    return d / s;
}

int main() {
    double worst = 0;
    bool ok = true;
    for (int n : {1, 2, 5, 12})
        for (int nj : {3, 4, 50, 500})
            for (bool xy : {false, true}) {
                BlockSystem s = random_system(n, nj, xy);
                std::vector<double> ref = solve(s), dc(n * nj);
                Factorization f = factor(s.view());
                ok = ok && f.status() == Status::ok;
                SolveInfo info = f.solve(s.G().data(), dc.data());
                ok = ok && info.status == Status::ok;
                worst = std::max(worst, rel_diff(dc, ref));
            }
    std::printf("max |factor+solve - solve| / max|solve| = %.2e\n", worst);
    check(ok, "factor + solve succeed on the whole sweep");
    check(worst < 1e-13, "factor + solve matches one-shot solve (sweep incl. X/Y and nj = 3)");

    {  // several right-hand sides with one factorization
        BlockSystem s = random_system(4, 80, true);
        Factorization f = factor(s.view());
        double w = 0;
        for (int r = 0; r < 5; ++r) {
            for (auto& g : s.G()) g = urand();
            std::vector<double> ref = solve(s), dc(4 * 80);
            f.solve(s.G().data(), dc.data());
            w = std::max(w, rel_diff(dc, ref));
        }
        check(w < 1e-13, "one factorization serves repeated right-hand sides");
        std::vector<double> dc2 = f.solve(s.G());
        check(rel_diff(dc2, solve(s)) < 1e-13, "vector convenience solve");
    }
    {  // failures
        BlockSystem s = random_system(2, 9, false);
        for (int k = 0; k < 2; ++k) { s.A(5, 0, k) = s.A(5, 1, k) = 0; s.B(5, 1, k) = 2 * s.B(5, 0, k); }
        Factorization f = factor(s.view());
        check(f.status() == Status::singular && f.fail_node() == 5, "singular block reported at factor time (node 5)");
        std::vector<double> dc(18);
        check(f.solve(s.G().data(), dc.data()).status == Status::singular, "solving with a failed factorization reports it");
        BlockSystem t = random_system(2, 5, false);
        t.B(2, 0, 0) = std::nan("");
        check(factor(t.view()).status() == Status::non_finite, "non-finite blocks rejected");
        SystemView v = t.view(); v.nj = 2;
        check(factor(v).status() == Status::invalid_argument, "nj < 3 rejected");
        BlockSystem u = random_system(2, 5, false);
        Factorization fu = factor(u.view());
        std::vector<double> g(10, 1.0); g[3] = std::nan("");
        check(fu.solve(g.data(), dc.data()).status == Status::non_finite, "non-finite right-hand side rejected");
    }
    return failures ? 1 : 0;
}
