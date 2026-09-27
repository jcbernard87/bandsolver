// Native (no Python) timing of the BAND block solve, C++ core and Fortran core (C ABI).
// Usage: native_bench [--quick] > native.csv
// Output CSV: backend,n,nj,reps,median_s,min_s,backward_error
#include <bandsolver/band.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

#include "bandsolver_f.h"

using namespace bandsolver;
using Clock = std::chrono::steady_clock;

static BlockSystem make_system(int n, int nj, unsigned seed) {
    std::mt19937_64 g(seed);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    BlockSystem s(n, nj);
    for (auto* v : {&s.A(), &s.B(), &s.D(), &s.G(), &s.X(), &s.Y()})
        for (auto& x : *v) x = u(g);
    for (int j = 0; j < nj; ++j)
        for (int i = 0; i < n; ++i) s.B(j, i, i) += 3.0 * n + 3.0;
    return s;
}

static double backward_error(BlockSystem& s, const std::vector<double>& dc) {
    const int n = s.n(), nj = s.nj();
    double err = 0, kn = 0, xn = 0, gn = 0;
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
            err = std::max(err, std::abs(r)); kn = std::max(kn, rs);
            xn = std::max(xn, std::abs(dc[j * n + i])); gn = std::max(gn, std::abs(s.G(j, i)));
        }
    return err / (kn * xn + gn);
}

template <class F>
static void time_it(F&& f, double min_total, int min_reps, int& reps, double& med, double& mn) {
    std::vector<double> t;
    double total = 0;
    while (static_cast<int>(t.size()) < min_reps || total < min_total) {
        auto t0 = Clock::now();
        f();
        double dt = std::chrono::duration<double>(Clock::now() - t0).count();
        t.push_back(dt); total += dt;
        if (t.size() >= 100000) break;
    }
    std::sort(t.begin(), t.end());
    reps = static_cast<int>(t.size()); med = t[t.size() / 2]; mn = t.front();
}

int main(int argc, char** argv) {
    const bool quick = argc > 1 && std::strcmp(argv[1], "--quick") == 0;
    std::vector<int> ns = quick ? std::vector<int>{1, 3} : std::vector<int>{1, 3, 5, 10, 20, 30};
    std::vector<int> njs = quick ? std::vector<int>{25, 50} : std::vector<int>{25, 50, 100, 200, 500, 1000, 2000};
    const double min_total = quick ? 0.01 : 0.2;
    std::printf("backend,n,nj,reps,median_s,min_s,backward_error\n");
    for (int n : ns)
        for (int nj : njs) {
            BlockSystem s = make_system(n, nj, 1000u * n + nj);
            std::vector<double> dc(static_cast<std::size_t>(n) * nj);
            int reps; double med, mn;
            time_it([&] { solve(s.view(), dc.data(), Pivot::partial); }, min_total, 5, reps, med, mn);
            std::printf("band_cpp_native,%d,%d,%d,%.6e,%.6e,%.3e\n", n, nj, reps, med, mn, backward_error(s, dc));
            int fnode; double mrp;
            time_it([&] {
                bandsolver_f_solve(n, nj, s.A().data(), s.B().data(), s.D().data(), s.G().data(), s.X().data(),
                                   s.Y().data(), 0, dc.data(), &fnode, &mrp);
            }, min_total, 5, reps, med, mn);
            std::printf("band_fortran_native,%d,%d,%d,%.6e,%.6e,%.3e\n", n, nj, reps, med, mn, backward_error(s, dc));
            time_it([&] {
                bandsolver_f_solve_kernel(n, nj, s.A().data(), s.B().data(), s.D().data(), s.G().data(), s.X().data(),
                                          s.Y().data(), 0, BANDSOLVER_KERNEL_REFERENCE, dc.data(), &fnode, &mrp);
            }, min_total, 5, reps, med, mn);
            std::printf("band_fortran_reference_native,%d,%d,%d,%.6e,%.6e,%.3e\n", n, nj, reps, med, mn,
                        backward_error(s, dc));
            std::fflush(stdout);
        }
    return 0;
}
