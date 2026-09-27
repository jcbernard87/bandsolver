// Shared nonlinear test problem for the finite-difference Jacobian tests (n unknowns per
// node, full neighbour coupling, nonlinear X/Y endpoint terms, exact analytic Jacobian).
#ifndef BANDSOLVER_TEST_FD_PROBLEM_HPP
#define BANDSOLVER_TEST_FD_PROBLEM_HPP

#include <bandsolver/band.hpp>

#include <cmath>
#include <random>
#include <vector>

using bandsolver::BlockSystem;

// F_j(i) = sum_k Bm c_j + Am c_{j-1} + Dm c_{j+1} + 0.5 sin(c_j[i]) c_j[(i+1)%n] - b_i (j+1)
//          + [j=0] sum_k Xm c_2[k]^2 + [j=nj-1] sum_k Ym c_{nj-3}[k]^3
struct Problem {
    int n, nj;
    std::vector<double> Am, Bm, Dm, Xm, Ym, b;
    Problem(int n_, int nj_, unsigned seed) : n(n_), nj(nj_) {
        std::mt19937_64 g(seed);
        std::uniform_real_distribution<double> u(-1, 1);
        for (auto* v : {&Am, &Bm, &Dm, &Xm, &Ym}) {
            v->resize(n * n);
            for (auto& x : *v) x = 0.3 * u(g);
        }
        for (int i = 0; i < n; ++i) Bm[i * n + i] += 3.0 + n;
        b.resize(n);
        for (auto& x : b) x = u(g);
    }
    void residual(const double* c, double* F) const {
        for (int j = 0; j < nj; ++j)
            for (int i = 0; i < n; ++i) {
                const double* cj = c + j * n;
                double f = 0.5 * std::sin(cj[i]) * cj[(i + 1) % n] - b[i] * (j + 1);
                for (int k = 0; k < n; ++k) {
                    f += Bm[i * n + k] * cj[k];
                    if (j > 0) f += Am[i * n + k] * cj[k - n];
                    if (j < nj - 1) f += Dm[i * n + k] * cj[k + n];
                    if (j == 0) f += Xm[i * n + k] * c[2 * n + k] * c[2 * n + k];
                    if (j == nj - 1) f += Ym[i * n + k] * std::pow(c[(nj - 3) * n + k], 3);
                }
                F[j * n + i] = f;
            }
    }
    void fill(const double* c, BlockSystem& s) const {
        std::vector<double> F(n * nj);
        residual(c, F.data());
        for (int j = 0; j < nj; ++j)
            for (int i = 0; i < n; ++i) {
                const double* cj = c + j * n;
                s.G(j, i) = -F[j * n + i];
                for (int k = 0; k < n; ++k) {
                    s.B(j, i, k) = Bm[i * n + k];
                    if (j > 0) s.A(j, i, k) = Am[i * n + k];
                    if (j < nj - 1) s.D(j, i, k) = Dm[i * n + k];
                }
                s.B(j, i, i) += 0.5 * std::cos(cj[i]) * cj[(i + 1) % n];
                s.B(j, i, (i + 1) % n) += 0.5 * std::sin(cj[i]);
            }
        for (int i = 0; i < n; ++i)
            for (int k = 0; k < n; ++k) {
                s.X(i, k) = 2 * Xm[i * n + k] * c[2 * n + k];
                s.Y(i, k) = 3 * Ym[i * n + k] * std::pow(c[(nj - 3) * n + k], 2);
            }
    }
};

static std::vector<double> start(int n, int nj) {
    std::vector<double> c(n * nj);
    for (int i = 0; i < n * nj; ++i) c[i] = 0.3 + 0.01 * i;
    return c;
}

#endif
