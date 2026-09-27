#include "bandsolver/factor.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace bandsolver {

namespace {

bool finite(const double* p, std::size_t count) {
    for (std::size_t i = 0; i < count; ++i)
        if (!std::isfinite(p[i])) return false;
    return true;
}

// In-place LU with row partial pivoting of an n x n row-major block. Same singularity rule as
// the one-shot solver: a pivot with |p| <= n*eps*max|block| is singular.
bool lu_factor(int n, double* a, int* perm, double& rel) {
    const double eps = std::numeric_limits<double>::epsilon();
    double scale = 0;
    for (int i = 0; i < n * n; ++i) scale = std::max(scale, std::abs(a[i]));
    rel = std::numeric_limits<double>::max();
    for (int i = 0; i < n; ++i) perm[i] = i;
    if (scale == 0) return false;
    for (int k = 0; k < n; ++k) {
        int p = k;
        for (int i = k + 1; i < n; ++i)
            if (std::abs(a[i * n + k]) > std::abs(a[p * n + k])) p = i;
        if (std::abs(a[p * n + k]) <= n * eps * scale) return false;
        rel = std::min(rel, std::abs(a[p * n + k]) / scale);
        if (p != k) {
            std::swap_ranges(a + p * n, a + p * n + n, a + k * n);
            std::swap(perm[p], perm[k]);
        }
        const double inv = 1.0 / a[k * n + k];
        for (int i = k + 1; i < n; ++i) {
            const double l = a[i * n + k] * inv;
            a[i * n + k] = l;
            for (int c = k + 1; c < n; ++c) a[i * n + c] -= l * a[k * n + c];
        }
    }
    return true;
}

// Solve (LU) x = b for one vector: b is permuted into x, then forward and back substitution.
void lu_solve(int n, const double* lu, const int* perm, const double* b, double* x) {
    for (int i = 0; i < n; ++i) x[i] = b[perm[i]];
    for (int i = 1; i < n; ++i)
        for (int k = 0; k < i; ++k) x[i] -= lu[i * n + k] * x[k];
    for (int i = n - 1; i >= 0; --i) {
        for (int k = i + 1; k < n; ++k) x[i] -= lu[i * n + k] * x[k];
        x[i] /= lu[i * n + i];
    }
}

// out = -(LU)^{-1} M for an n x n row-major M (column by column).
void lu_solve_neg_matrix(int n, const double* lu, const int* perm, const double* M, double* out) {
    std::vector<double> col(n), x(n);
    for (int c = 0; c < n; ++c) {
        for (int i = 0; i < n; ++i) col[i] = M[i * n + c];
        lu_solve(n, lu, perm, col.data(), x.data());
        for (int i = 0; i < n; ++i) out[i * n + c] = -x[i];
    }
}

// C += A * B for n x n row-major blocks.
void gemm_add(int n, const double* A, const double* B, double* C) {
    for (int i = 0; i < n; ++i)
        for (int l = 0; l < n; ++l) {
            const double a = A[i * n + l];
            for (int k = 0; k < n; ++k) C[i * n + k] += a * B[l * n + k];
        }
}

}  // namespace

Factorization factor(const SystemView& s) {
    Factorization f;
    const int n = s.n, nj = s.nj;
    if (n < 1 || nj < 3 || !s.A || !s.B || !s.D) return f;   // status invalid_argument
    const std::size_t nn = static_cast<std::size_t>(n) * n, nb = nn * nj;
    if (!finite(s.A, nb) || !finite(s.B, nb) || !finite(s.D, nb) || (s.X && !finite(s.X, nn)) ||
        (s.Y && !finite(s.Y, nn))) {
        f.status_ = Status::non_finite;
        return f;
    }
    f.n_ = n; f.nj_ = nj;
    f.lu_.assign(nb, 0.0); f.E_.assign(nb, 0.0); f.Aeff_.assign(nb, 0.0);
    f.perm_.assign(static_cast<std::size_t>(n) * nj, 0);
    f.Xp_.assign(nn, 0.0); f.Y_.assign(nn, 0.0);
    if (s.Y) std::copy(s.Y, s.Y + nn, f.Y_.begin());
    auto blk = [nn](const double* base, int j) { return base + nn * j; };
    auto mblk = [nn](std::vector<double>& v, int j) { return v.data() + nn * j; };
    double rel = 0;
    f.min_rel_pivot_ = std::numeric_limits<double>::max();
    auto fail = [&](int node) {
        f.status_ = Status::singular;
        f.fail_node_ = node;
        return f;
    };

    // Node 0: LU(B_0); E_0 = -B_0^{-1} D_0; X' = -B_0^{-1} X.
    std::copy(blk(s.B, 0), blk(s.B, 0) + nn, mblk(f.lu_, 0));
    if (!lu_factor(n, mblk(f.lu_, 0), f.perm_.data(), rel)) return fail(0);
    f.min_rel_pivot_ = std::min(f.min_rel_pivot_, rel);
    lu_solve_neg_matrix(n, mblk(f.lu_, 0), f.perm_.data(), blk(s.D, 0), mblk(f.E_, 0));
    if (s.X) lu_solve_neg_matrix(n, mblk(f.lu_, 0), f.perm_.data(), s.X, f.Xp_.data());

    std::vector<double> Dm(nn);
    for (int j = 1; j < nj; ++j) {
        double* Aj = mblk(f.Aeff_, j);
        double* Bh = mblk(f.lu_, j);
        std::copy(blk(s.A, j), blk(s.A, j) + nn, Aj);
        std::copy(blk(s.B, j), blk(s.B, j) + nn, Bh);
        std::copy(blk(s.D, j), blk(s.D, j) + nn, Dm.begin());
        if (j == 1) gemm_add(n, Aj, f.Xp_.data(), Dm.data());            // D_1 += A_1 X'
        if (j == nj - 1) {
            gemm_add(n, f.Y_.data(), mblk(f.E_, j - 2), Aj);              // A += Y E_{nj-3}
            if (nj == 3) gemm_add(n, f.Y_.data(), f.Xp_.data(), Bh);      // B += Y X' (nj = 3)
        }
        gemm_add(n, Aj, mblk(f.E_, j - 1), Bh);                           // B^ = B + A' E_{j-1}
        if (!lu_factor(n, Bh, f.perm_.data() + static_cast<std::size_t>(n) * j, rel)) return fail(j);
        f.min_rel_pivot_ = std::min(f.min_rel_pivot_, rel);
        if (j < nj - 1)
            lu_solve_neg_matrix(n, Bh, f.perm_.data() + static_cast<std::size_t>(n) * j, Dm.data(), mblk(f.E_, j));
    }
    f.status_ = Status::ok;
    return f;
}

SolveInfo Factorization::solve(const double* G, double* dc) const {
    SolveInfo info;
    info.min_rel_pivot = min_rel_pivot_;
    if (status_ != Status::ok) {
        info.status = status_;
        info.fail_node = fail_node_;
        return info;
    }
    const int n = n_, nj = nj_;
    const std::size_t nn = static_cast<std::size_t>(n) * n, nv = static_cast<std::size_t>(n) * nj;
    if (!G || !dc) {
        info.status = Status::invalid_argument;
        return info;
    }
    if (!finite(G, nv)) {
        std::fill(dc, dc + nv, 0.0);
        info.status = Status::non_finite;
        return info;
    }
    // Forward sweep: e_j stored in dc.
    std::vector<double> g(n);
    lu_solve(n, lu_.data(), perm_.data(), G, dc);
    for (int j = 1; j < nj; ++j) {
        std::copy(G + static_cast<std::size_t>(j) * n, G + static_cast<std::size_t>(j) * n + n, g.begin());
        if (j == nj - 1) {
            const double* e2 = dc + static_cast<std::size_t>(j - 2) * n;
            for (int i = 0; i < n; ++i)
                for (int l = 0; l < n; ++l) g[i] -= Y_[i * n + l] * e2[l];
        }
        const double* A = Aeff_.data() + nn * j;
        const double* e1 = dc + static_cast<std::size_t>(j - 1) * n;
        for (int i = 0; i < n; ++i)
            for (int l = 0; l < n; ++l) g[i] -= A[i * n + l] * e1[l];
        lu_solve(n, lu_.data() + nn * j, perm_.data() + static_cast<std::size_t>(n) * j, g.data(),
                 dc + static_cast<std::size_t>(j) * n);
    }
    // Back substitution: dc_j = E_j dc_{j+1} + e_j, then dc_0 += X' dc_2.
    for (int j = nj - 2; j >= 0; --j) {
        const double* E = E_.data() + nn * j;
        double* x = dc + static_cast<std::size_t>(j) * n;
        const double* xn = x + n;
        for (int i = 0; i < n; ++i)
            for (int l = 0; l < n; ++l) x[i] += E[i * n + l] * xn[l];
    }
    for (int i = 0; i < n; ++i)
        for (int l = 0; l < n; ++l) dc[i] += Xp_[i * n + l] * dc[2 * n + l];
    info.status = finite(dc, nv) ? Status::ok : Status::non_finite;
    return info;
}

std::vector<double> Factorization::solve(const std::vector<double>& G) const {
    std::vector<double> dc(static_cast<std::size_t>(n_) * nj_);
    if (G.size() != dc.size()) throw Error(Status::invalid_argument, -1, "Factorization::solve: wrong G size");
    const SolveInfo info = solve(G.data(), dc.data());
    if (info.status != Status::ok)
        throw Error(info.status, info.fail_node,
                    std::string("Factorization::solve: ") + to_string(info.status));
    return dc;
}

}  // namespace bandsolver
