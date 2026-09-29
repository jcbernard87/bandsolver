#include "bandsolver/band.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace bandsolver {

const char* to_string(Status s) noexcept {
    switch (s) {
    case Status::ok: return "ok";
    case Status::singular: return "singular block";
    case Status::invalid_argument: return "invalid argument";
    case Status::not_converged: return "not converged";
    case Status::non_finite: return "non-finite value";
    case Status::callback_error: return "callback error";
    }
    return "unknown status";
}

Error::Error(Status status, int node, const std::string& what)
    : std::runtime_error(what), status_(status), node_(node) {}

BlockSystem::BlockSystem(int n, int nj) : n_(n), nj_(nj) {
    if (n < 1 || nj < 3) throw Error(Status::invalid_argument, -1, "BlockSystem requires n >= 1 and nj >= 3");
    const std::size_t nb = static_cast<std::size_t>(n) * n * nj;
    A_.assign(nb, 0.0);
    B_.assign(nb, 0.0);
    D_.assign(nb, 0.0);
    G_.assign(static_cast<std::size_t>(n) * nj, 0.0);
    X_.assign(static_cast<std::size_t>(n) * n, 0.0);
    Y_.assign(static_cast<std::size_t>(n) * n, 0.0);
}

void BlockSystem::set_zero() {
    for (auto* v : {&A_, &B_, &D_, &G_, &X_, &Y_}) std::fill(v->begin(), v->end(), 0.0);
}

SystemView BlockSystem::view() const noexcept {
    return SystemView{n_, nj_, A_.data(), B_.data(), D_.data(), G_.data(), X_.data(), Y_.data()};
}

namespace {

bool all_finite(const double* p, std::size_t count) {
    for (std::size_t i = 0; i < count; ++i)
        if (!std::isfinite(p[i])) return false;
    return true;
}

// Solve Bm * S = R in place; Bm is n x n, R is n x m, both row-major. Pivots with
// |p| <= tol*max|Bm| are reported singular: tol = n*eps (Singular::relative), or 0
// (Singular::exact, as the archival MATINV, which only stops on zero).
Status solve_partial(int n, int m, double* Bm, double* R, double& rel, double tol) {
    double bscale = 0;
    for (int i = 0; i < n * n; ++i) bscale = std::max(bscale, std::abs(Bm[i]));
    rel = std::numeric_limits<double>::max();
    if (bscale == 0) return Status::singular;
    for (int k = 0; k < n; ++k) {
        int p = k;
        for (int i = k + 1; i < n; ++i)
            if (std::abs(Bm[i * n + k]) > std::abs(Bm[p * n + k])) p = i;
        if (std::abs(Bm[p * n + k]) <= tol * bscale) return Status::singular;
        rel = std::min(rel, std::abs(Bm[p * n + k]) / bscale);
        if (p != k) {
            std::swap_ranges(Bm + p * n, Bm + p * n + n, Bm + k * n);
            std::swap_ranges(R + p * m, R + p * m + m, R + k * m);
        }
        double f = 1.0 / Bm[k * n + k];
        for (int c = k; c < n; ++c) Bm[k * n + c] *= f;
        for (int c = 0; c < m; ++c) R[k * m + c] *= f;
        for (int i = 0; i < n; ++i) {
            if (i == k) continue;
            f = Bm[i * n + k];
            if (f == 0) continue;
            for (int c = k; c < n; ++c) Bm[i * n + c] -= f * Bm[k * n + c];
            for (int c = 0; c < m; ++c) R[i * m + c] -= f * R[k * m + c];
        }
    }
    return Status::ok;
}

// Archival MATINV pivot heuristic (Appendix C.4), operation order preserved.
Status solve_legacy(int n, int m, double* Bm, double* R, double& rel, double tol) {
    const double bmax0 = static_cast<double>(1.1f);  // default-real literal in the source
    double bscale = 0;
    for (int i = 0; i < n * n; ++i) bscale = std::max(bscale, std::abs(Bm[i]));
    rel = std::numeric_limits<double>::max();
    std::vector<char> used(n, 0);
    int irow = 0, jcol = 0, jc = 0;
    for (int nn = 0; nn < n; ++nn) {
        double bmax = bmax0;
        bool found = false;
        for (int i = 0; i < n; ++i) {
            if (used[i]) continue;
            double bnext = 0, btry = 0;
            for (int j = 0; j < n; ++j) {
                if (used[j]) continue;
                const double a = std::abs(Bm[i * n + j]);
                if (a <= bnext) continue;
                bnext = a;
                if (bnext <= btry) continue;
                bnext = btry;
                btry = a;
                jc = j;
                found = true;
            }
            if (bnext >= bmax * btry) continue;
            bmax = bnext / btry;
            irow = i;
            jcol = jc;
        }
        if (!found) return Status::singular;
        used[jcol] = 1;
        if (jcol != irow) {
            std::swap_ranges(Bm + irow * n, Bm + irow * n + n, Bm + jcol * n);
            std::swap_ranges(R + irow * m, R + irow * m + m, R + jcol * m);
        }
        if (std::abs(Bm[jcol * n + jcol]) <= tol * bscale) return Status::singular;
        rel = std::min(rel, std::abs(Bm[jcol * n + jcol]) / bscale);
        double f = 1.0 / Bm[jcol * n + jcol];
        for (int j = 0; j < n; ++j) Bm[jcol * n + j] *= f;
        for (int k = 0; k < m; ++k) R[jcol * m + k] *= f;
        for (int i = 0; i < n; ++i) {
            if (i == jcol) continue;
            f = Bm[i * n + jcol];
            for (int j = 0; j < n; ++j) Bm[i * n + j] -= f * Bm[jcol * n + j];
            for (int k = 0; k < m; ++k) R[i * m + k] -= f * R[jcol * m + k];
        }
    }
    return Status::ok;
}

}  // namespace

SolveInfo solve(const SystemView& s, double* dc, Pivot pivot) {
    return solve(s, dc, pivot, Singular::relative);
}

SolveInfo solve(const SystemView& s, double* dc, Pivot pivot, Singular singular) {
    SolveInfo info;
    const int n = s.n, nj = s.nj;
    if (n < 1 || nj < 3 || !s.A || !s.B || !s.D || !s.G || !dc ||
        (pivot != Pivot::partial && pivot != Pivot::legacy) ||
        (singular != Singular::relative && singular != Singular::exact)) {
        info.status = Status::invalid_argument;
        return info;
    }
    const std::size_t nn = static_cast<std::size_t>(n) * n, nb = nn * nj, nv = static_cast<std::size_t>(n) * nj;
    std::fill(dc, dc + nv, 0.0);
    if (!all_finite(s.A, nb) || !all_finite(s.B, nb) || !all_finite(s.D, nb) || !all_finite(s.G, nv) ||
        (s.X && !all_finite(s.X, nn)) || (s.Y && !all_finite(s.Y, nn))) {
        info.status = Status::non_finite;
        return info;
    }

    {
        const int np1 = n + 1;
        auto blk = [n](const double* base, int j) { return base + static_cast<std::size_t>(j) * n * n; };
        // E[j] is n x (n+1): columns 0..n-1 hold E_j, column n holds e_j.
        std::vector<double> E(static_cast<std::size_t>(nj) * n * np1);
        auto Ej = [&](int j) { return E.data() + static_cast<std::size_t>(j) * n * np1; };
        std::vector<double> Xp(nn, 0.0), Yw(nn, 0.0), Am(nn), Bm(nn), Gm(n), R(static_cast<std::size_t>(n) * (2 * n + 1));
        if (s.X) std::copy(s.X, s.X + nn, Xp.begin());
        if (s.Y) std::copy(s.Y, s.Y + nn, Yw.begin());
        double rel = 0;
        const double tol = singular == Singular::exact ? 0.0 : n * std::numeric_limits<double>::epsilon();
        info.min_rel_pivot = std::numeric_limits<double>::max();

        auto block_solve = [&](int m, int node) {
            const Status st = pivot == Pivot::legacy ? solve_legacy(n, m, Bm.data(), R.data(), rel, tol)
                                                     : solve_partial(n, m, Bm.data(), R.data(), rel, tol);
            if (st != Status::ok) {
                info.status = st;
                info.fail_node = node;
                return false;
            }
            info.min_rel_pivot = std::min(info.min_rel_pivot, rel);
            return true;
        };

        // Node 0: B_0 [S_D | S_X | s_G] = [D_0 | X | G_0], R is n x (2n+1).
        {
            const int m = 2 * n + 1;
            std::copy(blk(s.B, 0), blk(s.B, 0) + nn, Bm.begin());
            for (int i = 0; i < n; ++i) {
                for (int k = 0; k < n; ++k) {
                    R[i * m + k] = blk(s.D, 0)[i * n + k];
                    R[i * m + n + k] = Xp[i * n + k];
                }
                R[i * m + 2 * n] = s.G[i];
            }
            if (!block_solve(m, 0)) {
                std::fill(dc, dc + nv, 0.0);
                return info;
            }
            double* E0 = Ej(0);
            for (int k = 0; k < n; ++k) {
                E0[k * np1 + n] = R[k * m + 2 * n];
                for (int l = 0; l < n; ++l) {
                    E0[k * np1 + l] = -R[k * m + l];
                    Xp[k * n + l] = -R[k * m + n + l];
                }
            }
        }

        for (int j = 1; j < nj; ++j) {
            const int m = np1;
            std::copy(blk(s.A, j), blk(s.A, j) + nn, Am.begin());
            std::copy(blk(s.B, j), blk(s.B, j) + nn, Bm.begin());
            std::copy(s.G + static_cast<std::size_t>(j) * n, s.G + static_cast<std::size_t>(j) * n + n, Gm.begin());
            for (int i = 0; i < n; ++i)
                for (int k = 0; k < n; ++k) R[i * m + k] = blk(s.D, j)[i * n + k];
            if (j == 1) {  // D_1 += A_1 X'
                for (int i = 0; i < n; ++i)
                    for (int k = 0; k < n; ++k)
                        for (int l = 0; l < n; ++l) R[i * m + k] += Am[i * n + l] * Xp[l * n + k];
            }
            if (j == nj - 1) {  // eliminate Y dc_{nj-3}
                const double* E2 = Ej(j - 2);
                for (int i = 0; i < n; ++i)
                    for (int l = 0; l < n; ++l) {
                        Gm[i] -= Yw[i * n + l] * E2[l * np1 + n];
                        for (int q = 0; q < n; ++q) Am[i * n + l] += Yw[i * n + q] * E2[q * np1 + l];
                    }
                if (nj == 3) {  // node 0's X' dc_2 lands on the last node (absent in the archival kernel)
                    for (int i = 0; i < n; ++i)
                        for (int k = 0; k < n; ++k)
                            for (int l = 0; l < n; ++l) Bm[i * n + k] += Yw[i * n + l] * Xp[l * n + k];
                }
            }
            const double* E1 = Ej(j - 1);
            for (int i = 0; i < n; ++i) {
                R[i * m + n] = -Gm[i];
                for (int l = 0; l < n; ++l) {
                    R[i * m + n] += Am[i * n + l] * E1[l * np1 + n];
                    for (int k = 0; k < n; ++k) Bm[i * n + k] += Am[i * n + l] * E1[l * np1 + k];
                }
            }
            if (!block_solve(m, j)) {
                std::fill(dc, dc + nv, 0.0);
                return info;
            }
            double* Ec = Ej(j);
            for (int k = 0; k < n; ++k)
                for (int q = 0; q < np1; ++q) Ec[k * np1 + q] = -R[k * m + q];
        }

        // Back substitution.
        for (int k = 0; k < n; ++k) dc[static_cast<std::size_t>(nj - 1) * n + k] = Ej(nj - 1)[k * np1 + n];
        for (int j = nj - 2; j >= 0; --j) {
            const double* Ec = Ej(j);
            double* x = dc + static_cast<std::size_t>(j) * n;
            const double* xn = x + n;
            for (int k = 0; k < n; ++k) {
                x[k] = Ec[k * np1 + n];
                for (int l = 0; l < n; ++l) x[k] += Ec[k * np1 + l] * xn[l];
            }
        }
        for (int l = 0; l < n; ++l)
            for (int k = 0; k < n; ++k) dc[k] += Xp[k * n + l] * dc[2 * n + l];
    }

    info.status = all_finite(dc, nv) ? Status::ok : Status::non_finite;
    return info;
}

std::vector<double> solve(const BlockSystem& sys, Pivot pivot) {
    return solve(sys, pivot, Singular::relative);
}

std::vector<double> solve(const BlockSystem& sys, Pivot pivot, Singular singular) {
    std::vector<double> dc(static_cast<std::size_t>(sys.n()) * sys.nj());
    const SolveInfo info = solve(sys.view(), dc.data(), pivot, singular);
    if (info.status != Status::ok) {
        std::string msg = std::string("bandsolver::solve: ") + to_string(info.status);
        if (info.fail_node >= 0) msg += " at node " + std::to_string(info.fail_node);
        throw Error(info.status, info.fail_node, msg);
    }
    return dc;
}

}  // namespace bandsolver
