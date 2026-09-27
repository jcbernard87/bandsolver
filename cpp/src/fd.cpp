#include "bandsolver/fd.hpp"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace bandsolver {

long fd_jacobian(int n, int nj, const ResidualFunction& residual, const double* c, BlockSystem& sys,
                 const FdOptions& opts) {
    if (n < 1 || nj < 3 || !residual || !c || sys.n() != n || sys.nj() != nj || !(opts.rel_step > 0) ||
        !(opts.typical > 0))
        throw Error(Status::invalid_argument, -1, "fd_jacobian: invalid sizes or options");
    const std::size_t nv = static_cast<std::size_t>(n) * nj;
    sys.set_zero();
    std::vector<double> F0(nv), Fp(nv), cp(c, c + nv), h(nj);
    residual(c, F0.data());
    long evals = 1;
    for (std::size_t i = 0; i < nv; ++i) sys.G()[i] = -F0[i];

    auto dF = [&](int row, int i, int m) { return (Fp[static_cast<std::size_t>(row) * n + i] -
                                                   F0[static_cast<std::size_t>(row) * n + i]) / h[m]; };
    for (int r = 0; r < 3; ++r)
        for (int k = 0; k < n; ++k) {
            for (int m = r; m < nj; m += 3) {
                const double v = c[static_cast<std::size_t>(m) * n + k];
                const double step = opts.rel_step * std::max(std::abs(v), opts.typical) * (v < 0 ? -1.0 : 1.0);
                cp[static_cast<std::size_t>(m) * n + k] = v + step;
                h[m] = cp[static_cast<std::size_t>(m) * n + k] - v;  // the step actually taken
            }
            residual(cp.data(), Fp.data());
            ++evals;
            for (int m = r; m < nj; m += 3) {
                for (int i = 0; i < n; ++i) {
                    sys.B(m, i, k) = dF(m, i, m);
                    if (m >= 1) sys.D(m - 1, i, k) = dF(m - 1, i, m);
                    if (m + 1 <= nj - 1) sys.A(m + 1, i, k) = dF(m + 1, i, m);
                    if (m == 2) sys.X(i, k) = dF(0, i, m);
                    if (m == nj - 3) sys.Y(i, k) = dF(nj - 1, i, m);
                }
                cp[static_cast<std::size_t>(m) * n + k] = c[static_cast<std::size_t>(m) * n + k];
            }
        }
    return evals;
}

FillFunction fd_fill(int n, int nj, ResidualFunction residual, FdOptions opts, long* evaluations) {
    return [n, nj, residual = std::move(residual), opts, evaluations](const double* c, BlockSystem& s) {
        const long e = fd_jacobian(n, nj, residual, c, s, opts);
        if (evaluations) *evaluations += e;
    };
}

NewtonResult newton_fd(int n, int nj, const ResidualFunction& residual, double* c, const NewtonOptions& opts,
                       const FdOptions& fd) {
    long evals = 0;
    if (!residual) {
        NewtonResult r;
        r.status = Status::invalid_argument;
        return r;
    }
    ResidualFunction counted = [&](const double* x, double* F) {
        ++evals;
        residual(x, F);
    };
    NewtonResult r = newton(n, nj, fd_fill(n, nj, residual, fd, &evals), c, opts, counted);
    r.residual_evaluations = evals;
    return r;
}

double JacobianCheck::max_error() const noexcept {
    return std::max({A.error, B.error, D.error, X.error, Y.error});
}

JacobianCheck check_jacobian(int n, int nj, const FillFunction& fill, const double* c, const FdOptions& opts) {
    if (n < 1 || nj < 3 || !fill || !c) throw Error(Status::invalid_argument, -1, "check_jacobian: invalid arguments");
    BlockSystem user(n, nj), scratch(n, nj), fd(n, nj);
    fill(c, user);
    const std::size_t nv = static_cast<std::size_t>(n) * nj;
    auto residual = [&](const double* x, double* F) {
        scratch.set_zero();
        fill(x, scratch);
        for (std::size_t i = 0; i < nv; ++i) F[i] = -scratch.G()[i];
    };
    fd_jacobian(n, nj, residual, c, fd, opts);

    // Scale of each equation row (node j, row i): its largest Jacobian entry in any block.
    std::vector<double> rowscale(nv, 0.0);
    auto upd = [&](int j, int i, double a, double b) {
        double& r = rowscale[static_cast<std::size_t>(j) * n + i];
        r = std::max({r, std::abs(a), std::abs(b)});
    };
    for (int j = 0; j < nj; ++j)
        for (int i = 0; i < n; ++i)
            for (int k = 0; k < n; ++k) {
                if (j > 0) upd(j, i, user.A(j, i, k), fd.A(j, i, k));
                upd(j, i, user.B(j, i, k), fd.B(j, i, k));
                if (j < nj - 1) upd(j, i, user.D(j, i, k), fd.D(j, i, k));
                if (j == 0) upd(j, i, user.X(i, k), fd.X(i, k));
                if (j == nj - 1) upd(j, i, user.Y(i, k), fd.Y(i, k));
            }

    auto consider = [&](JacobianMismatch& m, int node, int i, int k, double u, double f) {
        const double floor = std::max(1e-3 * rowscale[static_cast<std::size_t>(node) * n + i], 1e-300);
        const double e = std::abs(u - f) / std::max({std::abs(u), std::abs(f), floor});
        if (m.node < 0 || e > m.error) m = JacobianMismatch{e, node, i, k, u, f};
    };
    JacobianCheck out;
    for (int j = 0; j < nj; ++j)
        for (int i = 0; i < n; ++i)
            for (int k = 0; k < n; ++k) {
                if (j > 0) consider(out.A, j, i, k, user.A(j, i, k), fd.A(j, i, k));
                consider(out.B, j, i, k, user.B(j, i, k), fd.B(j, i, k));
                if (j < nj - 1) consider(out.D, j, i, k, user.D(j, i, k), fd.D(j, i, k));
            }
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < n; ++k) {
            consider(out.X, 0, i, k, user.X(i, k), fd.X(i, k));
            consider(out.Y, nj - 1, i, k, user.Y(i, k), fd.Y(i, k));
        }
    return out;
}

}  // namespace bandsolver
