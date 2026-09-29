// pybind11 extension `bandsolver._core`: exposes the C++ core and the Fortran library
// (through its C ABI) with identical call signatures. Shape validation and exception
// mapping live in the Python package; this layer checks only what memory safety needs.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <exception>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <bandsolver/band.hpp>
#include <bandsolver/newton.hpp>
#include <bandsolver/fd.hpp>
#include <bandsolver/factor.hpp>
#include <bandsolver/integrate.hpp>

#include "bandsolver_f.h"

namespace py = pybind11;
using Arr = py::array_t<double, py::array::c_style | py::array::forcecast>;

namespace {

void require_size(const Arr& a, py::ssize_t expected, const char* name) {
    if (a.size() != expected)
        throw std::invalid_argument(std::string(name) + " has " + std::to_string(a.size()) + " values, expected " +
                                    std::to_string(expected));
}

const double* opt_ptr(const std::optional<Arr>& a, py::ssize_t expected, const char* name) {
    if (!a) return nullptr;
    require_size(*a, expected, name);
    return a->data();
}

py::tuple solve(int n, int nj, const Arr& A, const Arr& B, const Arr& D, const Arr& G,
                const std::optional<Arr>& X, const std::optional<Arr>& Y, int pivot, const std::string& backend,
                int kernel, int singular) {
    if (n < 1 || nj < 3) throw std::invalid_argument("require n >= 1 and nj >= 3");
    const py::ssize_t nb = static_cast<py::ssize_t>(n) * n * nj, nv = static_cast<py::ssize_t>(n) * nj;
    require_size(A, nb, "A"); require_size(B, nb, "B"); require_size(D, nb, "D"); require_size(G, nv, "G");
    const double* xp = opt_ptr(X, static_cast<py::ssize_t>(n) * n, "X");
    const double* yp = opt_ptr(Y, static_cast<py::ssize_t>(n) * n, "Y");
    Arr dc({static_cast<py::ssize_t>(nj), static_cast<py::ssize_t>(n)});
    int status, fail_node;
    double min_rel_pivot;
    if (backend == "cpp") {
        bandsolver::SystemView v{n, nj, A.data(), B.data(), D.data(), G.data(), xp, yp};
        bandsolver::SolveInfo info;
        {
            py::gil_scoped_release release;
            info = bandsolver::solve(v, dc.mutable_data(), static_cast<bandsolver::Pivot>(pivot),
                                     static_cast<bandsolver::Singular>(singular));
        }
        status = static_cast<int>(info.status);
        fail_node = info.fail_node;
        min_rel_pivot = info.min_rel_pivot;
    } else if (backend == "fortran") {
        {
            py::gil_scoped_release release;
            status = bandsolver_f_solve_ex(n, nj, A.data(), B.data(), D.data(), G.data(), xp, yp, pivot, kernel,
                                           singular, dc.mutable_data(), &fail_node, &min_rel_pivot);
        }
        fail_node = fail_node > 0 ? fail_node - 1 : -1;  // 1-based -> 0-based
    } else {
        throw std::invalid_argument("backend must be 'cpp' or 'fortran'");
    }
    return py::make_tuple(dc, status, fail_node, min_rel_pivot);
}

// Calls the Python fill(c) and copies its (A, B, D, G[, X, Y]) into the output buffers.
void call_fill(const py::function& fill, int n, int nj, const double* c, double* A, double* B, double* D,
               double* G, double* X, double* Y) {
    Arr cv({static_cast<py::ssize_t>(nj), static_cast<py::ssize_t>(n)});
    std::copy(c, c + static_cast<std::size_t>(n) * nj, cv.mutable_data());
    py::object out = fill(cv);
    auto t = out.cast<py::sequence>();
    if (t.size() != 4 && t.size() != 6)
        throw std::invalid_argument("fill(c) must return (A, B, D, G) or (A, B, D, G, X, Y)");
    const py::ssize_t nb = static_cast<py::ssize_t>(n) * n * nj, nv = static_cast<py::ssize_t>(n) * nj,
                      nn = static_cast<py::ssize_t>(n) * n;
    const char* names[6] = {"A", "B", "D", "G", "X", "Y"};
    double* dst[6] = {A, B, D, G, X, Y};
    const py::ssize_t sizes[6] = {nb, nb, nb, nv, nn, nn};
    for (std::size_t i = 0; i < t.size(); ++i) {
        if (t[i].is_none() && i >= 4) continue;
        Arr a = Arr::ensure(t[i]);
        if (!a) throw std::invalid_argument(std::string("fill output ") + names[i] + " is not numeric");
        require_size(a, sizes[i], names[i]);
        std::copy(a.data(), a.data() + sizes[i], dst[i]);
    }
}

struct FortranCtx {
    const py::function* fill;
    std::exception_ptr error;
};

extern "C" int fortran_trampoline(int n, int nj, const double* c, double* A, double* B, double* D, double* G,
                                  double* X, double* Y, void* ctx) {
    auto* f = static_cast<FortranCtx*>(ctx);
    try {  // exceptions must not unwind through Fortran frames
        call_fill(*f->fill, n, nj, c, A, B, D, G, X, Y);
        return 0;
    } catch (...) {
        f->error = std::current_exception();
        return 1;
    }
}

struct FortranPairCtx {
    const py::function* fill;
    const py::function* residual;
    std::exception_ptr error;
};

extern "C" int pair_fill(int n, int nj, const double* c, double* A, double* B, double* D, double* G, double* X,
                         double* Y, void* ctx) {
    auto* f = static_cast<FortranPairCtx*>(ctx);
    try { call_fill(*f->fill, n, nj, c, A, B, D, G, X, Y); return 0; }
    catch (...) { f->error = std::current_exception(); return 1; }
}

void call_residual(const py::function& residual, int n, int nj, const double* c, double* F);

extern "C" int pair_residual(int n, int nj, const double* c, double* F, void* ctx) {
    auto* f = static_cast<FortranPairCtx*>(ctx);
    try { call_residual(*f->residual, n, nj, c, F); return 0; }
    catch (...) { f->error = std::current_exception(); return 1; }
}

py::dict newton(int n, int nj, const py::function& fill, const Arr& c0, double rtol, double atol, double damping,
                int max_iter, int pivot, bool require_convergence, const std::string& backend, int kernel,
                bool jacobian_reuse, int reuse_max_iter, double reuse_contraction,
                const std::optional<py::function>& residual) {
    if (n < 1 || nj < 3) throw std::invalid_argument("require n >= 1 and nj >= 3");
    if (max_iter < 1) throw std::invalid_argument("max_iter must be >= 1");
    require_size(c0, static_cast<py::ssize_t>(n) * nj, "c0");
    Arr c({static_cast<py::ssize_t>(nj), static_cast<py::ssize_t>(n)});
    std::copy(c0.data(), c0.data() + c0.size(), c.mutable_data());
    py::dict r;
    if (backend == "cpp") {
        bandsolver::NewtonOptions o;
        o.rtol = rtol; o.atol = atol; o.damping = damping; o.max_iter = max_iter;
        o.pivot = static_cast<bandsolver::Pivot>(pivot); o.require_convergence = require_convergence;
        o.jacobian_reuse = jacobian_reuse; o.reuse_max_iter = reuse_max_iter; o.reuse_contraction = reuse_contraction;
        bandsolver::ResidualFunction rf;
        if (residual) rf = [&](const double* x, double* F) { call_residual(*residual, n, nj, x, F); };
        auto res = bandsolver::newton(
            n, nj,
            [&](const double* cc, bandsolver::BlockSystem& s) {
                call_fill(fill, n, nj, cc, s.A().data(), s.B().data(), s.D().data(), s.G().data(), s.X().data(),
                          s.Y().data());
            },
            c.mutable_data(), o, rf);
        if (res.callback_exception) std::rethrow_exception(res.callback_exception);
        r["status"] = static_cast<int>(res.status);
        r["iterations"] = res.iterations;
        r["converged"] = res.converged;
        r["fail_node"] = res.fail_node;
        r["update_norm"] = res.update_norm;
        r["step_norm"] = res.step_norm;
        r["residual_norm"] = res.residual_norm;
        r["jacobian_evaluations"] = res.jacobian_evaluations;
        r["factorizations"] = res.factorizations;
        r["residual_evaluations"] = res.residual_evaluations;
    } else if (backend == "fortran") {
        bandsolver_newton_options o{rtol, atol, damping, max_iter, pivot, require_convergence ? 1 : 0, kernel,
                                    jacobian_reuse ? 1 : 0, reuse_max_iter, reuse_contraction};
        bandsolver_newton_result res;
        std::vector<double> hu(max_iter), hs(max_iter), hr(max_iter);
        const py::function* rp = residual ? &*residual : nullptr;
        FortranPairCtx ctx{&fill, rp, nullptr};
        bandsolver_f_newton_ex(n, nj, pair_fill, rp ? pair_residual : nullptr, &ctx, c.mutable_data(), &o, &res,
                               hu.data(), hs.data(), hr.data());
        if (ctx.error) std::rethrow_exception(ctx.error);
        // Histories: residual is recorded before each solve, step/update after a successful one.
        const int nres = std::min(res.iterations, max_iter);
        const int nstep = (res.status == BANDSOLVER_OK || res.status == BANDSOLVER_NOT_CONVERGED) ? nres
                                                                                                    : std::max(nres - 1, 0);
        r["status"] = res.status;
        r["iterations"] = res.iterations;
        r["converged"] = res.converged != 0;
        r["fail_node"] = res.fail_node > 0 ? res.fail_node - 1 : -1;
        r["update_norm"] = std::vector<double>(hu.begin(), hu.begin() + nstep);
        r["step_norm"] = std::vector<double>(hs.begin(), hs.begin() + nstep);
        r["residual_norm"] = std::vector<double>(hr.begin(), hr.begin() + nres);
        r["jacobian_evaluations"] = res.jacobian_evaluations;
        r["factorizations"] = res.factorizations;
        r["residual_evaluations"] = res.residual_evaluations;
    } else {
        throw std::invalid_argument("backend must be 'cpp' or 'fortran'");
    }
    r["c"] = c;
    return r;
}

// Calls the Python residual(c) and copies F into the output buffer.
void call_residual(const py::function& residual, int n, int nj, const double* c, double* F) {
    Arr cv({static_cast<py::ssize_t>(nj), static_cast<py::ssize_t>(n)});
    std::copy(c, c + static_cast<std::size_t>(n) * nj, cv.mutable_data());
    Arr f = Arr::ensure(residual(cv));
    if (!f) throw std::invalid_argument("residual(c) must return a numeric array");
    require_size(f, static_cast<py::ssize_t>(n) * nj, "F");
    std::copy(f.data(), f.data() + f.size(), F);
}

struct FortranResidualCtx {
    const py::function* residual;
    std::exception_ptr error;
};

extern "C" int fortran_residual_trampoline(int n, int nj, const double* c, double* F, void* ctx) {
    auto* f = static_cast<FortranResidualCtx*>(ctx);
    try {  // exceptions must not unwind through Fortran frames
        call_residual(*f->residual, n, nj, c, F);
        return 0;
    } catch (...) {
        f->error = std::current_exception();
        return 1;
    }
}

bandsolver::FdOptions fd_opts(double rel_step, double typical) {
    if (!(rel_step > 0) || !(typical > 0)) throw std::invalid_argument("rel_step and typical must be > 0");
    bandsolver::FdOptions o;
    o.rel_step = rel_step;
    o.typical = typical;
    return o;
}

Arr blocks(int n, int nj) {
    return Arr({static_cast<py::ssize_t>(nj), static_cast<py::ssize_t>(n), static_cast<py::ssize_t>(n)});
}

py::tuple fd_jacobian(int n, int nj, const py::function& residual, const Arr& c, double rel_step, double typical,
                      const std::string& backend) {
    if (n < 1 || nj < 3) throw std::invalid_argument("require n >= 1 and nj >= 3");
    require_size(c, static_cast<py::ssize_t>(n) * nj, "c");
    const bandsolver::FdOptions o = fd_opts(rel_step, typical);
    Arr A = blocks(n, nj), B = blocks(n, nj), D = blocks(n, nj);
    Arr G({static_cast<py::ssize_t>(nj), static_cast<py::ssize_t>(n)});
    Arr X({static_cast<py::ssize_t>(n), static_cast<py::ssize_t>(n)}), Y({static_cast<py::ssize_t>(n), static_cast<py::ssize_t>(n)});
    long evals = 0;
    if (backend == "cpp") {
        bandsolver::BlockSystem s(n, nj);
        evals = bandsolver::fd_jacobian(n, nj, [&](const double* x, double* F) { call_residual(residual, n, nj, x, F); },
                                        c.data(), s, o);
        std::copy(s.A().begin(), s.A().end(), A.mutable_data());
        std::copy(s.B().begin(), s.B().end(), B.mutable_data());
        std::copy(s.D().begin(), s.D().end(), D.mutable_data());
        std::copy(s.G().begin(), s.G().end(), G.mutable_data());
        std::copy(s.X().begin(), s.X().end(), X.mutable_data());
        std::copy(s.Y().begin(), s.Y().end(), Y.mutable_data());
    } else if (backend == "fortran") {
        FortranResidualCtx ctx{&residual, nullptr};
        bandsolver_fd_options fo{o.rel_step, o.typical};
        int st = bandsolver_f_fd_jacobian(n, nj, fortran_residual_trampoline, &ctx, c.data(), &fo, A.mutable_data(),
                                          B.mutable_data(), D.mutable_data(), G.mutable_data(), X.mutable_data(),
                                          Y.mutable_data(), &evals);
        if (ctx.error) std::rethrow_exception(ctx.error);
        if (st != 0) throw std::runtime_error("bandsolver_f_fd_jacobian failed with status " + std::to_string(st));
    } else {
        throw std::invalid_argument("backend must be 'cpp' or 'fortran'");
    }
    return py::make_tuple(A, B, D, G, X, Y, evals);
}

py::dict newton_fd(int n, int nj, const py::function& residual, const Arr& c0, double rtol, double atol,
                   double damping, int max_iter, int pivot, bool require_convergence, double rel_step, double typical,
                   const std::string& backend, int kernel, bool jacobian_reuse, int reuse_max_iter,
                   double reuse_contraction) {
    if (n < 1 || nj < 3) throw std::invalid_argument("require n >= 1 and nj >= 3");
    if (max_iter < 1) throw std::invalid_argument("max_iter must be >= 1");
    require_size(c0, static_cast<py::ssize_t>(n) * nj, "c0");
    const bandsolver::FdOptions o = fd_opts(rel_step, typical);
    Arr c({static_cast<py::ssize_t>(nj), static_cast<py::ssize_t>(n)});
    std::copy(c0.data(), c0.data() + c0.size(), c.mutable_data());
    py::dict r;
    if (backend == "cpp") {
        bandsolver::NewtonOptions no;
        no.rtol = rtol; no.atol = atol; no.damping = damping; no.max_iter = max_iter;
        no.pivot = static_cast<bandsolver::Pivot>(pivot); no.require_convergence = require_convergence;
        no.jacobian_reuse = jacobian_reuse; no.reuse_max_iter = reuse_max_iter; no.reuse_contraction = reuse_contraction;
        auto res = bandsolver::newton_fd(
            n, nj, [&](const double* x, double* F) { call_residual(residual, n, nj, x, F); }, c.mutable_data(), no, o);
        if (res.callback_exception) std::rethrow_exception(res.callback_exception);
        r["status"] = static_cast<int>(res.status);
        r["iterations"] = res.iterations;
        r["converged"] = res.converged;
        r["fail_node"] = res.fail_node;
        r["update_norm"] = res.update_norm;
        r["step_norm"] = res.step_norm;
        r["residual_norm"] = res.residual_norm;
        r["residual_evaluations"] = res.residual_evaluations;
        r["jacobian_evaluations"] = res.jacobian_evaluations;
        r["factorizations"] = res.factorizations;
    } else if (backend == "fortran") {
        bandsolver_newton_options no{rtol, atol, damping, max_iter, pivot, require_convergence ? 1 : 0, kernel,
                                     jacobian_reuse ? 1 : 0, reuse_max_iter, reuse_contraction};
        bandsolver_fd_options fo{o.rel_step, o.typical};
        bandsolver_newton_result res;
        std::vector<double> hu(max_iter), hs(max_iter), hr(max_iter);
        long evals = 0;
        FortranResidualCtx ctx{&residual, nullptr};
        bandsolver_f_newton_fd(n, nj, fortran_residual_trampoline, &ctx, c.mutable_data(), &no, &fo, &res, hu.data(),
                               hs.data(), hr.data(), &evals);
        if (ctx.error) std::rethrow_exception(ctx.error);
        const int nres = std::min(res.iterations, max_iter);
        const int nstep = (res.status == BANDSOLVER_OK || res.status == BANDSOLVER_NOT_CONVERGED) ? nres
                                                                                                    : std::max(nres - 1, 0);
        r["status"] = res.status;
        r["iterations"] = res.iterations;
        r["converged"] = res.converged != 0;
        r["fail_node"] = res.fail_node > 0 ? res.fail_node - 1 : -1;
        r["update_norm"] = std::vector<double>(hu.begin(), hu.begin() + nstep);
        r["step_norm"] = std::vector<double>(hs.begin(), hs.begin() + nstep);
        r["residual_norm"] = std::vector<double>(hr.begin(), hr.begin() + nres);
        r["residual_evaluations"] = evals;
        r["jacobian_evaluations"] = res.jacobian_evaluations;
        r["factorizations"] = res.factorizations;
    } else {
        throw std::invalid_argument("backend must be 'cpp' or 'fortran'");
    }
    r["c"] = c;
    return r;
}

py::tuple mismatch(double error, int node, int row, int col, double user, double fd) {
    return py::make_tuple(error, node, row, col, user, fd);
}

// Returns {block: (error, node, row, col, user, fd)} with 0-based indices (-1 if empty).
py::dict check_jacobian(int n, int nj, const py::function& fill, const Arr& c, double rel_step, double typical,
                        const std::string& backend) {
    if (n < 1 || nj < 3) throw std::invalid_argument("require n >= 1 and nj >= 3");
    require_size(c, static_cast<py::ssize_t>(n) * nj, "c");
    const bandsolver::FdOptions o = fd_opts(rel_step, typical);
    py::dict r;
    if (backend == "cpp") {
        auto chk = bandsolver::check_jacobian(
            n, nj,
            [&](const double* x, bandsolver::BlockSystem& s) {
                call_fill(fill, n, nj, x, s.A().data(), s.B().data(), s.D().data(), s.G().data(), s.X().data(),
                          s.Y().data());
            },
            c.data(), o);
        const char* names[5] = {"A", "B", "D", "X", "Y"};
        const bandsolver::JacobianMismatch* ms[5] = {&chk.A, &chk.B, &chk.D, &chk.X, &chk.Y};
        for (int b = 0; b < 5; ++b)
            r[names[b]] = mismatch(ms[b]->error, ms[b]->node, ms[b]->row, ms[b]->col, ms[b]->user, ms[b]->fd);
    } else if (backend == "fortran") {
        FortranCtx ctx{&fill, nullptr};
        bandsolver_fd_options fo{o.rel_step, o.typical};
        bandsolver_jacobian_check chk;
        int st = bandsolver_f_check_jacobian(n, nj, fortran_trampoline, &ctx, c.data(), &fo, &chk);
        if (ctx.error) std::rethrow_exception(ctx.error);
        if (st != 0) throw std::runtime_error("bandsolver_f_check_jacobian failed with status " + std::to_string(st));
        const char* names[5] = {"A", "B", "D", "X", "Y"};
        const bandsolver_jacobian_mismatch* ms[5] = {&chk.A, &chk.B, &chk.D, &chk.X, &chk.Y};
        for (int b = 0; b < 5; ++b)  // 1-based -> 0-based (0 = no entries -> -1)
            r[names[b]] = mismatch(ms[b]->error, ms[b]->node - 1, ms[b]->row - 1, ms[b]->col - 1, ms[b]->user, ms[b]->fd);
    } else {
        throw std::invalid_argument("backend must be 'cpp' or 'fortran'");
    }
    return r;
}

// ---- DAE integration (C++ core) -------------------------------------------------------
Arr state_array(int n, int nj, const double* src) {
    Arr a({static_cast<py::ssize_t>(nj), static_cast<py::ssize_t>(n)});
    std::copy(src, src + static_cast<std::size_t>(n) * nj, a.mutable_data());
    return a;
}

py::dict integrate(int n, int nj, const py::function& residual, const std::optional<py::function>& jacobian,
                   double t0, const Arr& c0, const std::optional<Arr>& cdot0, const std::vector<double>& t_out,
                   const std::optional<py::array_t<bool, py::array::c_style | py::array::forcecast>>& algebraic,
                   const py::dict& opt) {
    const py::ssize_t N = static_cast<py::ssize_t>(n) * nj;
    require_size(c0, N, "c0");
    const double* cd0 = opt_ptr(cdot0, N, "cdot0");
    std::vector<char> alg;
    if (algebraic) {
        if (algebraic->size() != N) throw std::invalid_argument("algebraic must have shape (nj, n)");
        alg.assign(algebraic->data(), algebraic->data() + N);
    }
    bandsolver::IntegratorOptions o;
    o.adaptive = opt["adaptive"].cast<bool>();
    o.max_order = opt["max_order"].cast<int>();
    o.rtol = opt["rtol"].cast<double>();
    o.atol = opt["atol"].cast<double>();
    o.dt = opt["dt"].cast<double>();
    o.dt0 = opt["dt0"].cast<double>();
    o.dt_min = opt["dt_min"].cast<double>();
    o.dt_max = opt["dt_max"].cast<double>();
    o.jacobian_reuse = opt["jacobian_reuse"].cast<bool>();
    o.reuse_alpha_change = opt["reuse_alpha_change"].cast<double>();
    o.max_newton_iter = opt["max_newton_iter"].cast<int>();
    o.newton_tol = opt["newton_tol"].cast<double>();
    o.max_steps = opt["max_steps"].cast<int>();
    o.fd.rel_step = opt["rel_step"].cast<double>();
    o.fd.typical = opt["typical"].cast<double>();

    bandsolver::DaeResidual R = [&](double t, const double* c, const double* cd, double* F) {
        Arr f = Arr::ensure(residual(t, state_array(n, nj, c), state_array(n, nj, cd)));
        if (!f) throw std::invalid_argument("residual(t, c, cdot) must return a numeric array");
        require_size(f, N, "F");
        std::copy(f.data(), f.data() + N, F);
    };
    bandsolver::DaeJacobian Jf;
    if (jacobian) {
        Jf = [&](double t, const double* c, const double* cd, double alpha, bandsolver::BlockSystem& s) {
            py::object out = (*jacobian)(t, state_array(n, nj, c), state_array(n, nj, cd), alpha);
            auto seq = out.cast<py::sequence>();
            if (seq.size() != 3 && seq.size() != 5)
                throw std::invalid_argument("jacobian(t, c, cdot, alpha) must return (A, B, D) or (A, B, D, X, Y)");
            const py::ssize_t nb = N * n, nn = static_cast<py::ssize_t>(n) * n;
            double* dst[5] = {s.A().data(), s.B().data(), s.D().data(), s.X().data(), s.Y().data()};
            const py::ssize_t sizes[5] = {nb, nb, nb, nn, nn};
            const char* names[5] = {"A", "B", "D", "X", "Y"};
            for (std::size_t i = 0; i < seq.size(); ++i) {
                if (seq[i].is_none()) continue;
                Arr a = Arr::ensure(seq[i]);
                if (!a) throw std::invalid_argument(std::string("jacobian output ") + names[i] + " is not numeric");
                require_size(a, sizes[i], names[i]);
                std::copy(a.data(), a.data() + sizes[i], dst[i]);
            }
        };
    }
    auto r = bandsolver::integrate(n, nj, R, Jf, t0, c0.data(), cd0, t_out, alg, o);
    if (r.callback_exception) std::rethrow_exception(r.callback_exception);
    Arr Y({static_cast<py::ssize_t>(r.y.size()), static_cast<py::ssize_t>(nj), static_cast<py::ssize_t>(n)});
    for (std::size_t k = 0; k < r.y.size(); ++k)
        std::copy(r.y[k].begin(), r.y[k].end(), Y.mutable_data() + k * N);
    py::dict d;
    d["status"] = static_cast<int>(r.status);
    d["message"] = r.message;
    d["t"] = r.t;
    d["y"] = Y;
    d["t_reached"] = r.t_reached;
    py::dict st;
    st["steps"] = r.stats.steps;
    st["rejected_error"] = r.stats.rejected_error;
    st["rejected_newton"] = r.stats.rejected_newton;
    st["newton_iterations"] = r.stats.newton_iterations;
    st["jacobian_evaluations"] = r.stats.jacobian_evaluations;
    st["factorizations"] = r.stats.factorizations;
    st["residual_evaluations"] = r.stats.residual_evaluations;
    d["stats"] = st;
    return d;
}

// One factorization, from either backend.
class PyFactorization {
public:
    PyFactorization(int n, int nj, const Arr& A, const Arr& B, const Arr& D, const std::optional<Arr>& X,
                    const std::optional<Arr>& Y, const std::string& backend)
        : n_(n), nj_(nj), backend_(backend) {
        if (n < 1 || nj < 3) throw std::invalid_argument("require n >= 1 and nj >= 3");
        const py::ssize_t nb = static_cast<py::ssize_t>(n) * n * nj;
        require_size(A, nb, "A"); require_size(B, nb, "B"); require_size(D, nb, "D");
        const double* xp = opt_ptr(X, static_cast<py::ssize_t>(n) * n, "X");
        const double* yp = opt_ptr(Y, static_cast<py::ssize_t>(n) * n, "Y");
        if (backend == "cpp") {
            bandsolver::SystemView v{n, nj, A.data(), B.data(), D.data(), nullptr, xp, yp};
            cpp_ = bandsolver::factor(v);
            status_ = static_cast<int>(cpp_.status());
            fail_node_ = cpp_.fail_node();
        } else if (backend == "fortran") {
            int fnode = 0;
            status_ = bandsolver_f_factor(n, nj, A.data(), B.data(), D.data(), xp, yp, &handle_, &fnode);
            fail_node_ = fnode > 0 ? fnode - 1 : -1;
        } else {
            throw std::invalid_argument("backend must be 'cpp' or 'fortran'");
        }
    }
    ~PyFactorization() { if (handle_) bandsolver_f_factor_free(handle_); }
    PyFactorization(const PyFactorization&) = delete;
    PyFactorization& operator=(const PyFactorization&) = delete;

    py::tuple solve(const Arr& G) const {
        require_size(G, static_cast<py::ssize_t>(n_) * nj_, "G");
        Arr dc({static_cast<py::ssize_t>(nj_), static_cast<py::ssize_t>(n_)});
        int status;
        {
            py::gil_scoped_release release;
            if (backend_ == "cpp")
                status = static_cast<int>(cpp_.solve(G.data(), dc.mutable_data()).status);
            else
                status = bandsolver_f_factor_solve(handle_, G.data(), dc.mutable_data());
        }
        return py::make_tuple(dc, status);
    }
    int status() const { return status_; }
    int fail_node() const { return fail_node_; }
    int n() const { return n_; }
    int nj() const { return nj_; }
    const std::string& backend() const { return backend_; }

private:
    int n_, nj_;
    std::string backend_;
    bandsolver::Factorization cpp_;
    void* handle_ = nullptr;
    int status_ = 2, fail_node_ = -1;
};

}  // namespace

PYBIND11_MODULE(_core, m) {
    m.doc() = "Native backends for bandsolver (C++ core and Fortran library via C ABI).";
    m.def("solve", &solve, py::arg("n"), py::arg("nj"), py::arg("A"), py::arg("B"), py::arg("D"), py::arg("G"),
          py::arg("X"), py::arg("Y"), py::arg("pivot"), py::arg("backend"), py::arg("kernel") = 0,
          py::arg("singular") = 0);
    m.def("newton", &newton, py::arg("n"), py::arg("nj"), py::arg("fill"), py::arg("c0"), py::arg("rtol"),
          py::arg("atol"), py::arg("damping"), py::arg("max_iter"), py::arg("pivot"), py::arg("require_convergence"),
          py::arg("backend"), py::arg("kernel") = 0, py::arg("jacobian_reuse") = false,
          py::arg("reuse_max_iter") = 5, py::arg("reuse_contraction") = 0.5, py::arg("residual") = py::none());
    py::class_<PyFactorization>(m, "Factorization")
        .def(py::init<int, int, const Arr&, const Arr&, const Arr&, const std::optional<Arr>&,
                      const std::optional<Arr>&, const std::string&>(),
             py::arg("n"), py::arg("nj"), py::arg("A"), py::arg("B"), py::arg("D"), py::arg("X"), py::arg("Y"),
             py::arg("backend"))
        .def("solve", &PyFactorization::solve, py::arg("G"))
        .def_property_readonly("status", &PyFactorization::status)
        .def_property_readonly("fail_node", &PyFactorization::fail_node)
        .def_property_readonly("n", &PyFactorization::n)
        .def_property_readonly("nj", &PyFactorization::nj)
        .def_property_readonly("backend", &PyFactorization::backend);
    m.def("integrate", &integrate, py::arg("n"), py::arg("nj"), py::arg("residual"), py::arg("jacobian"),
          py::arg("t0"), py::arg("c0"), py::arg("cdot0"), py::arg("t_out"), py::arg("algebraic"), py::arg("options"));
    m.def("fd_jacobian", &fd_jacobian, py::arg("n"), py::arg("nj"), py::arg("residual"), py::arg("c"),
          py::arg("rel_step"), py::arg("typical"), py::arg("backend"));
    m.def("newton_fd", &newton_fd, py::arg("n"), py::arg("nj"), py::arg("residual"), py::arg("c0"), py::arg("rtol"),
          py::arg("atol"), py::arg("damping"), py::arg("max_iter"), py::arg("pivot"), py::arg("require_convergence"),
          py::arg("rel_step"), py::arg("typical"), py::arg("backend"), py::arg("kernel") = 0,
          py::arg("jacobian_reuse") = false, py::arg("reuse_max_iter") = 5, py::arg("reuse_contraction") = 0.5);
    m.def("check_jacobian", &check_jacobian, py::arg("n"), py::arg("nj"), py::arg("fill"), py::arg("c"),
          py::arg("rel_step"), py::arg("typical"), py::arg("backend"));
}
