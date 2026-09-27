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
                int kernel) {
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
            info = bandsolver::solve(v, dc.mutable_data(), static_cast<bandsolver::Pivot>(pivot));
        }
        status = static_cast<int>(info.status);
        fail_node = info.fail_node;
        min_rel_pivot = info.min_rel_pivot;
    } else if (backend == "fortran") {
        {
            py::gil_scoped_release release;
            status = bandsolver_f_solve_kernel(n, nj, A.data(), B.data(), D.data(), G.data(), xp, yp, pivot, kernel,
                                               dc.mutable_data(), &fail_node, &min_rel_pivot);
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

py::dict newton(int n, int nj, const py::function& fill, const Arr& c0, double rtol, double atol, double damping,
                int max_iter, int pivot, bool require_convergence, const std::string& backend, int kernel) {
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
        auto res = bandsolver::newton(
            n, nj,
            [&](const double* cc, bandsolver::BlockSystem& s) {
                call_fill(fill, n, nj, cc, s.A().data(), s.B().data(), s.D().data(), s.G().data(), s.X().data(),
                          s.Y().data());
            },
            c.mutable_data(), o);
        if (res.callback_exception) std::rethrow_exception(res.callback_exception);
        r["status"] = static_cast<int>(res.status);
        r["iterations"] = res.iterations;
        r["converged"] = res.converged;
        r["fail_node"] = res.fail_node;
        r["update_norm"] = res.update_norm;
        r["step_norm"] = res.step_norm;
        r["residual_norm"] = res.residual_norm;
    } else if (backend == "fortran") {
        bandsolver_newton_options o{rtol, atol, damping, max_iter, pivot, require_convergence ? 1 : 0, kernel};
        bandsolver_newton_result res;
        std::vector<double> hu(max_iter), hs(max_iter), hr(max_iter);
        FortranCtx ctx{&fill, nullptr};
        bandsolver_f_newton(n, nj, fortran_trampoline, &ctx, c.mutable_data(), &o, &res, hu.data(), hs.data(),
                            hr.data());
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
                   const std::string& backend, int kernel) {
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
    } else if (backend == "fortran") {
        bandsolver_newton_options no{rtol, atol, damping, max_iter, pivot, require_convergence ? 1 : 0, kernel};
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

}  // namespace

PYBIND11_MODULE(_core, m) {
    m.doc() = "Native backends for bandsolver (C++ core and Fortran library via C ABI).";
    m.def("solve", &solve, py::arg("n"), py::arg("nj"), py::arg("A"), py::arg("B"), py::arg("D"), py::arg("G"),
          py::arg("X"), py::arg("Y"), py::arg("pivot"), py::arg("backend"), py::arg("kernel") = 0);
    m.def("newton", &newton, py::arg("n"), py::arg("nj"), py::arg("fill"), py::arg("c0"), py::arg("rtol"),
          py::arg("atol"), py::arg("damping"), py::arg("max_iter"), py::arg("pivot"), py::arg("require_convergence"),
          py::arg("backend"), py::arg("kernel") = 0);
    m.def("fd_jacobian", &fd_jacobian, py::arg("n"), py::arg("nj"), py::arg("residual"), py::arg("c"),
          py::arg("rel_step"), py::arg("typical"), py::arg("backend"));
    m.def("newton_fd", &newton_fd, py::arg("n"), py::arg("nj"), py::arg("residual"), py::arg("c0"), py::arg("rtol"),
          py::arg("atol"), py::arg("damping"), py::arg("max_iter"), py::arg("pivot"), py::arg("require_convergence"),
          py::arg("rel_step"), py::arg("typical"), py::arg("backend"), py::arg("kernel") = 0);
    m.def("check_jacobian", &check_jacobian, py::arg("n"), py::arg("nj"), py::arg("fill"), py::arg("c"),
          py::arg("rel_step"), py::arg("typical"), py::arg("backend"));
}
