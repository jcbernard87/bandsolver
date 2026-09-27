// pybind11 extension `bandsolver._core`: exposes the C++ core and the Fortran library
// (through its C ABI) with identical call signatures. Shape validation and exception
// mapping live in the Python package; this layer checks only what memory safety needs.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

#include <bandsolver/band.hpp>
#include <bandsolver/newton.hpp>

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
                const std::optional<Arr>& X, const std::optional<Arr>& Y, int pivot, const std::string& backend) {
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
            status = bandsolver_f_solve(n, nj, A.data(), B.data(), D.data(), G.data(), xp, yp, pivot,
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
                int max_iter, int pivot, bool require_convergence, const std::string& backend) {
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
        bandsolver_newton_options o{rtol, atol, damping, max_iter, pivot, require_convergence ? 1 : 0};
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

}  // namespace

PYBIND11_MODULE(_core, m) {
    m.doc() = "Native backends for bandsolver (C++ core and Fortran library via C ABI).";
    m.def("solve", &solve, py::arg("n"), py::arg("nj"), py::arg("A"), py::arg("B"), py::arg("D"), py::arg("G"),
          py::arg("X"), py::arg("Y"), py::arg("pivot"), py::arg("backend"));
    m.def("newton", &newton, py::arg("n"), py::arg("nj"), py::arg("fill"), py::arg("c0"), py::arg("rtol"),
          py::arg("atol"), py::arg("damping"), py::arg("max_iter"), py::arg("pivot"), py::arg("require_convergence"),
          py::arg("backend"));
}
