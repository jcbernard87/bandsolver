#include "bandsolver/newton.hpp"

#include <algorithm>
#include <cmath>

namespace bandsolver {

NewtonResult newton(int n, int nj, const FillFunction& fill, double* c, const NewtonOptions& opts) {
    NewtonResult res;
    if (n < 1 || nj < 3 || !fill || !c || opts.max_iter < 1 || !(opts.damping > 0 && opts.damping <= 1) ||
        opts.rtol < 0 || opts.atol < 0 || (opts.rtol == 0 && opts.atol == 0)) {
        res.status = Status::invalid_argument;
        return res;
    }
    const std::size_t nv = static_cast<std::size_t>(n) * nj;
    BlockSystem sys(n, nj);
    std::vector<double> dc(nv);

    for (int it = 1; it <= opts.max_iter; ++it) {
        res.iterations = it;
        sys.set_zero();
        try {
            fill(c, sys);
        } catch (...) {
            res.status = Status::callback_error;
            res.callback_exception = std::current_exception();
            return res;
        }
        double rnorm = 0;
        for (double g : sys.G()) rnorm = std::max(rnorm, std::abs(g));
        res.residual_norm.push_back(rnorm);

        const SolveInfo info = solve(sys.view(), dc.data(), opts.pivot);
        if (info.status != Status::ok) {
            res.status = info.status;
            res.fail_node = info.fail_node;
            return res;
        }
        double snorm = 0, unorm = 0;
        for (std::size_t i = 0; i < nv; ++i) {
            c[i] += opts.damping * dc[i];
            snorm = std::max(snorm, std::abs(dc[i]));
            unorm = std::max(unorm, std::abs(dc[i]) / (opts.atol + opts.rtol * std::abs(c[i])));
        }
        res.step_norm.push_back(snorm);
        res.update_norm.push_back(unorm);
        if (unorm <= 1) {
            res.converged = true;
            break;
        }
    }
    if (!res.converged && opts.require_convergence) res.status = Status::not_converged;
    return res;
}

}  // namespace bandsolver
