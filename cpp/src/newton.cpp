#include "bandsolver/newton.hpp"

#include "bandsolver/factor.hpp"

#include <algorithm>
#include <cmath>

namespace bandsolver {

NewtonResult newton(int n, int nj, const FillFunction& fill, double* c, const NewtonOptions& opts,
                    const ResidualFunction& residual) {
    NewtonResult res;
    if (n < 1 || nj < 3 || !fill || !c || opts.max_iter < 1 || !(opts.damping > 0 && opts.damping <= 1) ||
        opts.rtol < 0 || opts.atol < 0 || (opts.rtol == 0 && opts.atol == 0) ||
        (opts.jacobian_reuse && (opts.reuse_max_iter < 1 || !(opts.reuse_contraction > 0)))) {
        res.status = Status::invalid_argument;
        return res;
    }
    const std::size_t nv = static_cast<std::size_t>(n) * nj;
    BlockSystem sys(n, nj);
    std::vector<double> dc(nv), F(nv);
    Factorization fac;
    bool refresh = true;
    int uses = 0;
    double prev_step = 0;

    for (int it = 1; it <= opts.max_iter; ++it) {
        res.iterations = it;
        const bool need_jacobian = !opts.jacobian_reuse || refresh;
        const double* G = nullptr;
        try {
            if (need_jacobian || !residual) {
                sys.set_zero();
                fill(c, sys);
                ++res.jacobian_evaluations;
                G = sys.G().data();
            } else {
                residual(c, F.data());
                ++res.residual_evaluations;
                for (auto& v : F) v = -v;
                G = F.data();
            }
        } catch (...) {
            res.status = Status::callback_error;
            res.callback_exception = std::current_exception();
            return res;
        }
        double rnorm = 0;
        for (std::size_t i = 0; i < nv; ++i) rnorm = std::max(rnorm, std::abs(G[i]));
        res.residual_norm.push_back(rnorm);

        SolveInfo info;
        if (!opts.jacobian_reuse) {
            info = solve(sys.view(), dc.data(), opts.pivot);
            ++res.factorizations;
        } else {
            if (need_jacobian) {
                fac = factor(sys.view());
                ++res.factorizations;
                uses = 0;
                refresh = false;
                if (!fac.ok()) {
                    res.status = fac.status();
                    res.fail_node = fac.fail_node();
                    return res;
                }
            }
            info = fac.solve(G, dc.data());
            ++uses;
        }
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
        if (opts.jacobian_reuse)
            refresh = uses >= opts.reuse_max_iter || (uses > 1 && snorm > opts.reuse_contraction * prev_step);
        prev_step = snorm;
    }
    if (!res.converged && opts.require_convergence) res.status = Status::not_converged;
    return res;
}

}  // namespace bandsolver
