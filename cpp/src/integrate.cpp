#include "bandsolver/integrate.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <utility>

#include "bandsolver/factor.hpp"

namespace bandsolver {

namespace {

struct Point {
    double t;
    std::vector<double> y;
};

// Lagrange extrapolation through the given points, evaluated at time T.
void extrapolate(const std::deque<Point>& pts, std::size_t count, double T, std::vector<double>& out) {
    const std::size_t N = out.size(), m = pts.size();
    std::fill(out.begin(), out.end(), 0.0);
    for (std::size_t i = m - count; i < m; ++i) {
        double L = 1.0;
        for (std::size_t k = m - count; k < m; ++k)
            if (k != i) L *= (T - pts[k].t) / (pts[i].t - pts[k].t);
        for (std::size_t q = 0; q < N; ++q) out[q] += L * pts[i].y[q];
    }
}

}  // namespace

IntegrationResult integrate(int n, int nj, const DaeResidual& residual, const DaeJacobian& jacobian, double t0,
                            const double* c0, const double* cdot0, const std::vector<double>& t_out,
                            const std::vector<char>& algebraic, const IntegratorOptions& o) {
    IntegrationResult res;
    res.t_reached = t0;
    const std::size_t N = static_cast<std::size_t>(n) * nj;
    bool bad = n < 1 || nj < 3 || !residual || !c0 || t_out.empty() || (o.max_order != 1 && o.max_order != 2) ||
               !(o.rtol >= 0) || !(o.atol >= 0) || (o.rtol == 0 && o.atol == 0) || o.max_newton_iter < 1 ||
               !(o.newton_tol > 0) || (!o.adaptive && !(o.dt > 0)) || (!algebraic.empty() && algebraic.size() != N);
    for (std::size_t k = 0; !bad && k < t_out.size(); ++k)
        bad = !(t_out[k] > (k ? t_out[k - 1] : t0));
    if (bad) {
        res.status = Status::invalid_argument;
        res.message = "invalid arguments or options";
        return res;
    }
    const double t_end = t_out.back();
    const double span = t_end - t0;
    const double eps_t = 64 * std::numeric_limits<double>::epsilon() * std::max(std::abs(t0), std::abs(t_end));

    std::size_t n_diff = 0;
    for (std::size_t i = 0; i < N; ++i) n_diff += algebraic.empty() || !algebraic[i];

    std::deque<Point> hist;
    hist.push_back({t0, std::vector<double>(c0, c0 + N)});
    std::vector<double> ydot0(N, 0.0);
    if (cdot0) std::copy(cdot0, cdot0 + N, ydot0.begin());

    std::vector<double> c(N), pred(N), cdot(N), beta(N), F(N), dc(N), w(N);
    BlockSystem J(n, nj);
    Factorization fac;
    bool have_fac = false;
    double alpha_fac = 0;
    std::size_t next_out = 0;
    double h = o.adaptive ? (o.dt0 > 0 ? o.dt0 : span * 1e-4) : o.dt;
    if (o.dt_max > 0) h = std::min(h, o.dt_max);
    double h_prev = 0;

    auto wrms = [&](const std::vector<double>& v, bool diff_only) {
        double s = 0;
        std::size_t count = 0;
        for (std::size_t i = 0; i < N; ++i) {
            if (diff_only && !algebraic.empty() && algebraic[i]) continue;
            const double x = v[i] * w[i];
            s += x * x;
            ++count;
        }
        return count ? std::sqrt(s / count) : 0.0;
    };

    try {
        while (true) {
            const Point& cur = hist.back();
            if (cur.t >= t_end - eps_t) break;
            if (res.stats.steps >= o.max_steps) {
                res.status = Status::not_converged;
                res.message = "max_steps reached";
                break;
            }
            // Land exactly on the next output time.
            const double target = t_out[next_out];
            bool hits_output = false;
            double hstep = h;
            if (cur.t + hstep >= target - eps_t) {
                hstep = target - cur.t;
                hits_output = true;
            }
            const double tn = cur.t + hstep;

            // Order and BDF coefficients: cdot = alpha c + beta.
            int q = std::min<int>(o.max_order, static_cast<int>(hist.size()));
            const double omega = h_prev > 0 ? hstep / h_prev : 1.0;
            if (q == 2 && omega > 2.4) q = 1;   // variable-step BDF2 zero-stability limit (1 + sqrt 2)
            double alpha;
            if (q == 1) {
                alpha = 1.0 / hstep;
                for (std::size_t i = 0; i < N; ++i) beta[i] = -cur.y[i] / hstep;
            } else {
                const auto& prev = hist[hist.size() - 2].y;
                alpha = (1 + 2 * omega) / ((1 + omega) * hstep);
                for (std::size_t i = 0; i < N; ++i)
                    beta[i] = (-(1 + omega) * cur.y[i] + omega * omega / (1 + omega) * prev[i]) / hstep;
            }
            // Predictor (q+1 points when available) and its error constant.
            const std::size_t npts = std::min<std::size_t>(q + 1, hist.size());
            double K;
            if (npts == static_cast<std::size_t>(q) + 1) {
                extrapolate(hist, npts, tn, pred);
                K = q == 1 ? 1.0 / 3.0 : 2.0 / 11.0;
            } else {
                for (std::size_t i = 0; i < N; ++i) pred[i] = cur.y[i] + hstep * ydot0[i];
                K = 0.5;
            }
            for (std::size_t i = 0; i < N; ++i) w[i] = 1.0 / (o.atol + o.rtol * std::abs(cur.y[i]));

            // Newton on the step.
            c = pred;
            bool converged = false, fresh = false;
            if (!o.jacobian_reuse || !have_fac || std::abs(alpha / alpha_fac - 1) > o.reuse_alpha_change)
                have_fac = false;
            for (int attempt = 0; attempt < 2 && !converged; ++attempt) {
                if (attempt == 1) {
                    if (fresh) break;           // a fresh Jacobian already failed: shrink the step
                    have_fac = false;           // retry once with a fresh Jacobian
                    c = pred;
                }
                double del_prev = 0;
                for (int k = 0; k < o.max_newton_iter; ++k) {
                    for (std::size_t i = 0; i < N; ++i) cdot[i] = alpha * c[i] + beta[i];
                    if (!have_fac || !o.jacobian_reuse) {
                        J.set_zero();
                        if (jacobian) {
                            jacobian(tn, c.data(), cdot.data(), alpha, J);
                        } else {
                            auto r = [&](const double* x, double* Fx) {
                                std::vector<double> xd(N);
                                for (std::size_t i = 0; i < N; ++i) xd[i] = alpha * x[i] + beta[i];
                                residual(tn, x, xd.data(), Fx);
                            };
                            res.stats.residual_evaluations += fd_jacobian(n, nj, r, c.data(), J, o.fd);
                        }
                        ++res.stats.jacobian_evaluations;
                        fac = factor(J.view());
                        ++res.stats.factorizations;
                        if (!fac.ok()) {
                            res.status = fac.status();
                            res.message = "singular iteration matrix";
                            res.t_reached = cur.t;
                            return res;
                        }
                        have_fac = true;
                        fresh = true;
                        alpha_fac = alpha;
                    }
                    residual(tn, c.data(), cdot.data(), F.data());
                    ++res.stats.residual_evaluations;
                    ++res.stats.newton_iterations;
                    for (auto& v : F) v = -v;
                    const SolveInfo si = fac.solve(F.data(), dc.data());
                    if (si.status != Status::ok) break;
                    // A factorization made at alpha_fac approximates alpha: rescale (as in IDA).
                    const double scale = 2.0 / (1.0 + alpha / alpha_fac);
                    for (std::size_t i = 0; i < N; ++i) c[i] += scale * dc[i];
                    for (auto& v : dc) v *= scale;
                    const double del = wrms(dc, false);
                    if (k == 0) {
                        if (del <= 0.1 * o.newton_tol) { converged = true; break; }
                    } else {
                        const double rate = del / del_prev;
                        if (rate > 0.9) break;
                        if (rate / (1 - rate) * del <= o.newton_tol) { converged = true; break; }
                    }
                    del_prev = del;
                    if (!std::isfinite(del)) break;
                }
                if (!o.jacobian_reuse) break;
            }
            if (!converged) {
                ++res.stats.rejected_newton;
                have_fac = false;
                h = hstep * 0.25;
                if ((o.dt_min > 0 && h < o.dt_min) || h < eps_t || !o.adaptive) {
                    res.status = Status::not_converged;
                    res.message = "Newton failed to converge";
                    res.t_reached = cur.t;
                    return res;
                }
                continue;
            }

            // Local error estimate and step control.
            for (std::size_t i = 0; i < N; ++i) dc[i] = K * (c[i] - pred[i]);
            const double err = n_diff ? wrms(dc, true) : 0.0;
            const double expo = 1.0 / (q + 1);
            if (o.adaptive && err > 1.0) {
                ++res.stats.rejected_error;
                h = hstep * std::max(0.2, 0.9 * std::pow(err, -expo));
                if ((o.dt_min > 0 && h < o.dt_min) || h < eps_t) {
                    res.status = Status::not_converged;
                    res.message = "step size underflow";
                    res.t_reached = cur.t;
                    return res;
                }
                continue;
            }
            // Accept.
            for (std::size_t i = 0; i < N; ++i) ydot0[i] = alpha * c[i] + beta[i];
            hist.push_back({hits_output ? target : tn, c});
            if (hist.size() > 3) hist.pop_front();
            h_prev = hstep;
            ++res.stats.steps;
            res.t_reached = hist.back().t;
            if (hits_output) {
                res.t.push_back(target);
                res.y.push_back(c);
                ++next_out;
            }
            if (o.adaptive) {
                const double grow = err > 0 ? std::min(2.0, std::max(0.2, 0.9 * std::pow(err, -expo))) : 2.0;
                h = (hits_output ? std::max(h, hstep) : hstep) * grow;
                if (o.dt_max > 0) h = std::min(h, o.dt_max);
            } else {
                h = o.dt;
            }
        }
    } catch (...) {
        res.status = Status::callback_error;
        res.callback_exception = std::current_exception();
        res.message = "exception in a user callback";
        return res;
    }
    return res;
}

}  // namespace bandsolver
