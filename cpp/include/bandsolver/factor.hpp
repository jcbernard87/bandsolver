// Factor once, solve many: the BAND elimination split into a factorization of the block
// matrix K (depends only on A, B, D, X, Y) and a cheap solve for any right-hand side G.
//
// Per node the factor stores the partially pivoted LU of the reduced pivot block
// B^_j = B_j + A'_j E_{j-1}, the elimination blocks E_j = -B^_j^{-1} D'_j, and the effective
// lower blocks A'_j (including the Y correction at the last node). A solve is then a forward
// sweep e_j = B^_j^{-1} (G_j - A'_j e_{j-1}) plus back substitution: O(nj n^2) instead of the
// O(nj n^3) elimination. Used for Jacobian reuse (modified Newton) and in time integration.
#ifndef BANDSOLVER_FACTOR_HPP
#define BANDSOLVER_FACTOR_HPP

#include <vector>

#include "bandsolver/band.hpp"

namespace bandsolver {

class Factorization {
public:
    Factorization() = default;

    Status status() const noexcept { return status_; }
    int fail_node() const noexcept { return fail_node_; }   // 0-based, or -1
    double min_rel_pivot() const noexcept { return min_rel_pivot_; }
    int n() const noexcept { return n_; }
    int nj() const noexcept { return nj_; }
    bool ok() const noexcept { return status_ == Status::ok; }

    // Solve K dc = G (both n*nj values). Never throws.
    SolveInfo solve(const double* G, double* dc) const;
    // Throwing convenience.
    std::vector<double> solve(const std::vector<double>& G) const;

private:
    friend Factorization factor(const SystemView& sys);
    int n_ = 0, nj_ = 0, fail_node_ = -1;
    Status status_ = Status::invalid_argument;
    double min_rel_pivot_ = 0.0;
    std::vector<double> lu_;    // [nj][n][n] LU of B^_j (unit lower L below the diagonal)
    std::vector<int> perm_;     // [nj][n] row permutation of each LU
    std::vector<double> E_;     // [nj][n][n]
    std::vector<double> Aeff_;  // [nj][n][n]
    std::vector<double> Xp_;    // [n][n]  X' = -B^_0^{-1} X
    std::vector<double> Y_;     // [n][n]
};

// Factor the block matrix of `sys` (sys.G is ignored and may be null). Check status().
Factorization factor(const SystemView& sys);

}  // namespace bandsolver

#endif
