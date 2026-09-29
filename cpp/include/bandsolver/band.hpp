// Newman BAND block solver (Electrochemical Systems, Appendix C) — C++17 core.
//
// Solves, for dc[nj][n] with nj >= 3,
//   j = 0:        B_0 dc_0 + D_0 dc_1 + X dc_2                      = G_0
//   0 < j < nj-1: A_j dc_{j-1} + B_j dc_j + D_j dc_{j+1}             = G_j
//   j = nj-1:     Y dc_{nj-3} + A_{nj-1} dc_{nj-2} + B_{nj-1} dc_{nj-1} = G_{nj-1}
// Storage is C row-major: blocks [nj][n][n] (A[j][i][k] couples equation i of node j to
// unknown k of node j-1), vectors [nj][n], X/Y [n][n]. Node indices are 0-based here.
#ifndef BANDSOLVER_BAND_HPP
#define BANDSOLVER_BAND_HPP

#include <stdexcept>
#include <string>
#include <vector>

namespace bandsolver {

enum class Status : int {
    ok = 0,
    singular = 1,
    invalid_argument = 2,
    not_converged = 3,
    non_finite = 4,
    callback_error = 5,
};

enum class Pivot : int {
    partial = 0,  // row partial pivoting (default)
    legacy = 1,   // archival MATINV pivot heuristic, for historical comparison
};

// When a pivot block counts as singular.
enum class Singular : int {
    relative = 0,  // |pivot| <= n*eps*max|block| (default)
    exact = 1,     // only an exactly zero pivot, as in the archival MATINV; with Pivot::legacy
                   // this reproduces the archival kernel on nearly singular blocks too
};

const char* to_string(Status s) noexcept;

// Thrown by the throwing convenience APIs; carries the status and failing node (or -1).
class Error : public std::runtime_error {
public:
    Error(Status status, int node, const std::string& what);
    Status status() const noexcept { return status_; }
    int node() const noexcept { return node_; }

private:
    Status status_;
    int node_;
};

// Non-owning view of a block system. X and Y may be null (treated as zero).
struct SystemView {
    int n = 0;
    int nj = 0;
    const double* A = nullptr;
    const double* B = nullptr;
    const double* D = nullptr;
    const double* G = nullptr;
    const double* X = nullptr;
    const double* Y = nullptr;
};

struct SolveInfo {
    Status status = Status::ok;
    int fail_node = -1;           // 0-based node whose pivot block was singular, else -1
    double min_rel_pivot = 0.0;   // min |pivot| / max|block| over all node factorizations
};

// Solve K dc = G. dc must hold n*nj values. Inputs are not modified. Errors are reported
// through SolveInfo::status; only std::bad_alloc (workspace allocation) can propagate.
SolveInfo solve(const SystemView& sys, double* dc, Pivot pivot = Pivot::partial);
SolveInfo solve(const SystemView& sys, double* dc, Pivot pivot, Singular singular);

// Owning, zero-initialized storage for one system, with element accessors.
class BlockSystem {
public:
    BlockSystem(int n, int nj);

    int n() const noexcept { return n_; }
    int nj() const noexcept { return nj_; }

    double& A(int j, int i, int k) { return A_[blk(j, i, k)]; }
    double& B(int j, int i, int k) { return B_[blk(j, i, k)]; }
    double& D(int j, int i, int k) { return D_[blk(j, i, k)]; }
    double& G(int j, int i) { return G_[static_cast<std::size_t>(j) * n_ + i]; }
    double& X(int i, int k) { return X_[static_cast<std::size_t>(i) * n_ + k]; }
    double& Y(int i, int k) { return Y_[static_cast<std::size_t>(i) * n_ + k]; }

    std::vector<double>& A() noexcept { return A_; }
    std::vector<double>& B() noexcept { return B_; }
    std::vector<double>& D() noexcept { return D_; }
    std::vector<double>& G() noexcept { return G_; }
    std::vector<double>& X() noexcept { return X_; }
    std::vector<double>& Y() noexcept { return Y_; }

    void set_zero();
    SystemView view() const noexcept;

private:
    std::size_t blk(int j, int i, int k) const {
        return (static_cast<std::size_t>(j) * n_ + i) * n_ + k;
    }
    int n_, nj_;
    std::vector<double> A_, B_, D_, G_, X_, Y_;
};

// Throwing convenience: returns dc (size n*nj) or throws bandsolver::Error.
std::vector<double> solve(const BlockSystem& sys, Pivot pivot = Pivot::partial);
std::vector<double> solve(const BlockSystem& sys, Pivot pivot, Singular singular);

}  // namespace bandsolver

#endif
