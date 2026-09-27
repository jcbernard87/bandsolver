/* C interface to the Fortran bandsolver library.
 *
 * Layout is C row-major: blocks A, B, D are [nj][n][n] (A[j][i][k] couples equation i of
 * node j to unknown k of node j-1), G, c, dc are [nj][n], and X, Y are [n][n]. X and Y may
 * be NULL (treated as zero). Node indices reported through fail_node are 1-based.
 */
#ifndef BANDSOLVER_F_H
#define BANDSOLVER_F_H

#ifdef __cplusplus
extern "C" {
#endif

enum bandsolver_status {
    BANDSOLVER_OK = 0,
    BANDSOLVER_SINGULAR = 1,
    BANDSOLVER_INVALID_ARGUMENT = 2,
    BANDSOLVER_NOT_CONVERGED = 3,
    BANDSOLVER_NON_FINITE = 4,
    BANDSOLVER_CALLBACK_ERROR = 5
};

enum bandsolver_pivot { BANDSOLVER_PIVOT_PARTIAL = 0, BANDSOLVER_PIVOT_LEGACY = 1 };

/* Fortran loop organisation: FAST (default, column-major loops) and REFERENCE (archival
 * row-wise loops) give bit-identical results; REFERENCE exists for comparison/benchmarks. */
enum bandsolver_kernel { BANDSOLVER_KERNEL_FAST = 0, BANDSOLVER_KERNEL_REFERENCE = 1 };

/* Fill callback: evaluate Jacobian blocks and G = -F(c) at state c. Output arrays arrive
 * zeroed. Return nonzero to abort with BANDSOLVER_CALLBACK_ERROR. */
typedef int (*bandsolver_fill_fn)(int n, int nj, const double *c, double *A, double *B,
                                  double *D, double *G, double *X, double *Y, void *ctx);

typedef struct {
    double rtol;               /* default 1e-10 */
    double atol;               /* default 1e-12 */
    double damping;            /* in (0, 1], default 1 */
    int max_iter;              /* default 50 */
    int pivot;                 /* enum bandsolver_pivot */
    int require_convergence;   /* default 1; 0 with max_iter=1 = archival one-step use */
    int kernel;                /* enum bandsolver_kernel, default FAST (added in 0.1.2) */
} bandsolver_newton_options;

typedef struct {
    int status;
    int iterations;
    int converged;
    int fail_node;
    double update_norm;        /* last scaled update norm */
    double step_norm;          /* last max |dc| */
    double residual_norm;      /* last max |G| */
} bandsolver_newton_result;

void bandsolver_f_default_options(bandsolver_newton_options *opts);

int bandsolver_f_solve(int n, int nj, const double *A, const double *B, const double *D,
                       const double *G, const double *X, const double *Y, int pivot,
                       double *dc, int *fail_node, double *min_rel_pivot);

/* As bandsolver_f_solve, with an explicit kernel (enum bandsolver_kernel). */
int bandsolver_f_solve_kernel(int n, int nj, const double *A, const double *B, const double *D,
                              const double *G, const double *X, const double *Y, int pivot, int kernel,
                              double *dc, int *fail_node, double *min_rel_pivot);

/* c is updated in place. Histories (length >= max_iter) may be NULL. */
int bandsolver_f_newton(int n, int nj, bandsolver_fill_fn fill, void *ctx, double *c,
                        const bandsolver_newton_options *opts, bandsolver_newton_result *res,
                        double *update_history, double *step_history, double *residual_history);

/* ---- Factor once, solve many ---------------------------------------------------------- */

/* Factor the block matrix (A, B, D, X, Y; row-major as for bandsolver_f_solve). On return
 * *handle owns the factorization even on failure; free it with bandsolver_f_factor_free.
 * fail_node is the 1-based node of a singular block (0 otherwise). */
int bandsolver_f_factor(int n, int nj, const double *A, const double *B, const double *D,
                        const double *X, const double *Y, void **handle, int *fail_node);
/* Solve K dc = G ([nj][n]) with a factorization. */
int bandsolver_f_factor_solve(void *handle, const double *G, double *dc);
void bandsolver_f_factor_free(void *handle);

/* ---- Finite-difference Jacobians (see docs/math.md) ---------------------------------- */

/* Residual callback: evaluate F(c) ([nj][n]). Return nonzero to abort. */
typedef int (*bandsolver_residual_fn)(int n, int nj, const double *c, double *F, void *ctx);

typedef struct {
    double rel_step;           /* default sqrt(eps) = 1.4901161193847656e-8 */
    double typical;            /* default 1 */
} bandsolver_fd_options;

typedef struct {               /* error = |user-fd| / max(|user|, |fd|, 1e-3*rowscale) */
    double error;
    int node, row, col;        /* 1-based; 0 if the block has no entries */
    double user, fd;
} bandsolver_jacobian_mismatch;

typedef struct {
    bandsolver_jacobian_mismatch A, B, D, X, Y;
} bandsolver_jacobian_check;

void bandsolver_f_default_fd_options(bandsolver_fd_options *opts);

/* Blocks by finite differences (3n+1 residual evaluations) and G = -F(c). fd_opts and
 * evaluations may be NULL. Output layout as for bandsolver_f_solve. */
int bandsolver_f_fd_jacobian(int n, int nj, bandsolver_residual_fn residual, void *ctx, const double *c,
                             const bandsolver_fd_options *fd_opts, double *A, double *B, double *D,
                             double *G, double *X, double *Y, long *evaluations);

/* Newton with finite-difference Jacobians. fd_opts, histories and evaluations may be NULL. */
int bandsolver_f_newton_fd(int n, int nj, bandsolver_residual_fn residual, void *ctx, double *c,
                           const bandsolver_newton_options *opts, const bandsolver_fd_options *fd_opts,
                           bandsolver_newton_result *res, double *update_history, double *step_history,
                           double *residual_history, long *evaluations);

/* Compare a fill callback's blocks with finite differences of its own G. fd_opts may be NULL. */
int bandsolver_f_check_jacobian(int n, int nj, bandsolver_fill_fn fill, void *ctx, const double *c,
                                const bandsolver_fd_options *fd_opts, bandsolver_jacobian_check *check);

#ifdef __cplusplus
}
#endif
#endif
