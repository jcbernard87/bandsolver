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

/* c is updated in place. Histories (length >= max_iter) may be NULL. */
int bandsolver_f_newton(int n, int nj, bandsolver_fill_fn fill, void *ctx, double *c,
                        const bandsolver_newton_options *opts, bandsolver_newton_result *res,
                        double *update_history, double *step_history, double *residual_history);

#ifdef __cplusplus
}
#endif
#endif
