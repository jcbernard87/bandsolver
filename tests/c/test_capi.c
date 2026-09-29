/* T3: C-ABI round trip against the Fortran library (row-major layout, callbacks). */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "bandsolver_f.h"

static int failures = 0;
static void check(int cond, const char *msg) {
    printf("%s %s\n", cond ? "PASS" : "FAIL", msg);
    if (!cond) failures++;
}
static double urand(void) { return 2.0 * rand() / (double)RAND_MAX - 1.0; }

/* Backward error of the row-major system, computed directly from the blocks. */
static double backward_error(int n, int nj, const double *A, const double *B, const double *D,
                             const double *G, const double *X, const double *Y, const double *dc) {
    double err = 0, knorm = 0, xnorm = 0, gnorm = 0;
    for (int j = 0; j < nj; j++)
        for (int i = 0; i < n; i++) {
            double r = -G[j*n + i], rs = 0;
            for (int k = 0; k < n; k++) {
                size_t b = ((size_t)j*n + i)*n + k;
                r += B[b]*dc[j*n + k]; rs += fabs(B[b]);
                if (j > 0)      { r += A[b]*dc[(j-1)*n + k]; rs += fabs(A[b]); }
                if (j < nj - 1) { r += D[b]*dc[(j+1)*n + k]; rs += fabs(D[b]); }
                if (j == 0 && X)      { r += X[i*n + k]*dc[2*n + k]; rs += fabs(X[i*n + k]); }
                if (j == nj - 1 && Y) { r += Y[i*n + k]*dc[(nj-3)*n + k]; rs += fabs(Y[i*n + k]); }
            }
            err = fmax(err, fabs(r)); knorm = fmax(knorm, rs);
            xnorm = fmax(xnorm, fabs(dc[j*n + i])); gnorm = fmax(gnorm, fabs(G[j*n + i]));
        }
    return err / (knorm*xnorm + gnorm);
}

typedef struct { int fail; } diff_ctx;

/* Linear -c'' = pi^2 sin(pi x), c(0)=c(1)=0 on nj nodes. */
static int fill_linear(int n, int nj, const double *c, double *A, double *B, double *D,
                       double *G, double *X, double *Y, void *ctx) {
    (void)n; (void)X; (void)Y;
    if (((diff_ctx *)ctx)->fail) return 7;
    double h = 1.0/(nj - 1), pi = acos(-1.0);
    B[0] = 1; G[0] = -c[0];
    B[nj-1] = 1; G[nj-1] = -c[nj-1];
    for (int j = 1; j < nj - 1; j++) {
        A[j] = -1/(h*h); B[j] = 2/(h*h); D[j] = -1/(h*h);
        double F = -(c[j+1] - 2*c[j] + c[j-1])/(h*h) - pi*pi*sin(pi*j*h);
        G[j] = -F;
    }
    return 0;
}

int main(void) {
    srand(99);
    const int n = 3, nj = 40;
    size_t nb = (size_t)n*n*nj;
    double *A = malloc(nb*sizeof *A), *B = malloc(nb*sizeof *B), *D = malloc(nb*sizeof *D);
    double *G = malloc(n*nj*sizeof *G), *dc = malloc(n*nj*sizeof *dc), X[9], Y[9];
    for (size_t i = 0; i < nb; i++) { A[i] = urand(); B[i] = urand(); D[i] = urand(); }
    for (int j = 0; j < nj; j++) for (int i = 0; i < n; i++) B[((size_t)j*n + i)*n + i] += 3*n + 3;
    for (int i = 0; i < n*nj; i++) G[i] = urand();
    for (int i = 0; i < 9; i++) { X[i] = urand(); Y[i] = urand(); }

    int fail_node = -1; double mrp = 0;
    for (int piv = 0; piv < 2; piv++) {
        int st = bandsolver_f_solve(n, nj, A, B, D, G, X, Y, piv, dc, &fail_node, &mrp);
        double be = backward_error(n, nj, A, B, D, G, X, Y, dc);
        printf("pivot=%d backward error %.2e, min rel pivot %.2e\n", piv, be, mrp);
        check(st == BANDSOLVER_OK && be < 1e-13, "C-ABI solve with X/Y (row-major, nonsymmetric blocks)");
    }
    int st = bandsolver_f_solve(n, nj, A, B, D, G, NULL, NULL, 0, dc, &fail_node, &mrp);
    check(st == BANDSOLVER_OK && backward_error(n, nj, A, B, D, G, NULL, NULL, dc) < 1e-13, "C-ABI solve with NULL X/Y");
    for (int k = 0; k < n*n; k++) B[5*n*n + k] = 0; for (int k = 0; k < n*n; k++) A[5*n*n + k] = 0;
    st = bandsolver_f_solve(n, nj, A, B, D, G, NULL, NULL, 0, dc, &fail_node, &mrp);
    check(st == BANDSOLVER_SINGULAR && fail_node == 6, "C-ABI singular block reports 1-based node 6");
    st = bandsolver_f_solve(n, 2, A, B, D, G, NULL, NULL, 0, dc, &fail_node, &mrp);
    check(st == BANDSOLVER_INVALID_ARGUMENT, "C-ABI rejects nj < 3");
    st = bandsolver_f_solve_ex(n, nj, A, B, D, G, NULL, NULL, BANDSOLVER_PIVOT_LEGACY, BANDSOLVER_KERNEL_FAST,
                               BANDSOLVER_SINGULAR_EXACT, dc, &fail_node, &mrp);
    check(st == BANDSOLVER_SINGULAR && fail_node == 6, "C-ABI exact singular rule still reports a zero block");
    st = bandsolver_f_solve_ex(n, nj, A, B, D, G, NULL, NULL, 0, 0, 3, dc, &fail_node, &mrp);
    check(st == BANDSOLVER_INVALID_ARGUMENT, "C-ABI rejects an unknown singular rule");

    const int m = 51;
    double c[51] = {0}, hist_u[50], hist_s[50], hist_r[50];
    diff_ctx ctx = {0};
    bandsolver_newton_options opts; bandsolver_newton_result res;
    bandsolver_f_default_options(&opts);
    check(opts.max_iter == 50 && opts.damping == 1.0 && opts.require_convergence == 1, "default options");
    st = bandsolver_f_newton(1, m, fill_linear, &ctx, c, &opts, &res, hist_u, hist_s, hist_r);
    double err = 0, pi = acos(-1.0);
    for (int j = 0; j < m; j++) err = fmax(err, fabs(c[j] - sin(pi*j/(m - 1.0))));
    printf("newton: status %d iterations %d, residual history %.2e %.2e, error %.2e\n",
           st, res.iterations, hist_r[0], hist_r[1], err);
    check(st == BANDSOLVER_OK && res.converged && res.iterations == 2, "C-ABI Newton converges");
    check(err < 1e-3, "C-ABI Newton solution accurate");
    ctx.fail = 1;
    st = bandsolver_f_newton(1, m, fill_linear, &ctx, c, &opts, &res, NULL, NULL, NULL);
    check(st == BANDSOLVER_CALLBACK_ERROR, "C callback error propagates");

    free(A); free(B); free(D); free(G); free(dc);
    return failures ? 1 : 0;
}
