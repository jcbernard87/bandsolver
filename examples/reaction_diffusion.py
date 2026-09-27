"""Coupled differential-algebraic system with three unknowns per node (n = 3).

    -u'' + w         = f1      (diffusion with a source that depends on w)
    -v'' + (v - u^2) = f2      (diffusion with a nonlinear exchange term)
     w - u v         = 0       (algebraic constraint, like an electroneutrality row)

on [0,1] with Dirichlet data for u and v from the manufactured solution
u* = sin(pi x), v* = 1 + x(1 - x), w* = u* v*. The Jacobian blocks are full 3x3 at the
nodes and couple unknowns across species, which exercises the block (not scalar) path.
"""
import numpy as np

import bandsolver as bs


def exact(x):
    u = np.sin(np.pi * x)
    v = 1 + x * (1 - x)
    return u, v, u * v


def make_fill():
    def fill(c):
        nj = c.shape[0]
        h = 1.0 / (nj - 1)
        x = np.linspace(0, 1, nj)
        us, vs, _ = exact(x)
        f1 = np.pi**2 * us + us * vs
        f2 = 2 + (vs - us**2)
        u, v, w = c[:, 0], c[:, 1], c[:, 2]
        A = np.zeros((nj, 3, 3)); B = np.zeros((nj, 3, 3)); D = np.zeros((nj, 3, 3)); G = np.zeros((nj, 3))
        i = np.arange(1, nj - 1)
        G[i, 0] = -(-(u[i + 1] - 2 * u[i] + u[i - 1]) / h**2 + w[i] - f1[i])
        G[i, 1] = -(-(v[i + 1] - 2 * v[i] + v[i - 1]) / h**2 + v[i] - u[i]**2 - f2[i])
        A[i, 0, 0] = D[i, 0, 0] = A[i, 1, 1] = D[i, 1, 1] = -1 / h**2
        B[i, 0, 0] = 2 / h**2
        B[i, 0, 2] = 1
        B[i, 1, 1] = 2 / h**2 + 1
        B[i, 1, 0] = -2 * u[i]
        # algebraic row at every node, including the boundaries
        G[:, 2] = -(w - u * v)
        B[:, 2, 2] = 1
        B[:, 2, 0] = -v
        B[:, 2, 1] = -u
        for j in (0, nj - 1):
            B[j, 0, 0] = B[j, 1, 1] = 1
            G[j, 0], G[j, 1] = -(u[j] - us[j]), -(v[j] - vs[j])
        return A, B, D, G

    return fill


def run(backend="cpp", grids=(21, 41, 81, 161, 321)):
    rows = []
    for nj in grids:
        x = np.linspace(0, 1, nj)
        c0 = np.zeros((nj, 3)); c0[:, 1] = 1
        r = bs.newton(make_fill(), c0, backend=backend)
        err = np.abs(r.c - np.stack(exact(x), axis=1)).max()
        rows.append(dict(nj=nj, iterations=r.iterations, error=err, steps=r.step_norm))
    for a, b in zip(rows[:-1], rows[1:]):
        b["order"] = np.log2(a["error"] / b["error"])
    return rows


if __name__ == "__main__":
    for backend in bs.BACKENDS:
        print(f"backend={backend}")
        for row in run(backend):
            print("  nj={nj:4d} iterations={iterations} max error={error:.3e} order={o}".format(
                o=f"{row['order']:.3f}" if "order" in row else "  -  ", **row))
        print("  step norms (nj=81):", " ".join(f"{s:.1e}" for s in run(backend, (81,))[0]["steps"]))
