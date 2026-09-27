"""Why the X / Y endpoint blocks exist: second-order Neumann boundaries.

    -c'' + c = (pi^2 + 1) cos(pi x) on [0,1],  c'(0) = c'(1) = 0,   exact c* = cos(pi x).

A second-order one-sided flux boundary, (-3 c_0 + 4 c_1 - c_2) / 2h = 0, couples node 0
to node 2, which is the X block; the mirror condition at the right end uses Y. Dropping
to the first-order closure (c_1 - c_0)/h = 0 needs no X/Y but degrades the whole solution
to first order. The system is linear, so a single bandsolver.solve gives the answer.
"""
import numpy as np

import bandsolver as bs


def system(nj, second_order):
    h = 1.0 / (nj - 1)
    x = np.linspace(0, 1, nj)
    A = np.zeros((nj, 1, 1)); B = np.zeros((nj, 1, 1)); D = np.zeros((nj, 1, 1)); G = np.zeros((nj, 1))
    i = np.arange(1, nj - 1)
    A[i, 0, 0] = D[i, 0, 0] = -1 / h**2
    B[i, 0, 0] = 2 / h**2 + 1
    G[i, 0] = (np.pi**2 + 1) * np.cos(np.pi * x[i])
    if second_order:
        B[0, 0, 0], D[0, 0, 0], X = -3 / (2 * h), 4 / (2 * h), np.array([[-1 / (2 * h)]])
        B[-1, 0, 0], A[-1, 0, 0], Y = 3 / (2 * h), -4 / (2 * h), np.array([[1 / (2 * h)]])
    else:
        B[0, 0, 0], D[0, 0, 0], X = -1 / h, 1 / h, None
        B[-1, 0, 0], A[-1, 0, 0], Y = 1 / h, -1 / h, None
    return x, (A, B, D, G, X, Y)


def run(backend="cpp", grids=(21, 41, 81, 161, 321)):
    out = {}
    for second_order in (True, False):
        rows = []
        for nj in grids:
            x, sys_ = system(nj, second_order)
            c = bs.solve(*sys_, backend=backend)
            rows.append(dict(nj=nj, error=np.abs(c[:, 0] - np.cos(np.pi * x)).max()))
        for a, b in zip(rows[:-1], rows[1:]):
            b["order"] = np.log2(a["error"] / b["error"])
        out["second_order_XY" if second_order else "first_order"] = rows
    return out


if __name__ == "__main__":
    for backend in bs.BACKENDS:
        for name, rows in run(backend).items():
            print(f"backend={backend} closure={name}")
            for row in rows:
                print("  nj={nj:4d} max error={error:.3e} order={o}".format(
                    o=f"{row['order']:.3f}" if "order" in row else "  -  ", **row))
