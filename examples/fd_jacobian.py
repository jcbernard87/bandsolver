"""Residual-only solving and Jacobian checking with finite differences.

The coupled differential-algebraic system of ``reaction_diffusion.py`` (n = 3 unknowns per
node) is solved twice: with its hand-written Jacobian (``newton``) and from its residual
alone (``newton_fd``). The finite-difference Jacobian costs 3n + 1 = 10 residual calls per
iteration regardless of the grid size, and gives the same solution and convergence.
``check_jacobian`` then confirms the hand-written Jacobian and pinpoints a planted mistake.
"""
import numpy as np

import bandsolver as bs
import reaction_diffusion as rd


def residual(c):
    """F(c) for the reaction_diffusion system (the fill's G is -F)."""
    return -rd.make_fill()(c)[3]


def buggy_fill(c):
    """The hand-written Jacobian with a typical mistake: the sign of dF_v/du flipped."""
    A, B, D, G = rd.make_fill()(c)
    B = B.copy()
    B[1:-1, 1, 0] *= -1
    return A, B, D, G


def run(backend="cpp", grids=(21, 41, 81, 161, 321)):
    rows = []
    for nj in grids:
        x = np.linspace(0, 1, nj)
        c0 = np.zeros((nj, 3)); c0[:, 1] = 1
        ra = bs.newton(rd.make_fill(), c0, backend=backend)
        rf = bs.newton_fd(residual, c0, backend=backend)
        err = np.abs(rf.c - np.stack(rd.exact(x), axis=1)).max()
        rows.append(dict(nj=nj, it_analytic=ra.iterations, it_fd=rf.iterations, evals=rf.residual_evaluations,
                         diff=np.abs(rf.c - ra.c).max(), error=err))
    for a, b in zip(rows[:-1], rows[1:]):
        b["order"] = np.log2(a["error"] / b["error"])
    c = np.zeros((81, 3)); c[:, 1] = 1 + np.linspace(0, 1, 81) * 0.1; c[:, 0] = 0.5
    good = bs.check_jacobian(rd.make_fill(), c, backend=backend)
    bad = bs.check_jacobian(buggy_fill, c, backend=backend)
    return rows, good, bad


if __name__ == "__main__":
    for backend in bs.BACKENDS:
        rows, good, bad = run(backend)
        print(f"backend={backend}")
        for r in rows:
            print("  nj={nj:4d} newton {it_analytic} it | newton_fd {it_fd} it, {evals} residual calls | "
                  "|c_fd - c_analytic| = {diff:.1e} | error {error:.3e} order {o}".format(
                      o=f"{r['order']:.3f}" if "order" in r else "  -  ", **r))
        name, m = bad.worst()
        print(f"  check_jacobian: correct fill max error {good.max_error:.1e}; buggy fill worst = block {name} "
              f"node {m.node} row {m.row} col {m.col}, error {m.error:.2f} (user {m.user:.3g}, fd {m.fd:.3g})")
