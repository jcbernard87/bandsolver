"""Steady nonlinear diffusion -(D(c) c')' = f on [0,1], D = 1 + c^2, c(0) = c(1) = 0.

Manufactured solution c* = sin(pi x). Finite-volume (face-averaged D) discretization,
second order in space. Demonstrates the Newton driver and grid convergence.
"""
import numpy as np

import bandsolver as bs


def make_fill():
    def fill(c):
        nj = c.shape[0]
        u = c[:, 0]
        h = 1.0 / (nj - 1)
        x = np.linspace(0, 1, nj)
        s, sx, sxx = np.sin(np.pi * x), np.pi * np.cos(np.pi * x), -np.pi**2 * np.sin(np.pi * x)
        f = -(2 * s * sx**2 + (1 + s**2) * sxx)
        A = np.zeros((nj, 1, 1)); B = np.zeros((nj, 1, 1)); D = np.zeros((nj, 1, 1)); G = np.zeros((nj, 1))
        i = np.arange(1, nj - 1)
        cm, cp = (u[i] + u[i - 1]) / 2, (u[i] + u[i + 1]) / 2
        Dm, Dp = 1 + cm**2, 1 + cp**2
        G[i, 0] = (Dp * (u[i + 1] - u[i]) - Dm * (u[i] - u[i - 1])) / h**2 + f[i]
        A[i, 0, 0] = (cm * (u[i] - u[i - 1]) - Dm) / h**2
        D[i, 0, 0] = -(Dp + cp * (u[i + 1] - u[i])) / h**2
        B[i, 0, 0] = (Dp + Dm - cp * (u[i + 1] - u[i]) + cm * (u[i] - u[i - 1])) / h**2
        B[0, 0, 0] = B[-1, 0, 0] = 1
        G[0, 0], G[-1, 0] = -u[0], -u[-1]
        return A, B, D, G

    return fill


def run(backend="cpp", grids=(21, 41, 81, 161, 321)):
    rows = []
    for nj in grids:
        r = bs.newton(make_fill(), np.zeros((nj, 1)), backend=backend)
        x = np.linspace(0, 1, nj)
        err = np.abs(r.c[:, 0] - np.sin(np.pi * x)).max()
        rows.append(dict(nj=nj, iterations=r.iterations, error=err, last_step=r.step_norm[-1]))
    for a, b in zip(rows[:-1], rows[1:]):
        b["order"] = np.log2(a["error"] / b["error"])
    return rows


if __name__ == "__main__":
    for backend in bs.BACKENDS:
        print(f"backend={backend}")
        for row in run(backend):
            print("  nj={nj:4d} iterations={iterations} max error={error:.3e} order={o}".format(
                o=f"{row['order']:.3f}" if "order" in row else "  -  ", **row))
