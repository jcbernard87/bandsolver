"""Transient diffusion c_t = c_xx on [0,1], c(0,t) = c(1,t) = 0, c(x,0) = sin(pi x).

Exact solution exp(-pi^2 t) sin(pi x). Implicit Euler in time; each step is one call to
bandsolver.newton (the problem is linear, so Newton converges at iteration 2). This is
how a time-dependent PDE is driven with the library: the stepper lives in user code.
"""
import numpy as np

import bandsolver as bs


def step_fill(c_old, dt):
    def fill(c):
        nj = c.shape[0]
        u, h = c[:, 0], 1.0 / (nj - 1)
        A = np.zeros((nj, 1, 1)); B = np.zeros((nj, 1, 1)); D = np.zeros((nj, 1, 1)); G = np.zeros((nj, 1))
        i = np.arange(1, nj - 1)
        F = (u[i] - c_old[i, 0]) / dt - (u[i + 1] - 2 * u[i] + u[i - 1]) / h**2
        G[i, 0] = -F
        A[i, 0, 0] = D[i, 0, 0] = -1 / h**2
        B[i, 0, 0] = 1 / dt + 2 / h**2
        B[0, 0, 0] = B[-1, 0, 0] = 1
        G[0, 0], G[-1, 0] = -u[0], -u[-1]
        return A, B, D, G

    return fill


def simulate(nj, dt, t_end, backend="cpp"):
    x = np.linspace(0, 1, nj)
    c = np.sin(np.pi * x)[:, None]
    for _ in range(int(round(t_end / dt))):
        c = bs.newton(step_fill(c, dt), c, backend=backend).c
    return np.abs(c[:, 0] - np.exp(-np.pi**2 * t_end) * np.sin(np.pi * x)).max()


def run(backend="cpp", nj=801, t_end=0.1, dts=(0.01, 0.005, 0.0025, 0.00125)):
    rows = [dict(dt=dt, error=simulate(nj, dt, t_end, backend)) for dt in dts]
    for a, b in zip(rows[:-1], rows[1:]):
        b["order"] = np.log2(a["error"] / b["error"])
    return rows


if __name__ == "__main__":
    for backend in bs.BACKENDS:
        print(f"backend={backend}")
        for row in run(backend):
            print("  dt={dt:.5f} max error={error:.3e} order={o}".format(
                o=f"{row['order']:.3f}" if "order" in row else "  -  ", **row))
