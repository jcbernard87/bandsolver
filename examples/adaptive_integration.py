"""Adaptive time integration of a coupled DAE with bandsolver.integrate.

The binary-electrolyte model of notebook 2 (Nernst-Planck transport of concentration c and
potential phi between two electrodes at constant current) written as a DAE
F(t, c, cdot) = 0 on a finite-volume mesh. phi is algebraic (no time derivative), so it is
flagged in the `algebraic` mask. The same model is integrated

* with a fixed step (adaptive=False, the previous "user-stepping" standard), and
* adaptively (adaptive=True), with and without Jacobian reuse,

and compared with the analytic Fourier-series solution of the salt concentration.
"""
import time

import numpy as np

import bandsolver as bs

F_, DP, DM = 96485.33212, 1.0e-9, 2.0e-9
L, C0, I = 100e-6, 100.0, 100.0
D_SALT = 2 * DP * DM / (DP + DM)
Q = (1 - DP / (DP + DM)) * I / F_


def c_analytic(x, t, terms=400):
    n = np.arange(1, 2 * terms, 2)[:, None]
    return (C0 + (Q / D_SALT) * (L / 2 - x)
            - np.sum(4 * Q * L / (D_SALT * n**2 * np.pi**2) * np.cos(n * np.pi * x / L)
                     * np.exp(-n**2 * np.pi**2 * D_SALT * t / L**2), axis=0))


class BinaryElectrolyte:
    def __init__(self, nj=81):
        self.nj, self.h = nj, L / (nj - 1)
        self.V = np.full(nj, self.h)
        self.V[0] = self.V[-1] = self.h / 2
        self.x = np.linspace(0, L, nj)

    def _faces(self, c, phi, z, Dz):
        h = self.h
        cbar, dphi = 0.5 * (c[:-1] + c[1:]), np.diff(phi)
        N = -Dz * (np.diff(c) / h + z * cbar * dphi / h)
        dcl = -Dz * (-1 / h + z * dphi / (2 * h))
        dcr = -Dz * (1 / h + z * dphi / (2 * h))
        dpl = Dz * z * cbar / h
        return N, dcl, dcr, dpl, -dpl

    def residual(self, t, u, udot):
        c, phi, cdot = u[:, 0], u[:, 1], udot[:, 0]
        R = np.zeros((self.nj, 2))
        for row, (z, Dz, wall) in enumerate([(-1, DM, 0.0), (+1, DP, I / F_)]):
            N = self._faces(c, phi, z, Dz)[0]
            R[:, row] = self.V * cdot
            R[:-1, row] += N
            R[1:, row] -= N
            R[0, row] -= wall
            R[-1, row] += wall
        R[-1, 1] = phi[-1]                      # reference potential (algebraic row)
        return R

    def jacobian(self, t, u, udot, alpha):
        nj = self.nj
        c, phi = u[:, 0], u[:, 1]
        A = np.zeros((nj, 2, 2)); B = np.zeros((nj, 2, 2)); D = np.zeros((nj, 2, 2))
        for row, (z, Dz) in enumerate([(-1, DM), (+1, DP)]):
            _, dcl, dcr, dpl, dpr = self._faces(c, phi, z, Dz)
            B[:, row, 0] += alpha * self.V
            B[:-1, row, 0] += dcl; B[:-1, row, 1] += dpl; D[:-1, row, 0] += dcr; D[:-1, row, 1] += dpr
            B[1:, row, 0] -= dcr; B[1:, row, 1] -= dpr; A[1:, row, 0] -= dcl; A[1:, row, 1] -= dpl
        A[-1, 1, :] = 0; B[-1, 1, :] = [0.0, 1.0]
        return A, B, D

    def initial(self):
        u = np.zeros((self.nj, 2))
        u[:, 0] = C0
        return u

    def algebraic(self):
        mask = np.zeros((self.nj, 2), dtype=bool)
        mask[:, 1] = True                        # phi has no time derivative
        return mask


def run(t_end=5.0, **opts):
    m = BinaryElectrolyte()
    t0 = time.perf_counter()
    r = bs.integrate(m.residual, m.initial(), [t_end], jacobian=m.jacobian, algebraic=m.algebraic(), **opts)
    wall = time.perf_counter() - t0
    err = np.abs(r.y[-1][:, 0] - c_analytic(m.x, t_end)).max()
    return r, wall, err


if __name__ == "__main__":
    cases = [
        ("fixed step BDF2, dt = 0.01 s (no reuse)", dict(adaptive=False, dt=0.01, jacobian_reuse=False)),
        ("fixed step BDF2, dt = 0.01 s (reuse)", dict(adaptive=False, dt=0.01)),
        ("adaptive BDF1-2, rtol = 1e-6 (no reuse)", dict(rtol=1e-6, atol=1e-6, jacobian_reuse=False)),
        ("adaptive BDF1-2, rtol = 1e-6 (reuse)", dict(rtol=1e-6, atol=1e-6)),
    ]
    for label, opts in cases:
        r, wall, err = run(**opts)
        s = r.stats
        print(f"{label:42s} max |c - analytic| = {err:.2e} mol/m^3 | {wall * 1e3:6.1f} ms | "
              f"{s['steps']} steps, {s['factorizations']} factorizations, {s['residual_evaluations']} residual calls")
