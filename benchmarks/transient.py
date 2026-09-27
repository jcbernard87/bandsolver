"""Layer 2: a nonlinear transient DAE (binary electrolyte, notebook 2) with three solver stacks.

All three solve the *same* finite-volume discretization with nj nodes:

* bandsolver: Newton on each fixed time step, backward Euler (BE) or BDF2, stepped from
  user code, with the analytic Jacobian blocks;
* SUNDIALS IDA (scikit-sundae): variable-order adaptive BDF on the DAE in (c, phi), with a
  band linear solver and IDA's internal finite-difference Jacobian;
* SciPy solve_ivp(BDF): SciPy has no DAE support, so phi is eliminated exactly at every face
  (current conservation gives i = I at each face), leaving an ODE in c with a tridiagonal
  Jacobian sparsity. The discrete solution is identical to the DAE's.

Unknowns are interleaved y = [c_0, phi_0, c_1, phi_1, ...]; phi is dimensionless (F Phi / RT).
"""
import time

import numpy as np
import scipy.integrate
import scipy.sparse

import bandsolver as bs

FARADAY = 96485.33212
DP, DM = 1.0e-9, 2.0e-9          # cation / anion diffusivity, m^2/s
LENGTH = 100e-6                  # m
C0 = 100.0                       # mol/m^3
CURRENT = 100.0                  # A/m^2
T_END = 5.0                      # s (diffusion time L^2/D = 7.5 s)


class Model:
    def __init__(self, nj=81):
        self.nj = nj
        self.h = LENGTH / (nj - 1)
        self.V = np.full(nj, self.h)
        self.V[0] = self.V[-1] = self.h / 2
        self.wall = CURRENT / FARADAY

    # ---- spatial operator: flux differences R(c, phi), without accumulation ---------------
    def _face(self, c, phi, z, Dz):
        h = self.h
        cl, cr, pl, pr = c[:-1], c[1:], phi[:-1], phi[1:]
        cbar, dphi = 0.5 * (cl + cr), pr - pl
        N = -Dz * ((cr - cl) / h + z * cbar * dphi / h)
        return (N, -Dz * (-1 / h + z * dphi / (2 * h)), -Dz * (1 / h + z * dphi / (2 * h)),
                Dz * z * cbar / h, -Dz * z * cbar / h)

    def operator(self, c, phi, with_jacobian=False):
        """R (nj, 2) with rows (anion, cation); last cation row is the reference phi(L)."""
        nj = self.nj
        R = np.zeros((nj, 2))
        if with_jacobian:
            A = np.zeros((nj, 2, 2)); B = np.zeros((nj, 2, 2)); D = np.zeros((nj, 2, 2))
        for row, (z, Dz, wall) in enumerate([(-1, DM, 0.0), (+1, DP, self.wall)]):
            N, dcl, dcr, dpl, dpr = self._face(c, phi, z, Dz)
            R[:-1, row] += N
            R[1:, row] -= N
            R[0, row] -= wall
            R[-1, row] += wall
            if with_jacobian:
                B[:-1, row, 0] += dcl; B[:-1, row, 1] += dpl; D[:-1, row, 0] += dcr; D[:-1, row, 1] += dpr
                B[1:, row, 0] -= dcr; B[1:, row, 1] -= dpr; A[1:, row, 0] -= dcl; A[1:, row, 1] -= dpl
        R[-1, 1] = phi[-1]
        if with_jacobian:
            A[-1, 1, :] = 0; B[-1, 1, :] = [0, 1]; D[-1, 1, :] = 0
            return R, A, B, D
        return R

    def initial(self):
        y = np.zeros((self.nj, 2))
        y[:, 0] = C0
        return y

    # ---- bandsolver: fixed-step BE / BDF2 ------------------------------------------------------
    def _step_fill(self, a0, hist):
        """Implicit step with c_dot ~ a0*c + hist (BE: 1/dt, -c_old/dt; BDF2: 3/(2dt), ...)."""
        V = self.V

        def fill(u):
            c, phi = u[:, 0], u[:, 1]
            R, A, B, D = self.operator(c, phi, with_jacobian=True)
            acc = V * (a0 * c + hist)
            F = R.copy()
            F[:, 0] += acc
            F[:-1, 1] += acc[:-1]                    # last cation row is the reference
            B[:, 0, 0] += V * a0
            B[:-1, 1, 0] += V[:-1] * a0
            return A, B, D, -F
        return fill

    def run_bandsolver(self, dt, scheme="BE", backend="cpp"):
        steps = int(round(T_END / dt))
        u = self.initial()
        c_prev = None
        newton_iters = 0
        for k in range(steps):
            c_n = u[:, 0].copy()
            if scheme == "BE" or c_prev is None:        # BDF2 starts with one BE step
                fill = self._step_fill(1 / dt, -c_n / dt)
            else:
                fill = self._step_fill(1.5 / dt, (-4 * c_n + c_prev) / (2 * dt))
            r = bs.newton(fill, u, backend=backend)
            u, c_prev = r.c, c_n
            newton_iters += r.iterations
        return dict(y=u, steps=steps, fevals=newton_iters, jevals=newton_iters, newton_iters=newton_iters)

    # ---- SUNDIALS IDA on the DAE -----------------------------------------------------------
    def run_ida(self, rtol, max_steps=1_000_000):
        from sksundae.ida import IDA
        nj, V = self.nj, self.V

        def resfn(t, y, yp, res):
            c, phi = y[0::2], y[1::2]
            R = self.operator(c, phi)
            R[:, 0] += V * yp[0::2]
            R[:-1, 1] += V[:-1] * yp[:-2:2]
            res[:] = R.ravel()

        atol = np.empty(2 * nj)
        atol[0::2], atol[1::2] = rtol * C0, rtol * 1.0
        solver = IDA(resfn, linsolver="band", lband=3, uband=3, rtol=rtol, atol=atol,
                     algebraic_idx=np.arange(1, 2 * nj, 2), calc_initcond="yp0", calc_init_dt=1e-6,
                     max_num_steps=max_steps)
        y0 = self.initial().ravel()
        sol = solver.solve(np.array([0.0, T_END]), y0, np.zeros_like(y0))
        if not sol.success:
            raise RuntimeError(f"IDA failed: {sol.message}")
        return dict(y=sol.y[-1].reshape(nj, 2), steps=len(sol.t) - 1, fevals=sol.nfev, jevals=sol.njev,
                    newton_iters=None)

    # ---- SciPy BDF on the reduced ODE ------------------------------------------------------------
    def reduced_rhs(self, t, c):
        h = self.h
        dc = np.diff(c)
        cbar = 0.5 * (c[:-1] + c[1:])
        dphi = -(h * self.wall + (DP - DM) * dc) / ((DP + DM) * cbar)   # from i = I at each face
        N = -DM * (dc / h - cbar * dphi / h)                             # anion flux (z = -1)
        net = np.zeros(self.nj)
        net[:-1] += N
        net[1:] -= N
        return -net / self.V

    def phi_from_c(self, c):
        h = self.h
        dc = np.diff(c)
        cbar = 0.5 * (c[:-1] + c[1:])
        dphi = -(h * self.wall + (DP - DM) * dc) / ((DP + DM) * cbar)
        phi = np.concatenate([[0.0], np.cumsum(dphi)])
        return phi - phi[-1]                                             # reference phi(L) = 0

    def run_scipy(self, rtol):
        nj = self.nj
        sparsity = scipy.sparse.diags([1.0, 1.0, 1.0], [-1, 0, 1], shape=(nj, nj))
        sol = scipy.integrate.solve_ivp(self.reduced_rhs, (0.0, T_END), np.full(nj, C0), method="BDF",
                                        rtol=rtol, atol=rtol * C0, jac_sparsity=sparsity)
        if not sol.success:
            raise RuntimeError(f"SciPy BDF failed: {sol.message}")
        c = sol.y[:, -1]
        return dict(y=np.column_stack([c, self.phi_from_c(c)]), steps=len(sol.t) - 1, fevals=sol.nfev,
                    jevals=sol.njev, newton_iters=None)


def timed(f, repeats=3):
    best, out = np.inf, None
    for _ in range(repeats):
        t0 = time.perf_counter()
        out = f()
        best = min(best, time.perf_counter() - t0)
    out["wall_s"] = best
    return out


if __name__ == "__main__":
    # Consistency check: every stack converges to the same discrete solution.
    m = Model(81)
    ref = m.run_ida(1e-12)["y"]
    for name, run in [("IDA rtol=1e-10", lambda: m.run_ida(1e-10)),
                      ("SciPy BDF rtol=1e-10", lambda: m.run_scipy(1e-10)),
                      ("bandsolver BDF2 dt=1e-4", lambda: m.run_bandsolver(1e-4, "BDF2")),
                      ("bandsolver BE dt=1e-4", lambda: m.run_bandsolver(1e-4, "BE"))]:
        t0 = time.perf_counter(); r = run(); dt = time.perf_counter() - t0
        print(f"{name:26s} max|c - c_ref| = {np.abs(r['y'][:, 0] - ref[:, 0]).max():.2e} mol/m^3   "
              f"max|phi - phi_ref| = {np.abs(r['y'][:, 1] - ref[:, 1]).max():.2e}   ({dt:.2f} s, {r['steps']} steps)")
