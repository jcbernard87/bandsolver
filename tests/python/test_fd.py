"""Finite-difference Jacobians from Python on both backends (F3)."""
import numpy as np
import pytest

import bandsolver as bs


class Problem:
    """F_j = Bm c_j + Am c_{j-1} + Dm c_{j+1} + 0.5 sin(c_j) * roll(c_j) - b (j+1)
    + [j=0] Xm c_2^2 + [j=nj-1] Ym c_{nj-3}^3, with its exact Jacobian."""

    def __init__(self, n, nj, seed):
        rng = np.random.default_rng(seed)
        self.n, self.nj = n, nj
        self.Am, self.Bm, self.Dm, self.Xm, self.Ym = (0.3 * rng.uniform(-1, 1, (n, n)) for _ in range(5))
        self.Bm += (3 + n) * np.eye(n)
        self.b = rng.uniform(-1, 1, n)

    def residual(self, c):
        nxt = np.roll(c, -1, axis=1)
        F = c @ self.Bm.T + 0.5 * np.sin(c) * nxt - np.outer(np.arange(1, self.nj + 1), self.b)
        F[1:] += c[:-1] @ self.Am.T
        F[:-1] += c[1:] @ self.Dm.T
        F[0] += self.Xm @ c[2] ** 2
        F[-1] += self.Ym @ c[-3] ** 3
        return F

    def fill(self, c):
        n, nj = self.n, self.nj
        A = np.broadcast_to(self.Am, (nj, n, n)).copy(); A[0] = 0
        D = np.broadcast_to(self.Dm, (nj, n, n)).copy(); D[-1] = 0
        B = np.broadcast_to(self.Bm, (nj, n, n)).copy()
        i = np.arange(n); ip = (i + 1) % n
        for j in range(nj):
            B[j, i, i] += 0.5 * np.cos(c[j, i]) * c[j, ip]
            B[j, i, ip] += 0.5 * np.sin(c[j, i])
        X = 2 * self.Xm * c[2]
        Y = 3 * self.Ym * c[-3] ** 2
        return A, B, D, -self.residual(c), X, Y


def start(n, nj):
    return (0.3 + 0.01 * np.arange(n * nj)).reshape(nj, n)


@pytest.mark.parametrize("n,nj", [(1, 3), (3, 3), (3, 4), (3, 10)])
def test_fd_jacobian_matches_analytic(backend, n, nj):
    p = Problem(n, nj, 11 + nj)
    c = start(n, nj)
    fd = bs.fd_jacobian(p.residual, c, backend=backend)
    ex = p.fill(c)
    scale = max(np.abs(e).max() for e in ex[:3])
    for a, e in zip(fd, ex):
        assert np.abs(a - e).max() / scale < 1e-6


def test_newton_fd_matches_analytic_and_counts(backend):
    p = Problem(3, 40, 5)
    ra = bs.newton(p.fill, start(3, 40), backend=backend)
    rf = bs.newton_fd(p.residual, start(3, 40), backend=backend)
    assert rf.converged and rf.iterations == ra.iterations
    np.testing.assert_allclose(rf.c, ra.c, rtol=0, atol=1e-10)
    assert rf.residual_evaluations == (3 * 3 + 1) * rf.iterations


def test_backends_agree():
    """Bit-identical on the reference toolchains (gfortran + clang/gcc). Across compilers
    (e.g. ifx + MSVC on Windows) the solutions differ by a few ulp, and forward differences
    amplify an ulp in F to ~eps/h ~ 1e-8 relative in the blocks, hence the tolerances."""
    p = Problem(3, 20, 7)
    a = bs.newton_fd(p.residual, start(3, 20), backend="cpp")
    b = bs.newton_fd(p.residual, start(3, 20), backend="fortran")
    assert a.iterations == b.iterations and a.residual_evaluations == b.residual_evaluations
    np.testing.assert_allclose(a.c, b.c, rtol=0, atol=1e-13 * np.abs(a.c).max())
    fa = bs.fd_jacobian(p.residual, start(3, 20), backend="cpp")
    fb = bs.fd_jacobian(p.residual, start(3, 20), backend="fortran")
    scale = max(np.abs(x).max() for x in fa[:3])
    for x, y in zip(fa, fb):
        np.testing.assert_allclose(x, y, rtol=0, atol=1e-6 * scale)


def test_check_jacobian_clean_and_planted(backend):
    p = Problem(3, 10, 9)
    c = start(3, 10)
    ok = bs.check_jacobian(p.fill, c, backend=backend)
    assert ok.max_error < 1e-4

    def bad_fill(x):
        A, B, D, G, X, Y = p.fill(x)
        D = D.copy(); D[4, 1, 2] += 0.5
        return A, B, D, G, X, Y

    bad = bs.check_jacobian(bad_fill, c, backend=backend)
    assert (bad.D.node, bad.D.row, bad.D.col) == (4, 1, 2) and bad.D.error > 1e-2
    assert max(bad.A.error, bad.B.error, bad.X.error, bad.Y.error) < 1e-4
    assert bad.worst()[0] == "D"


def test_check_jacobian_accepts_four_tuple_fill(backend):
    p = Problem(2, 6, 3)
    r = bs.check_jacobian(lambda x: p.fill(x)[:4], start(2, 6), backend=backend)
    # X/Y omitted by the fill but nonzero in the residual: reported as mismatches.
    assert r.X.error > 1e-2 and r.Y.error > 1e-2


class Boom(Exception):
    pass


def test_errors_propagate(backend):
    def residual(c):
        raise Boom("residual failed")

    with pytest.raises(Boom):
        bs.fd_jacobian(residual, np.zeros((5, 2)), backend=backend)
    with pytest.raises(Boom):
        bs.newton_fd(residual, np.zeros((5, 2)), backend=backend)
    with pytest.raises(Boom):
        bs.check_jacobian(lambda c: (_ for _ in ()).throw(Boom("fill failed")), np.zeros((5, 2)), backend=backend)
    with pytest.raises(ValueError, match="F has"):
        bs.fd_jacobian(lambda c: np.zeros(3), np.zeros((5, 2)), backend=backend)


def test_option_validation():
    p = Problem(2, 5, 1)
    with pytest.raises(ValueError):
        bs.fd_jacobian(p.residual, start(2, 5), rel_step=0)
    with pytest.raises(ValueError):
        bs.newton_fd(p.residual, start(2, 5), typical=-1)
    with pytest.raises(ValueError):
        bs.fd_jacobian(p.residual, np.zeros(5))
