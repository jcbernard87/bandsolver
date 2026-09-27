import numpy as np
import pytest

import bandsolver as bs
from conftest import random_system


def bvp_fill(nonlinear=True):
    """-(D(c) c')' = f on [0,1], c(0)=c(1)=0, manufactured c* = sin(pi x), vectorized."""
    k = 1.0 if nonlinear else 0.0

    def fill(c):
        nj = c.shape[0]
        u = c[:, 0]
        h = 1.0 / (nj - 1)
        x = np.linspace(0, 1, nj)
        s, sx, sxx = np.sin(np.pi * x), np.pi * np.cos(np.pi * x), -np.pi**2 * np.sin(np.pi * x)
        f = -(2 * k * s * sx**2 + (1 + k * s**2) * sxx)
        A = np.zeros((nj, 1, 1)); B = np.zeros((nj, 1, 1)); D = np.zeros((nj, 1, 1)); G = np.zeros((nj, 1))
        i = np.arange(1, nj - 1)
        cm, cp = (u[i] + u[i - 1]) / 2, (u[i] + u[i + 1]) / 2
        Dm, Dp = 1 + k * cm**2, 1 + k * cp**2
        dDm, dDp = 2 * k * cm, 2 * k * cp
        F = -(Dp * (u[i + 1] - u[i]) - Dm * (u[i] - u[i - 1])) / h**2 - f[i]
        G[i, 0] = -F
        A[i, 0, 0] = (0.5 * dDm * (u[i] - u[i - 1]) - Dm) / h**2
        D[i, 0, 0] = -(Dp + 0.5 * dDp * (u[i + 1] - u[i])) / h**2
        B[i, 0, 0] = (Dp + Dm - 0.5 * dDp * (u[i + 1] - u[i]) + 0.5 * dDm * (u[i] - u[i - 1])) / h**2
        B[0, 0, 0] = B[-1, 0, 0] = 1
        G[0, 0], G[-1, 0] = -u[0], -u[-1]
        return A, B, D, G

    return fill


def test_nonlinear_bvp_quadratic_convergence(backend):
    nj = 101
    r = bs.newton(bvp_fill(), np.zeros((nj, 1)), backend=backend)
    assert r.converged and r.iterations <= 8
    s = r.step_norm
    for a, b in zip(s[:-1], s[1:]):
        if a < 0.1 and b > 1e-13:
            assert b <= 10 * a * a
    x = np.linspace(0, 1, nj)
    assert np.abs(r.c[:, 0] - np.sin(np.pi * x)).max() < 1e-3
    assert len(r.residual_norm) == r.iterations == len(r.step_norm)


def test_backends_identical_iterates():
    rc = bs.newton(bvp_fill(), np.zeros((101, 1)), backend="cpp")
    rf = bs.newton(bvp_fill(), np.zeros((101, 1)), backend="fortran")
    assert rc.iterations == rf.iterations
    np.testing.assert_allclose(rc.c, rf.c, rtol=0, atol=1e-13)
    np.testing.assert_allclose(rc.step_norm, rf.step_norm, rtol=1e-12)


def test_linear_system_with_xy_matches_solve(rng, backend):
    A, B, D, G, X, Y = random_system(rng, 3, 12, True)
    ref = bs.solve(A, B, D, G, X, Y)

    def fill(c):  # constant Jacobian, residual K c - G
        return A, B, D, G - _apply(A, B, D, X, Y, c), X, Y

    r = bs.newton(fill, np.zeros_like(G), backend=backend)
    assert r.converged and r.iterations == 2
    np.testing.assert_allclose(r.c, ref, rtol=0, atol=1e-13)


def _apply(A, B, D, X, Y, c):
    out = np.einsum("jik,jk->ji", B, c)
    out[1:] += np.einsum("jik,jk->ji", A[1:], c[:-1])
    out[:-1] += np.einsum("jik,jk->ji", D[:-1], c[1:])
    out[0] += X @ c[2]
    out[-1] += Y @ c[-3]
    return out


class Boom(Exception):
    pass


def test_callback_exception_propagates(backend):
    def fill(c):
        raise Boom("from fill")

    with pytest.raises(Boom, match="from fill"):
        bs.newton(fill, np.zeros((5, 1)), backend=backend)


def test_bad_fill_output(backend):
    with pytest.raises(ValueError, match="return"):
        bs.newton(lambda c: (1, 2, 3), np.zeros((5, 1)), backend=backend)
    with pytest.raises(ValueError, match="B has"):
        bs.newton(lambda c: (np.zeros((5, 1, 1)), np.zeros(3), np.zeros((5, 1, 1)), np.zeros((5, 1))),
                  np.zeros((5, 1)), backend=backend)


def test_not_converged_and_one_step(backend):
    with pytest.raises(bs.NotConvergedError) as e:
        bs.newton(bvp_fill(), np.zeros((51, 1)), max_iter=2, backend=backend)
    assert e.value.result.iterations == 2 and not e.value.result.converged
    r = bs.newton(bvp_fill(), np.zeros((51, 1)), max_iter=1, require_convergence=False, backend=backend)
    assert r.iterations == 1 and not r.converged


def test_singular_jacobian(backend):
    def fill(c):
        nj = c.shape[0]
        B = np.ones((nj, 1, 1)); B[3] = 0
        return np.zeros((nj, 1, 1)), B, np.zeros((nj, 1, 1)), np.ones((nj, 1))

    with pytest.raises(bs.SingularBlockError) as e:
        bs.newton(fill, np.zeros((6, 1)), backend=backend)
    assert e.value.node == 3


def test_option_validation():
    with pytest.raises(ValueError):
        bs.newton(bvp_fill(), np.zeros((5, 1)), damping=0)
    with pytest.raises(ValueError):
        bs.newton(bvp_fill(), np.zeros((5, 1)), rtol=0, atol=0)
    with pytest.raises(ValueError):
        bs.newton(bvp_fill(), np.zeros(5))
