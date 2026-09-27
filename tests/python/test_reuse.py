"""P4: Jacobian reuse (modified Newton) from Python, both backends."""
import numpy as np
import pytest

import bandsolver as bs
from test_fd import Problem, start


def test_reuse_same_root_fewer_factorizations(backend):
    p = Problem(3, 40, 5)
    full = bs.newton(p.fill, start(3, 40), backend=backend)
    reuse = bs.newton(p.fill, start(3, 40), jacobian_reuse=True, backend=backend)
    assert full.factorizations == full.iterations == full.jacobian_evaluations
    assert reuse.factorizations < reuse.iterations
    np.testing.assert_allclose(reuse.c, full.c, rtol=0, atol=1e-9)


def test_residual_callback_skips_blocks(backend):
    p = Problem(3, 40, 5)
    calls = {"fill": 0, "res": 0}

    def fill(c):
        calls["fill"] += 1
        return p.fill(c)

    def residual(c):
        calls["res"] += 1
        return p.residual(c)

    r = bs.newton(fill, start(3, 40), jacobian_reuse=True, residual=residual, backend=backend)
    assert r.converged
    assert calls["fill"] == r.factorizations == r.jacobian_evaluations
    assert calls["res"] == r.residual_evaluations == r.iterations - r.factorizations


def test_newton_fd_reuse(backend):
    p = Problem(3, 40, 5)
    a = bs.newton_fd(p.residual, start(3, 40), backend=backend)
    b = bs.newton_fd(p.residual, start(3, 40), jacobian_reuse=True, backend=backend)
    np.testing.assert_allclose(b.c, a.c, rtol=0, atol=1e-9)
    assert b.residual_evaluations == 10 * b.factorizations + (b.iterations - b.factorizations)
    assert b.residual_evaluations < a.residual_evaluations


def test_backends_agree_with_reuse():
    p = Problem(3, 40, 5)
    a = bs.newton(p.fill, start(3, 40), jacobian_reuse=True, residual=p.residual, backend="cpp")
    b = bs.newton(p.fill, start(3, 40), jacobian_reuse=True, residual=p.residual, backend="fortran")
    assert (a.iterations, a.factorizations) == (b.iterations, b.factorizations)
    np.testing.assert_allclose(a.c, b.c, rtol=0, atol=1e-12)


class Boom(Exception):
    pass


def test_residual_exception_propagates(backend):
    p = Problem(2, 10, 1)

    def residual(c):
        raise Boom("residual")

    with pytest.raises(Boom):
        bs.newton(p.fill, start(2, 10), jacobian_reuse=True, residual=residual, backend=backend)


def test_reuse_option_validation():
    p = Problem(2, 10, 1)
    with pytest.raises(ValueError):
        bs.newton(p.fill, start(2, 10), jacobian_reuse=True, reuse_max_iter=0)
    with pytest.raises(ValueError):
        bs.newton_fd(p.residual, start(2, 10), jacobian_reuse=True, reuse_contraction=0)
