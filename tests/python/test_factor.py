"""P3: factor once, solve many, on both backends."""
import numpy as np
import pytest

import bandsolver as bs
from conftest import random_system


@pytest.mark.parametrize("xy", [False, True])
@pytest.mark.parametrize("n,nj", [(1, 3), (2, 3), (3, 4), (6, 40), (4, 300)])
def test_factor_matches_solve(rng, backend, xy, n, nj):
    A, B, D, G, X, Y = random_system(rng, n, nj, xy)
    f = bs.factor(A, B, D, X, Y, backend=backend)
    assert f.backend == backend and (f.nj, f.n) == (nj, n)
    for _ in range(3):
        G = rng.uniform(-1, 1, (nj, n))
        ref = bs.solve(A, B, D, G, X, Y)
        np.testing.assert_allclose(f.solve(G), ref, rtol=0, atol=1e-13 * np.abs(ref).max())


def test_backends_agree(rng):
    A, B, D, G, X, Y = random_system(rng, 5, 60, True)
    a = bs.factor(A, B, D, X, Y, backend="cpp").solve(G)
    b = bs.factor(A, B, D, X, Y, backend="fortran").solve(G)
    np.testing.assert_allclose(a, b, rtol=0, atol=1e-13 * np.abs(a).max())


def test_factor_errors(rng, backend):
    A, B, D, G, _, _ = random_system(rng, 2, 9, False)
    A[5] = 0
    B[5, 1] = 2 * B[5, 0]
    with pytest.raises(bs.SingularBlockError) as e:
        bs.factor(A, B, D, backend=backend)
    assert e.value.node == 5
    A, B, D, G, _, _ = random_system(rng, 2, 5, False)
    f = bs.factor(A, B, D, backend=backend)
    with pytest.raises(ValueError, match="shape"):
        f.solve(G[:, :1])
    G[1, 0] = np.nan
    with pytest.raises(bs.NonFiniteError):
        f.solve(G)
