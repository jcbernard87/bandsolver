"""Newman BAND block-tridiagonal solver (Electrochemical Systems, Appendix C).

Two interchangeable native backends are available: ``"cpp"`` (C++17 core) and
``"fortran"`` (Fortran 2008 library via its C ABI).

System (0-based nodes, ``nj >= 3`` nodes, ``n`` unknowns per node)::

    j = 0:        B[0] dc[0] + D[0] dc[1] + X dc[2]                 = G[0]
    0 < j < nj-1: A[j] dc[j-1] + B[j] dc[j] + D[j] dc[j+1]          = G[j]
    j = nj-1:     Y dc[nj-3] + A[nj-1] dc[nj-2] + B[nj-1] dc[nj-1]  = G[nj-1]

Arrays: ``A, B, D`` have shape ``(nj, n, n)``; ``G`` and ``dc`` have shape ``(nj, n)``;
``X, Y`` have shape ``(n, n)``. ``A[0]`` and ``D[nj-1]`` are ignored.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence

import numpy as np

from . import _core

__all__ = [
    "BACKENDS",
    "BandError",
    "SingularBlockError",
    "NonFiniteError",
    "NotConvergedError",
    "NewtonResult",
    "solve",
    "newton",
]

BACKENDS = ("cpp", "fortran")
_PIVOTS = {"partial": 0, "legacy": 1}

_OK, _SINGULAR, _INVALID, _NOT_CONVERGED, _NON_FINITE, _CALLBACK = range(6)


class BandError(RuntimeError):
    """Base class for solver failures. ``status`` is the numeric status code."""

    def __init__(self, message: str, status: int, node: int = -1):
        super().__init__(message)
        self.status = status
        self.node = node


class SingularBlockError(BandError):
    """A node's pivot block was (numerically) singular; ``node`` is 0-based."""


class NonFiniteError(BandError):
    """NaN or Inf in the inputs or the solution."""


class NotConvergedError(BandError):
    """Newton reached ``max_iter`` without converging; ``result`` holds the last iterate."""

    def __init__(self, message: str, result: "NewtonResult"):
        super().__init__(message, _NOT_CONVERGED)
        self.result = result


@dataclass
class NewtonResult:
    c: np.ndarray
    converged: bool
    iterations: int
    status: int
    update_norm: np.ndarray = field(repr=False)
    step_norm: np.ndarray = field(repr=False)
    residual_norm: np.ndarray = field(repr=False)


def _raise_for(status: int, node: int, what: str) -> None:
    if status == _OK:
        return
    if status == _SINGULAR:
        raise SingularBlockError(f"{what}: singular pivot block at node {node}", status, node)
    if status == _NON_FINITE:
        raise NonFiniteError(f"{what}: non-finite value in inputs or solution", status)
    if status == _INVALID:
        raise ValueError(f"{what}: invalid argument")
    raise BandError(f"{what}: failed with status {status}", status, node)


def _check_backend(backend: str) -> None:
    if backend not in BACKENDS:
        raise ValueError(f"backend must be one of {BACKENDS}, got {backend!r}")


def _pivot_code(pivot: str) -> int:
    try:
        return _PIVOTS[pivot]
    except KeyError:
        raise ValueError(f"pivot must be one of {tuple(_PIVOTS)}, got {pivot!r}") from None


def _as_blocks(a, name: str, shape=None) -> np.ndarray:
    arr = np.ascontiguousarray(a, dtype=np.float64)
    if shape is not None and arr.shape != shape:
        raise ValueError(f"{name} must have shape {shape}, got {arr.shape}")
    return arr


def solve(A, B, D, G, X=None, Y=None, *, pivot: str = "partial", backend: str = "cpp") -> np.ndarray:
    """Solve the block system and return ``dc`` with shape ``(nj, n)``.

    Raises :class:`SingularBlockError`, :class:`NonFiniteError` or ``ValueError``.
    """
    _check_backend(backend)
    B = _as_blocks(B, "B")
    if B.ndim != 3 or B.shape[1] != B.shape[2]:
        raise ValueError(f"B must have shape (nj, n, n), got {B.shape}")
    nj, n = B.shape[0], B.shape[1]
    if nj < 3:
        raise ValueError(f"need at least 3 nodes, got nj={nj}")
    A = _as_blocks(A, "A", (nj, n, n))
    D = _as_blocks(D, "D", (nj, n, n))
    G = _as_blocks(G, "G", (nj, n))
    X = None if X is None else _as_blocks(X, "X", (n, n))
    Y = None if Y is None else _as_blocks(Y, "Y", (n, n))
    dc, status, node, _ = _core.solve(n, nj, A, B, D, G, X, Y, _pivot_code(pivot), backend)
    _raise_for(status, node, "solve")
    return dc


FillResult = Sequence[Optional[np.ndarray]]


def newton(
    fill: Callable[[np.ndarray], FillResult],
    c0,
    *,
    rtol: float = 1e-10,
    atol: float = 1e-12,
    damping: float = 1.0,
    max_iter: int = 50,
    pivot: str = "partial",
    require_convergence: bool = True,
    backend: str = "cpp",
) -> NewtonResult:
    """Newton iteration ``c <- c + damping * dc`` with ``K(c) dc = G(c)``.

    ``fill(c)`` receives the current state (shape ``(nj, n)``, a copy) and returns
    ``(A, B, D, G)`` or ``(A, B, D, G, X, Y)`` where ``G = -F(c)`` is the negative
    residual and ``A, B, D, X, Y`` are the Jacobian blocks. Convergence: after the
    update, ``max |dc| / (atol + rtol |c|) <= 1``. ``max_iter=1`` with
    ``require_convergence=False`` reproduces the archival one-correction-per-step use.

    Exceptions raised inside ``fill`` propagate unchanged. Raises
    :class:`NotConvergedError` (with ``.result``) when ``require_convergence`` is set
    and ``max_iter`` is exhausted.
    """
    _check_backend(backend)
    c0 = np.ascontiguousarray(c0, dtype=np.float64)
    if c0.ndim != 2:
        raise ValueError(f"c0 must have shape (nj, n), got {c0.shape}")
    nj, n = c0.shape
    if nj < 3 or n < 1:
        raise ValueError(f"need nj >= 3 and n >= 1, got shape {c0.shape}")
    if not (0 < damping <= 1):
        raise ValueError("damping must be in (0, 1]")
    if rtol < 0 or atol < 0 or (rtol == 0 and atol == 0):
        raise ValueError("rtol and atol must be >= 0 and not both zero")
    r = _core.newton(n, nj, fill, c0, rtol, atol, damping, int(max_iter), _pivot_code(pivot),
                     bool(require_convergence), backend)
    result = NewtonResult(
        c=r["c"],
        converged=bool(r["converged"]),
        iterations=int(r["iterations"]),
        status=int(r["status"]),
        update_norm=np.asarray(r["update_norm"]),
        step_norm=np.asarray(r["step_norm"]),
        residual_norm=np.asarray(r["residual_norm"]),
    )
    if result.status == _NOT_CONVERGED:
        raise NotConvergedError(
            f"newton: not converged after {result.iterations} iterations "
            f"(last scaled update {result.update_norm[-1]:.3e})", result)
    _raise_for(result.status, int(r["fail_node"]), "newton")
    return result
