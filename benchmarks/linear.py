"""Layer 1: linear block solves. BAND (Python and native) vs LAPACK banded, SuperLU, dense.

Every method solves the *same* Appendix C system (block tridiagonal plus X/Y endpoint
blocks). Only the factor+solve step is timed; converting the blocks into each library's
input format is excluded (and reported separately as `assembly_s`).

    python benchmarks/linear.py [--quick] [--native build-bench/benchmarks/native_bench]
"""
import argparse
import csv
import pathlib
import subprocess
import sys
import time

import numpy as np
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg

import bandsolver as bs

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import env  # noqa: E402

RESULTS = pathlib.Path(__file__).resolve().parent / "results"


def make_system(n, nj, seed):
    rng = np.random.default_rng(seed)
    A, B, D = (rng.uniform(-1, 1, (nj, n, n)) for _ in range(3))
    B += (3 * n + 3) * np.eye(n)
    G = rng.uniform(-1, 1, (nj, n))
    X, Y = rng.uniform(-1, 1, (n, n)), rng.uniform(-1, 1, (n, n))
    return A, B, D, G, X, Y


def to_coo(A, B, D, X, Y):
    """Global K as COO triplets (row, col, value)."""
    nj, n, _ = B.shape
    i, k = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    rows, cols, vals = [], [], []
    for off, blocks, js in ((0, B, np.arange(nj)), (-1, A, np.arange(1, nj)), (1, D, np.arange(nj - 1))):
        rows.append((js[:, None, None] * n + i).ravel())
        cols.append(((js + off)[:, None, None] * n + k).ravel())
        vals.append(blocks[js].ravel())
    rows += [i.ravel(), ((nj - 1) * n + i).ravel()]
    cols += [(2 * n + k).ravel(), ((nj - 3) * n + k).ravel()]
    vals += [X.ravel(), Y.ravel()]
    return np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)


def bandwidth(n):
    return 3 * n - 1          # X couples node 0 to node 2 (and Y the last node to nj-3)


def to_lapack_band(r, c, v, N, lu):
    ab = np.zeros((2 * lu + 1, N))
    np.add.at(ab, (lu + r - c, c), v)
    return ab


def backward_error(A, B, D, G, X, Y, dc):
    r = np.einsum("jik,jk->ji", B, dc) - G
    r[1:] += np.einsum("jik,jk->ji", A[1:], dc[:-1])
    r[:-1] += np.einsum("jik,jk->ji", D[:-1], dc[1:])
    r[0] += X @ dc[2]
    r[-1] += Y @ dc[-3]
    rows = np.abs(B).sum(2)
    rows[1:] += np.abs(A[1:]).sum(2)
    rows[:-1] += np.abs(D[:-1]).sum(2)
    rows[0] += np.abs(X).sum(1)
    rows[-1] += np.abs(Y).sum(1)
    return np.abs(r).max() / (rows.max() * np.abs(dc).max() + np.abs(G).max())


def time_call(f, min_total, min_reps=5):
    times = []
    while len(times) < min_reps or sum(times) < min_total:
        t0 = time.perf_counter()
        out = f()
        times.append(time.perf_counter() - t0)
        if len(times) >= 20000:
            break
    times.sort()
    return out, times[len(times) // 2], times[0], len(times)


def run_python(ns, njs, min_total, dense_limit=3000):
    rows = []
    for n in ns:
        for nj in njs:
            A, B, D, G, X, Y = make_system(n, nj, 1000 * n + nj)
            N, lu = n * nj, bandwidth(n)
            g = G.ravel()
            t0 = time.perf_counter()
            r, c, v = to_coo(A, B, D, X, Y)
            K = scipy.sparse.csc_matrix((v, (r, c)), shape=(N, N))
            t_csc = time.perf_counter() - t0
            t0 = time.perf_counter()
            ab = to_lapack_band(r, c, v, N, lu)
            t_band = time.perf_counter() - t0

            methods = {
                "band_cpp_py": (lambda: bs.solve(A, B, D, G, X, Y, backend="cpp"), 0.0,
                                nj * n * (n + 1) * 8),
                "band_fortran_py": (lambda: bs.solve(A, B, D, G, X, Y, backend="fortran"), 0.0,
                                    nj * n * (n + 1) * 8),
                "band_fortran_reference_py": (lambda: bs.solve(A, B, D, G, X, Y, backend="fortran",
                                                               kernel="reference"), 0.0, nj * n * (n + 1) * 8),
                "lapack_band": (lambda: scipy.linalg.solve_banded((lu, lu), ab, g, check_finite=False).reshape(nj, n),
                                t_band, (3 * lu + 1) * N * 8),
            }
            lu_nnz = {}
            for spec in ("NATURAL", "COLAMD"):
                def superlu(spec=spec):
                    lu_ = scipy.sparse.linalg.splu(K, permc_spec=spec)
                    lu_nnz[spec] = lu_.L.nnz + lu_.U.nnz
                    return lu_.solve(g).reshape(nj, n)
                methods[f"superlu_{spec.lower()}"] = (superlu, t_csc, spec)
            if N <= dense_limit:
                Kd = K.toarray()
                methods["dense_lapack"] = (lambda: np.linalg.solve(Kd, g).reshape(nj, n), 0.0, N * N * 8)

            for name, (f, t_asm, nbytes) in methods.items():
                dc, med, mn, reps = time_call(f, min_total)
                if isinstance(nbytes, str):             # SuperLU: L+U nonzeros (8-byte value + 4-byte index)
                    nbytes = lu_nnz[nbytes] * 12
                rows.append(dict(method=name, n=n, nj=nj, N=N, median_s=med, min_s=mn, reps=reps,
                                 backward_error=backward_error(A, B, D, G, X, Y, dc),
                                 factor_bytes=nbytes, assembly_s=t_asm))
            print(f"n={n:2d} nj={nj:5d} " + "  ".join(
                f"{r['method']}={r['median_s'] * 1e6:9.1f}us" for r in rows if r["n"] == n and r["nj"] == nj),
                flush=True)
    return rows


def run_native(exe, quick):
    out = subprocess.run([exe] + (["--quick"] if quick else []), capture_output=True, text=True, check=True).stdout
    rows = []
    for rec in csv.DictReader(out.splitlines()):
        n, nj = int(rec["n"]), int(rec["nj"])
        rows.append(dict(method=rec["backend"], n=n, nj=nj, N=n * nj, median_s=float(rec["median_s"]),
                         min_s=float(rec["min_s"]), reps=int(rec["reps"]),
                         backward_error=float(rec["backward_error"]), factor_bytes=nj * n * (n + 1) * 8,
                         assembly_s=0.0))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="tiny sweep for smoke testing")
    ap.add_argument("--native", default=None, help="path to the native_bench executable")
    ap.add_argument("--out", default=str(RESULTS / "linear.csv"))
    a = ap.parse_args()
    ns = [1, 3] if a.quick else [1, 3, 5, 10, 20, 30]
    njs = [25, 50] if a.quick else [25, 50, 100, 200, 500, 1000, 2000]
    min_total = 0.01 if a.quick else 0.2
    rows = run_python(ns, njs, min_total)
    if a.native:
        rows += run_native(a.native, a.quick)
    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    env.write(out.parent / "env.json")
    worst = max(r["backward_error"] for r in rows)
    print(f"wrote {len(rows)} rows to {out}; worst backward error {worst:.1e}")
    if worst > 1e-12:
        sys.exit("backward error check FAILED")


if __name__ == "__main__":
    main()
