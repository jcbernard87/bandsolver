"""Figures for the benchmark results (reads benchmarks/results/*.csv, writes PNGs there)."""
import csv
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RESULTS = pathlib.Path(__file__).resolve().parent / "results"
BLUE, ORANGE, AQUA, YELLOW, MAGENTA = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb",
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2, "lines.markersize": 5,
    "legend.frameon": False, "figure.dpi": 120, "font.size": 10, "axes.titlesize": 11,
})

LINEAR_SERIES = [  # (method, label, colour, linestyle) in fixed order
    ("band_cpp_py", "BAND C++ (from Python)", BLUE, "-"),
    ("band_cpp_native", "BAND C++ (native)", BLUE, "--"),
    ("band_fortran_py", "BAND Fortran (from Python)", ORANGE, "-"),
    ("lapack_band", "LAPACK banded (dgbsv)", AQUA, "-"),
    ("superlu_natural", "SciPy SuperLU", YELLOW, "-"),
    ("dense_lapack", "dense LAPACK (dgesv)", MAGENTA, "-"),
]


def read(name):
    return list(csv.DictReader((RESULTS / name).open()))


def plot_linear():
    rows = read("linear.csv")
    t = {(r["method"], int(r["n"]), int(r["nj"])): float(r["median_s"]) for r in rows}
    njs = sorted({int(r["nj"]) for r in rows})
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.9), sharey=False)
    for ax, n in zip(axes, [3, 10, 30]):
        for method, label, color, ls in LINEAR_SERIES:
            pts = [(nj, t[(method, n, nj)]) for nj in njs if (method, n, nj) in t]
            if pts:
                ax.loglog(*zip(*pts), ls, marker="o", color=color, label=label)
        ax.set_title(f"n = {n} unknowns per node")
        ax.set_xlabel("nodes nj"); ax.set_xticks(njs, [str(j) for j in njs]); ax.minorticks_off()
    axes[0].set_ylabel("median solve time (s)")
    axes[0].legend(fontsize=8, loc="upper left")
    fig.suptitle("Layer 1: one block-tridiagonal solve (with X/Y endpoint blocks)", y=1.0)
    fig.tight_layout(); fig.savefig(RESULTS / "linear_time.png", bbox_inches="tight"); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.7), sharey=True)
    for ax, other, title in [(axes[0], "lapack_band", "vs LAPACK banded"), (axes[1], "superlu_natural", "vs SciPy SuperLU")]:
        for n, color in zip([1, 3, 10, 20, 30], RAMP):
            pts = [(nj, t[(other, n, nj)] / t[("band_cpp_py", n, nj)]) for nj in njs]
            ax.semilogx(*zip(*pts), "o-", color=color, label=f"n = {n}")
        ax.axhline(1, color=INK2, linewidth=1, linestyle=":")
        ax.set_title(f"Speed-up of BAND C++ {title}")
        ax.set_xlabel("nodes nj"); ax.set_xticks(njs, [str(j) for j in njs]); ax.minorticks_off()
    axes[0].set_ylabel("time(other) / time(BAND)")
    axes[1].legend(fontsize=8, loc="upper right")
    fig.tight_layout(); fig.savefig(RESULTS / "linear_speedup.png", bbox_inches="tight"); plt.close(fig)
    return ["linear_time.png", "linear_speedup.png"]


if __name__ == "__main__":
    which = sys.argv[1:] or ["linear"]
    for w in which:
        print(globals()[f"plot_{w}"]())
