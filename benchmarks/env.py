"""Record the machine, toolchain and package versions a benchmark run used."""
import json
import pathlib
import platform
import subprocess
import sys
from importlib import metadata


def _run(cmd):
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=20).stdout.strip().splitlines()[0]
    except Exception:
        return None


def cpu_name():
    if sys.platform == "darwin":
        return _run(["sysctl", "-n", "machdep.cpu.brand_string"])
    try:
        for line in pathlib.Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or None


def collect():
    root = pathlib.Path(__file__).resolve().parents[1]
    pkgs = {}
    for name in ("bandsolver", "numpy", "scipy", "scikit-sundae", "matplotlib"):
        try:
            pkgs[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            pkgs[name] = None
    return {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu": cpu_name(),
        "python": sys.version.split()[0],
        "packages": pkgs,
        "gfortran": _run(["gfortran", "--version"]),
        "cxx": _run(["c++", "--version"]),
        "git_commit": _run(["git", "-C", str(root), "rev-parse", "--short", "HEAD"]),
    }


def write(path):
    info = collect()
    pathlib.Path(path).write_text(json.dumps(info, indent=2) + "\n")
    return info


if __name__ == "__main__":
    out = pathlib.Path(__file__).resolve().parent / "results" / "env.json"
    print(json.dumps(write(out), indent=2))
