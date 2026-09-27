"""Bundle the Intel Fortran runtime DLLs into a Windows wheel (cibuildwheel repair step)."""
import os
import shutil
import subprocess
import sys

dest_dir, wheel = sys.argv[1], sys.argv[2]
ifx = shutil.which("ifx")
if ifx is None:
    sys.exit("ifx not found on PATH; cannot locate the Intel runtime DLLs")
paths = {os.path.dirname(ifx)}
root = os.environ.get("ONEAPI_ROOT") or os.environ.get("CMPLR_ROOT")
for cand in (os.path.join(os.path.dirname(ifx), "..", "redist", "intel64_win", "compiler"),
             os.path.join(root or "", "compiler", "latest", "bin")):
    if os.path.isdir(cand):
        paths.add(os.path.abspath(cand))
subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "delvewheel"])
cmd = [sys.executable, "-m", "delvewheel", "repair", "-w", dest_dir, "--add-path", os.pathsep.join(sorted(paths)), wheel]
print(" ".join(cmd), flush=True)
subprocess.check_call(cmd)
