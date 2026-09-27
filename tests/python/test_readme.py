"""The README's Python quickstart must run and give the exact discrete solution."""
import pathlib
import re

import numpy as np


def test_readme_python_snippet():
    text = (pathlib.Path(__file__).resolve().parents[2] / "README.md").read_text(encoding="utf-8")
    code = re.search(r"```python\n(.*?)```", text, re.S).group(1)
    ns = {}
    exec(compile(code.replace("print(", "(lambda *a: None)("), "README", "exec"), ns)
    r, nj = ns["r"], ns["nj"]
    x = np.linspace(0, 1, nj)
    assert r.converged and np.abs(r.c[:, 0] - x * (1 - x) / 2).max() < 1e-14
