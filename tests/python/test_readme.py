"""The README's Python snippets must run (in order, sharing one namespace) and be correct."""
import pathlib
import re

import numpy as np


def test_readme_python_snippets():
    text = (pathlib.Path(__file__).resolve().parents[2] / "README.md").read_text(encoding="utf-8")
    blocks = re.findall(r"```python\n(.*?)```", text, re.S)
    assert len(blocks) >= 2
    ns = {}
    for i, code in enumerate(blocks):
        exec(compile(code.replace("print(", "(lambda *a: None)("), f"README[{i}]", "exec"), ns)
    r, nj = ns["r"], ns["nj"]
    x = np.linspace(0, 1, nj)
    assert r.converged and np.abs(r.c[:, 0] - x * (1 - x) / 2).max() < 1e-12
    assert r.residual_evaluations == (3 * 1 + 1) * r.iterations  # final r is from newton_fd, n = 1
