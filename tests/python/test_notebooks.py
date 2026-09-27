"""The tutorial notebooks must execute top to bottom without errors.

Skipped when the notebook tooling is not installed (e.g. in the wheel tests), or when the
notebooks are not present next to the tests.
"""
import pathlib

import pytest

nbformat = pytest.importorskip("nbformat")
nbclient = pytest.importorskip("nbclient")
pytest.importorskip("matplotlib")
pytest.importorskip("ipykernel")

NOTEBOOK_DIR = pathlib.Path(__file__).resolve().parents[2] / "notebooks"
NOTEBOOKS = sorted(NOTEBOOK_DIR.glob("*.ipynb"))


@pytest.mark.skipif(not NOTEBOOKS, reason="notebooks directory not available")
@pytest.mark.parametrize("path", NOTEBOOKS, ids=lambda p: p.stem)
def test_notebook_executes(path):
    nb = nbformat.read(path, as_version=4)
    client = nbclient.NotebookClient(nb, timeout=600, kernel_name="python3",
                                     resources={"metadata": {"path": str(path.parent)}})
    client.execute()   # raises CellExecutionError on the first failing cell
