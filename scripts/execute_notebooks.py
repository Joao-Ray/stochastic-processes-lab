"""Execute selected notebooks, or all notebooks, using this Python environment.

Run from the repository root: python scripts/execute_notebooks.py [notebook ...]
"""
from pathlib import Path
import sys

import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager


def main():
    root = Path(__file__).resolve().parents[1]
    paths = [Path(value).resolve() for value in sys.argv[1:]] if len(sys.argv) > 1 else sorted(
        (root / "notebooks").glob("*.ipynb"))
    for path in paths:
        notebook = nbformat.read(path, as_version=4)
        manager = KernelManager(kernel_name="python3")
        manager.kernel_spec.argv = [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"]
        print(f"Executing {path.name}", flush=True)
        NotebookClient(notebook, km=manager, timeout=300,
                       resources={"metadata": {"path": str(root)}}).execute()
        nbformat.write(notebook, path)
        for cell in notebook.cells:
            for output in cell.get("outputs", []):
                if output.output_type == "stream":
                    print(output.text, end="" if output.text.endswith("\n") else "\n")
        print(f"Saved {path.name}", flush=True)


if __name__ == "__main__":
    main()
