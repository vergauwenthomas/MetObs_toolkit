"""Run sanity checks embedded in module ``__main__`` blocks.

Any ``.py`` file under ``src/metobs_toolkit/`` that contains an
``if __name__ == "__main__":`` guard is automatically discovered and
executed as a subprocess.  A non-zero exit code makes the test fail.

This lets developers add quick self-tests inside a module while
keeping them part of the CI/CD pipeline.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

# --------------------------------------------------------------------------- #
# Discovery
# --------------------------------------------------------------------------- #

_TOOLKIT_ROOT = Path(__file__).resolve().parent.parent / "src" / "metobs_toolkit"


def _discover_main_modules() -> list[Path]:
    """Return every ``.py`` file under the toolkit that has a __main__ block."""
    modules: list[Path] = []
    for pyfile in sorted(_TOOLKIT_ROOT.rglob("*.py")):
        try:
            text = pyfile.read_text(encoding="utf-8")
        except Exception:
            continue
        if 'if __name__ == "__main__"' in text or "if __name__ == '__main__'" in text:
            modules.append(pyfile)
    return modules


_DISCOVERED = _discover_main_modules()


# --------------------------------------------------------------------------- #
# Parametrised test
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "module_path",
    _DISCOVERED,
    ids=[str(p.relative_to(_TOOLKIT_ROOT)) for p in _DISCOVERED],
)
def test_module_main_block(module_path: Path) -> None:
    """Execute a module's ``__main__`` block and assert it exits cleanly."""
    result = subprocess.run(
        [sys.executable, str(module_path)],
        capture_output=True,
        text=True,
        timeout=120,
    )
    if result.returncode != 0:
        # Show both stdout and stderr so failures are easy to diagnose
        msg = (
            f"Module {module_path.relative_to(_TOOLKIT_ROOT)} exited with "
            f"code {result.returncode}.\n"
            f"--- stdout ---\n{result.stdout}\n"
            f"--- stderr ---\n{result.stderr}"
        )
        pytest.fail(msg)
