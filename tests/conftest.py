import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pymol  # noqa: E402

# -c: no GUI, -k: ignore ~/.pymolrc, -q: quiet
pymol.finish_launching(['pymol', '-ckq'])

from pymol import cmd  # noqa: E402


@pytest.fixture(autouse=True)
def fresh_pymol():
    cmd.reinitialize()
    yield cmd
