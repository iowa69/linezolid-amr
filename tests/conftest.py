"""Test-suite wide fixtures.

The single most important thing here is neutralising the reference override
directory. `references.override_dir()` prefers
``$LINEZOLID_AMR_REFDIR`` / ``~/.local/share/linezolid-amr/references`` over the
data bundled in the package — correct at runtime, but catastrophic for testing:
a developer who has ever run ``fetch-references`` gets a green suite that
validates their local copy while the files actually shipped to users go
unchecked. That is exactly how a re-curation of the linezolid position list
once passed every test without ever reaching the package.

Tests must see what a fresh install sees.
"""

from __future__ import annotations

import pytest


@pytest.fixture(autouse=True)
def _use_bundled_references(monkeypatch, tmp_path_factory):
    """Point reference lookups at the package data for every test.

    Setting LINEZOLID_AMR_REFDIR to a directory that exists but is empty makes
    `_resolve` find nothing there and fall through to the bundled files, while
    also shadowing XDG_DATA_HOME and the home-directory default.
    """
    empty = tmp_path_factory.mktemp("no-reference-override")
    monkeypatch.setenv("LINEZOLID_AMR_REFDIR", str(empty))
    monkeypatch.setenv("XDG_DATA_HOME", str(empty))
    # MLST schemes resolve through a separate variable with the same hazard.
    monkeypatch.delenv("LINEZOLID_AMR_MLST_DIR", raising=False)
    yield
