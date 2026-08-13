"""M0 smoke tests for the packaged Python module."""

from __future__ import annotations


def test_version_string() -> None:
    import samrust

    assert isinstance(samrust.__version__, str)
    assert samrust.__version__
    assert samrust.version() == samrust.__version__
