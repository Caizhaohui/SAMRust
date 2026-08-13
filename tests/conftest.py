"""Shared pytest configuration.

``--require-oracles`` turns oracle-missing skips into failures. Use in CI /
release gates so a missing pysam/numpy/rubam/bcftools can never silently
green-light a parity run (v0.1.1 P3c).
"""

from __future__ import annotations

import pytest

_ORACLE_SKIP_MARKERS = ("could not import", "not on PATH")


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--require-oracles",
        action="store_true",
        default=False,
        help="fail tests that would skip because an oracle (pysam/numpy/rubam/bcftools) is missing",
    )


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item: pytest.Item, call: pytest.CallInfo):
    outcome = yield
    report = outcome.get_result()
    if (
        report.when == "call"
        and report.skipped
        and item.config.getoption("--require-oracles")
    ):
        reason = str(call.excinfo.value) if call.excinfo else ""
        longrepr = str(report.longrepr)
        text = f"{reason} {longrepr}"
        if any(marker in text for marker in _ORACLE_SKIP_MARKERS):
            report.outcome = "failed"
            report.longrepr = (
                f"{item.nodeid}: oracle missing but --require-oracles is set\n{text}"
            )
