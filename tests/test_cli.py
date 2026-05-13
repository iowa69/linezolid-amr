"""Smoke tests for the CLI."""

from __future__ import annotations

from click.testing import CliRunner

from linezolid_amr.cli import main


def test_help():
    result = CliRunner().invoke(main, ["--help"])
    assert result.exit_code == 0
    assert "linezolid-amr" in result.output.lower() or "Usage:" in result.output


def test_list_organisms():
    result = CliRunner().invoke(main, ["list-organisms"])
    assert result.exit_code == 0
    assert "Staphylococcus_aureus" in result.output
    assert "Mycobacterium_tuberculosis" in result.output


def test_version():
    result = CliRunner().invoke(main, ["--version"])
    assert result.exit_code == 0
