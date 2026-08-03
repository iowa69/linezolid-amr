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
    assert "Enterococcus_faecium" in result.output
    # Mtb was removed in v0.1.2 (domain-specific tools cover it)
    assert "Mycobacterium_tuberculosis" not in result.output


def test_version():
    result = CliRunner().invoke(main, ["--version"])
    assert result.exit_code == 0


def test_run_without_inputs_is_rejected():
    """Neither assembly nor reads means there is nothing to analyse."""
    result = CliRunner().invoke(main, ["run", "-o", "/tmp/does-not-matter"])
    assert result.exit_code != 0
    assert "nothing to analyse" in result.output.lower()


def test_r2_without_r1_is_rejected(tmp_path):
    r2 = tmp_path / "x_R2.fastq"
    r2.write_text("@r\nACGT\n+\n!!!!\n")
    result = CliRunner().invoke(main, ["run", "-2", str(r2), "-o", str(tmp_path / "out")])
    assert result.exit_code != 0


def test_skip_rrna23s_without_assembly_reports_the_real_problem(tmp_path):
    """The user should be told nothing will run, not how to configure a
    stage they just disabled."""
    r1 = tmp_path / "s_R1.fastq"
    r1.write_text("@r\nACGT\n+\n!!!!\n")
    result = CliRunner().invoke(
        main, ["run", "-1", str(r1), "-o", str(tmp_path / "out"), "--skip-rrna23s"]
    )
    assert result.exit_code != 0
    assert "nothing to do" in result.output.lower()
    assert "requires --organism" not in result.output


def test_evidence_tier_choice_is_validated(tmp_path):
    r1 = tmp_path / "s_R1.fastq"
    r1.write_text("@r\nACGT\n+\n!!!!\n")
    result = CliRunner().invoke(
        main,
        ["run", "-1", str(r1), "-O", "Enterococcus_faecium",
         "-o", str(tmp_path / "out"), "--evidence-tier", "made-up-tier"],
    )
    assert result.exit_code != 0
    assert "made-up-tier" in result.output
