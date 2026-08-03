"""CLI-level integration tests: simulated reads in, full report chain out.

The unit tests check each stage in isolation; these check that a
heteroresistant sample actually arrives at the user as a heteroresistant
result, through every output file the tool writes.
"""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from linezolid_amr import references as ref_mod
from linezolid_amr.cli import main
from tests.simulate import read_fasta_single, simulate_reads

pytestmark = pytest.mark.integration

requires_tools = pytest.mark.skipif(
    shutil.which("minimap2") is None or shutil.which("samtools") is None
    or shutil.which("bcftools") is None,
    reason="minimap2/samtools/bcftools not on PATH",
)

ORGANISM = "Enterococcus_faecium"
G2576_SPECIES_POS = 2590


def _make_reads(tmp_path: Path, sample: str, af: float, depth: int = 250, seed: int = 42):
    _name, seq = read_fasta_single(ref_mod.organism_fasta_path(ORGANISM))
    r1 = tmp_path / f"{sample}_R1.fastq.gz"
    r2 = tmp_path / f"{sample}_R2.fastq.gz"
    simulate_reads(
        seq, r1, r2,
        variants={G2576_SPECIES_POS: ("T", af)} if af else {},
        depth=depth, seed=seed,
    )
    return r1, r2


def _run_cli(args):
    result = CliRunner().invoke(main, args)
    if result.exit_code != 0:
        raise AssertionError(
            f"CLI failed (exit {result.exit_code}): {result.output}\n{result.exception}"
        )
    return result


@requires_tools
def test_heteroresistant_sample_reported_end_to_end(tmp_path):
    """A 2-of-6-operon G2576T must surface as POSITIVE + heteroresistant everywhere."""
    r1, r2 = _make_reads(tmp_path, "HET001", af=1 / 3)
    out = tmp_path / "out"
    result = _run_cli([
        "run", "-1", str(r1), "-2", str(r2),
        "-O", ORGANISM, "-o", str(out), "-t", "2",
    ])
    assert "POSITIVE" in result.output

    report = json.loads((out / "HET001.linezolid_amr.json").read_text())
    summary = report["summary"]
    assert summary["linezolid_resistance_call"] is True
    assert summary["heteroresistance_detected"] is True
    assert any("heteroresistant" in m for m in summary["resistance_mechanism"])

    positions = summary["linezolid_23s_positions"]
    assert len(positions) == 1
    g2576 = positions[0]
    assert g2576["ecoli_position"] == 2576
    assert g2576["mutation"] == "G2576T"
    assert g2576["zygosity"] == "heteroresistant"
    # 6 rrn operons in E. faecium; a third of them carries the allele.
    assert g2576["rrn_copy_number"] == 6
    assert g2576["est_mutated_operons"] == 2
    assert g2576["af_ci_low"] < g2576["max_resistance_af"] < g2576["af_ci_high"]

    text = (out / "HET001.linezolid_amr.txt").read_text()
    assert "HETERORESISTANCE" in text
    assert "G2576T" in text

    with (out / "HET001.summary_wide.csv").open() as fh:
        wide = next(csv.DictReader(fh))
    assert wide["linezolid_call"] == "POS"
    assert float(wide["G2576T"]) > 0.15

    evidence = out / "rrna23s" / "HET001.23S_lzd_evidence.tsv"
    with evidence.open() as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if r["ecoli_position"] == "2576" and r["alt_base"] == "T"]
    assert rows and rows[0]["filters"] == "PASS"
    assert int(rows[0]["alt_fwd"]) > 0 and int(rows[0]["alt_rev"]) > 0


@requires_tools
def test_fixed_resistance_is_not_labelled_heteroresistant(tmp_path):
    """An allele in every operon is fixed resistance, a different clinical picture."""
    r1, r2 = _make_reads(tmp_path, "FIX001", af=1.0, seed=7)
    out = tmp_path / "out"
    _run_cli(["run", "-1", str(r1), "-2", str(r2), "-O", ORGANISM, "-o", str(out), "-t", "2"])

    summary = json.loads((out / "FIX001.linezolid_amr.json").read_text())["summary"]
    assert summary["linezolid_resistance_call"] is True
    assert summary["heteroresistance_detected"] is False
    assert summary["linezolid_23s_positions"][0]["zygosity"] == "fixed"
    assert any("fixed" in m for m in summary["resistance_mechanism"])


@requires_tools
def test_wildtype_sample_is_negative_and_quiet(tmp_path):
    """A clean sample must be negative and must not report spurious rejections."""
    r1, r2 = _make_reads(tmp_path, "WT001", af=0.0, seed=3)
    out = tmp_path / "out"
    result = _run_cli(["run", "-1", str(r1), "-2", str(r2), "-O", ORGANISM, "-o", str(out), "-t", "2"])
    assert "negative" in result.output

    summary = json.loads((out / "WT001.linezolid_amr.json").read_text())["summary"]
    assert summary["linezolid_resistance_call"] is False
    assert summary["n_23s_positions_with_resistance_allele"] == 0
    # Sub-threshold sequencing noise must not be dressed up as a rejected finding.
    assert summary["n_23s_positions_rejected_by_filters"] == 0
    assert "rejected by filters" not in result.output


@requires_tools
def test_run_parameters_are_recorded_for_reproducibility(tmp_path):
    """Every threshold that shaped the call must be in the report."""
    r1, r2 = _make_reads(tmp_path, "PARAM01", af=0.4, seed=11)
    out = tmp_path / "out"
    _run_cli([
        "run", "-1", str(r1), "-2", str(r2), "-O", ORGANISM, "-o", str(out),
        "-t", "2", "--min-af", "0.20", "--min-mapq", "30", "--evidence-tier", "established",
    ])
    params = json.loads((out / "PARAM01.linezolid_amr.json").read_text())["parameters"]
    assert params["min_af"] == 0.20
    assert params["min_mapq"] == 30
    assert params["evidence_tier"] == "established"
    assert params["reads_only"] is True


@requires_tools
def test_folder_mode_reads_only_cohort(tmp_path):
    """Cohort mode must aggregate every sample into the shared summary."""
    indir = tmp_path / "in"
    indir.mkdir()
    _make_reads(indir, "S1", af=0.0, seed=1)
    _make_reads(indir, "S2", af=0.45, seed=2)
    _make_reads(indir, "S3", af=1.0, seed=3)

    out = tmp_path / "out"
    result = _run_cli([
        "folder", "-i", str(indir), "-o", str(out),
        "--reads-only", "-O", ORGANISM, "-t", "2",
    ])
    assert "Linezolid-positive samples: 2 / 3" in result.output

    with (out / "ALL_samples.summary_wide.csv").open() as fh:
        rows = {r["sample"]: r for r in csv.DictReader(fh)}
    assert set(rows) == {"S1", "S2", "S3"}
    assert rows["S1"]["linezolid_call"] == "neg"
    assert rows["S2"]["linezolid_call"] == "POS"
    assert rows["S3"]["linezolid_call"] == "POS"
    assert not (out / "ALL_samples.failed.tsv").exists()


@requires_tools
def test_evidence_tier_gate_changes_the_call(tmp_path):
    """--evidence-tier must actually move the bar, not just annotate.

    A2572U is an engineered substitution: reported, but not clinical evidence.
    """
    _name, seq = read_fasta_single(ref_mod.organism_fasta_path(ORGANISM))
    r1 = tmp_path / "TIER_R1.fastq.gz"
    r2 = tmp_path / "TIER_R2.fastq.gz"
    # Species coordinate of E. coli 2572 for E. faecium is 2586.
    simulate_reads(seq, r1, r2, variants={2586: ("T", 0.6)}, depth=250, seed=21)

    default_out = tmp_path / "default"
    _run_cli(["run", "-1", str(r1), "-2", str(r2), "-O", ORGANISM, "-o", str(default_out), "-t", "2"])
    default_summary = json.loads(
        (default_out / "TIER.linezolid_amr.json").read_text()
    )["summary"]
    assert default_summary["linezolid_resistance_call"] is False, (
        "an experimental-tier position drove a positive call at the default tier"
    )

    wide_out = tmp_path / "wide"
    _run_cli([
        "run", "-1", str(r1), "-2", str(r2), "-O", ORGANISM, "-o", str(wide_out),
        "-t", "2", "--evidence-tier", "experimental",
    ])
    wide_summary = json.loads((wide_out / "TIER.linezolid_amr.json").read_text())["summary"]
    assert wide_summary["linezolid_resistance_call"] is True


@requires_tools
def test_no_reads_and_no_assembly_is_rejected():
    result = CliRunner().invoke(main, ["run", "-o", "/tmp/nowhere"])
    assert result.exit_code != 0
    assert "nothing to analyse" in result.output.lower()
