"""Unit tests for folder-mode sample discovery."""

from __future__ import annotations

from pathlib import Path

from linezolid_amr.cli import discover_samples


def _touch(p: Path):
    p.write_bytes(b"")


def test_folder_discovery_illumina_full(tmp_path: Path):
    _touch(tmp_path / "S1.fasta")
    _touch(tmp_path / "S1_R1_001.fastq.gz")
    _touch(tmp_path / "S1_R2_001.fastq.gz")
    res = discover_samples(tmp_path)
    assert len(res) == 1
    assert res[0]["sample"] == "S1"
    assert res[0]["r1"].name == "S1_R1_001.fastq.gz"
    assert res[0]["r2"].name == "S1_R2_001.fastq.gz"


def test_folder_discovery_short_R(tmp_path: Path):
    _touch(tmp_path / "iso.fa")
    _touch(tmp_path / "iso_R1.fq.gz")
    _touch(tmp_path / "iso_R2.fq.gz")
    res = discover_samples(tmp_path)
    assert len(res) == 1
    assert res[0]["sample"] == "iso"
    assert res[0]["r1"].suffix == ".gz"


def test_folder_discovery_numeric(tmp_path: Path):
    _touch(tmp_path / "x.fasta")
    _touch(tmp_path / "x_1.fastq.gz")
    _touch(tmp_path / "x_2.fastq.gz")
    res = discover_samples(tmp_path)
    assert len(res) == 1
    assert res[0]["r1"].name == "x_1.fastq.gz"
    assert res[0]["r2"].name == "x_2.fastq.gz"


def test_folder_discovery_missing_reads_flagged(tmp_path: Path):
    _touch(tmp_path / "lonely.fasta")
    res = discover_samples(tmp_path)
    assert len(res) == 1
    assert res[0]["missing_reads"] is True


def test_folder_discovery_multi_sample(tmp_path: Path):
    for name in ("A", "B", "C"):
        _touch(tmp_path / f"{name}.fasta")
        _touch(tmp_path / f"{name}_R1.fastq.gz")
        _touch(tmp_path / f"{name}_R2.fastq.gz")
    res = discover_samples(tmp_path)
    names = {r["sample"] for r in res}
    assert names == {"A", "B", "C"}
