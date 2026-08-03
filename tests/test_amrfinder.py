"""Unit tests for amrfinder.py parsing logic."""

from __future__ import annotations

from pathlib import Path

from linezolid_amr import amrfinder as amr_mod


AMRFINDER_TSV = (
    "Name\tContig id\tStart\tStop\tStrand\tGene symbol\tSequence name\tScope\t"
    "Element type\tElement subtype\tClass\tSubclass\tMethod\tTarget length\t"
    "Reference sequence length\t% Coverage of reference sequence\t"
    "% Identity to reference sequence\tAlignment length\tAccession of closest sequence\t"
    "Name of closest sequence\tHMM id\tHMM description\n"
    "S1\tctg1\t1000\t2500\t+\tcfr\tphenicol-lincosamide-oxazolidinone-pleuromutilin-streptogramin A 23S rRNA methyltransferase Cfr\t"
    "core\tAMR\tAMR\tPHENICOL/LINCOSAMIDE/OXAZOLIDINONE/PLEUROMUTILIN/STREPTOGRAMIN\tCFR\tEXACTX\t"
    "1149\t1149\t100.00\t100.00\t1149\tNG_047815.1\tcfr\tNA\tNA\n"
    "S1\tctg2\t100\t1100\t+\toptrA\tABC-F type ribosomal protection protein OptrA\t"
    "core\tAMR\tAMR\tOXAZOLIDINONE/PHENICOL\tOPTRA\tEXACTX\t"
    "1971\t1971\t100.00\t100.00\t1971\tNG_048023.1\toptrA\tNA\tNA\n"
    "S1\tctg3\t500\t900\t+\tmecA\tPBP2a\t"
    "core\tAMR\tAMR\tBETA-LACTAM\tMETHICILLIN\tEXACTX\t"
    "2007\t2007\t100.00\t100.00\t2007\tNG_047902.1\tmecA\tNA\tNA\n"
)


def test_parse_amrfinder_tsv(tmp_path: Path):
    tsv = tmp_path / "amr.tsv"
    tsv.write_text(AMRFINDER_TSV)
    hits = amr_mod.parse_amrfinder_tsv(tsv)
    assert len(hits) == 3
    assert hits[0].gene_symbol == "cfr"
    assert hits[0].coverage_pct == 100.0


def test_linezolid_flagging(tmp_path: Path):
    tsv = tmp_path / "amr.tsv"
    tsv.write_text(AMRFINDER_TSV)
    hits = amr_mod.parse_amrfinder_tsv(tsv)
    lzd = amr_mod.linezolid_relevant_hits(hits)
    syms = {h.gene_symbol for h in lzd}
    assert "cfr" in syms
    assert "optrA" in syms
    assert "mecA" not in syms


def test_empty_tsv(tmp_path: Path):
    tsv = tmp_path / "missing.tsv"
    assert amr_mod.parse_amrfinder_tsv(tsv) == []


# --------------------------------------------------------------------------
# Database discovery
# --------------------------------------------------------------------------

def test_db_dir_is_read_from_amrfinder_not_guessed(monkeypatch):
    """The database path must come from `amrfinder -V`, not the binary's prefix.

    In a conda layout the database commonly sits under the base prefix while
    the binary lives in an environment. Guessing made every run believe the
    database was missing and re-trigger a ~150 MB update — once per sample in
    folder mode.
    """
    import subprocess as sp

    amr_mod._db_dir_cache = None
    monkeypatch.setattr(amr_mod.shutil, "which", lambda _n: "/opt/envs/amr/bin/amrfinder")

    def fake_run(cmd, **kw):
        assert cmd == ["amrfinder", "-V"]
        return sp.CompletedProcess(
            cmd, 0,
            stdout=(
                "Software directory: '/opt/envs/amr/bin/'\n"
                "Software version: 4.2.7\n"
                "Database directory: '/opt/base/share/amrfinderplus/data/2026-05-15.1'\n"
                "Database version: 2026-05-15.1\n"
            ),
            stderr="",
        )

    monkeypatch.setattr(amr_mod.subprocess, "run", fake_run)
    try:
        assert str(amr_mod._amrfinder_db_dir()) == "/opt/base/share/amrfinderplus/data/2026-05-15.1"
    finally:
        amr_mod._db_dir_cache = None


def test_db_dir_falls_back_when_version_output_is_unparsable(monkeypatch):
    import subprocess as sp

    amr_mod._db_dir_cache = None
    monkeypatch.setattr(amr_mod.shutil, "which", lambda _n: "/opt/envs/amr/bin/amrfinder")
    monkeypatch.setattr(
        amr_mod.subprocess, "run",
        lambda cmd, **kw: sp.CompletedProcess(cmd, 1, stdout="", stderr="boom"),
    )
    try:
        got = amr_mod._amrfinder_db_dir()
        assert got is not None and got.name == "latest"
    finally:
        amr_mod._db_dir_cache = None


def test_db_dir_is_none_without_amrfinder(monkeypatch):
    amr_mod._db_dir_cache = None
    monkeypatch.setattr(amr_mod.shutil, "which", lambda _n: None)
    try:
        assert amr_mod._amrfinder_db_dir() is None
        assert amr_mod.amrfinder_db_ready() is False
    finally:
        amr_mod._db_dir_cache = None


def test_failure_message_quotes_the_log(tmp_path):
    """A bare 'see the log' is useless in folder mode; the cause must be inline."""
    log = tmp_path / "amrfinder.log"
    log.write_text(
        "Running: amrfinder ...\n"
        "Database version: 2026-03-24.1\n"
        "*** ERROR ***\n"
        "Database requires software version at least 4.2.0\n"
    )
    tail = amr_mod._log_tail(log)
    assert "Database requires software version at least 4.2.0" in tail


def test_log_tail_handles_missing_file(tmp_path):
    assert "unavailable" in amr_mod._log_tail(tmp_path / "nope.log")
