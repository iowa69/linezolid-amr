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
